!===================================================================================================
! Module Name: module_stream_temperature_eq
!
! Abstract:
!   Simple equilibrium stream temperature model for use with the UDMP reach-based routing
!   (UDMP_OPT=1, water_temperature=2).  Implements a WEAP-style exponential relaxation
!   of stream temperature toward a local equilibrium temperature driven by air temperature
!   and solar radiation.
!
! Physics:
!   Equilibrium temperature:
!     T_e = T_air_smooth + SW_net / K_sw
!     T_air_smooth = exponential moving average of the reach air temperature
!                    (e-folding time Route_Link TauTair, default T_smooth_tau below)
!     SW_net = SW_in * (1 - albedo) * (1 - shade_factor)
!     K_sw   = shortwave equilibrium coefficient [W/m2/C]   (Route_Link Ksw, default K_sw below)
!     K_e    = bulk surface heat exchange coefficient [W/m2/C]  (Route_Link Kexch, default K_exchange below)
!
!   Temperature routing (exponential relaxation along each reach):
!     T_out = T_e + (T_in - T_e) * exp(-K_rate * travel_time)
!     K_rate       = K_e / (rho_w * c_p * depth)              [1/s]
!     travel_time  = reach_length / velocity                   [s]
!
!   Inflow temperature (volume-weighted mixing):
!     T_in = [sum(Q_j * T_j, upstream j) + T_land*Q_land + T_bucket*Q_gw] /
!            [sum(Q_j, upstream j)        + Q_land + Q_gw]
!     (T_land per surf_temp_opt: Noah-MP runoff temperature at the cells where it enters
!      the reach, or instantaneous T_air; floored at 0 °C. T_bucket the groundwater
!      bucket outflow temperature tracked in simp_gw_buck_nhd, falling back to the fixed
!      T_gw parameter below where no bucket temperature is available)
!
!   All water temperatures are clamped to [0, 40] °C.  There is no river ice
!   physics, so water is held at freezing when T_e or inflows would go below 0 °C.
!     Headwater (no upstream inflow): T_in = T_e
!
!   Lake outlets:
!     Reaches marked TYPEL=1 (reservoir pour points) receive their temperature from the
!     1-D CLM lake model: T = LakeOutflowTemp(k) - 273.15 [°C].
!
! Upstream inflow / MPI (pull pattern, mirrors drive_CHANNEL_RSL):
!   In UDMP mode TO_NODE holds the downstream ComID (read from Route_Link `to`),
!   NOT a local array index, so a push via TO_NODE is impossible.  Instead each
!   reach PULLS its upstream contributions exactly the way the flow routing does:
!
!     1. gbcastReal2 gathers the local per-reach outflow QLINK(:,1) and the
!        previous-timestep temperature T_QSUM_prev into the compacted global
!        to-node index space (size rt_domain%nToInd) that gtoNode points into.
!     2. For each local reach i (1..NLINKSL), the upstream links are
!        m = gtoNode(i, n+1), n = 1..gtoNode(i,1); their Q and Q*T are summed
!        into Q_acc(i)/QT_acc(i).
!     3. Upstream temperatures therefore carry a 1-timestep lag (T_QSUM_prev),
!        which is rank-count independent and negligible at hourly dt.
!
!   Serial mode uses the same 1-step lag with the O(n^2) LINKID==TO_NODE match
!   that drive_CHANNEL_RSL uses in its non-MPP branch.
!
! Coupling:
!   - Reach-scale met forcing (T2, downward SW) is the catchment area-weighted mean
!     from the UDMP spatialweights mapping (module_UDMAP::UDMP_grid2reach_mean),
!     the same mapping used to aggregate land runoff to reaches.
!   - Reads from and writes to rt_domain(did)%T_QSUM, which is already output to CHRTOUT.
!   - Enabled via hydro.namelist: water_temperature = 2
!
! References:
!   Edinger et al. (1974), J. Water Pollut. Control Fed.
!   WEAP Technical Reference, Stockholm Environment Institute.
!   Dingman (2002), Physical Hydrology, Prentice Hall.
!===================================================================================================

module module_stream_temperature_eq

   use module_lsm_forcing,        only: weather2d
   use config_base,               only: nlst
   use module_RT_data,            only: rt_domain
   use module_UDMAP,              only: UDMP_grid2reach_mean, UDMP_lateral_temp
   use module_hydro_stop,         only: HYDRO_stop
#ifdef MPP_LAND
   use module_mpp_reachls,        only: gbcastReal2
   use module_mpp_land,           only: left_id, down_id, sum_double, sum_real1d, my_id, io_id
   use module_mpp_reachls,        only: ReachLS_decomp
#endif
   use netcdf

   implicit none
   private
   public :: run_stream_temp_eq, init_stream_temp_params

   ! --- Physical constants ---
   real, parameter :: rho_w      = 1000.0   ! water density [kg/m3]
   real, parameter :: c_p_w      = 4182.0   ! water specific heat [J/kg/°C]

   ! --- Default model parameters ---
   ! K_exchange, K_sw and T_smooth_tau (and the shade fraction and surface
   ! albedo, default 0) are per-reach parameters read from optional Route_Link
   ! variables Kexch, Ksw, TauTair, Shade and Albedo (init_stream_temp_params);
   ! the values below are the defaults where a variable is absent or invalid.
   real, parameter :: K_exchange =  50.0    ! bulk surface heat exchange coeff [W/m2/°C]
                                            ! Controls relaxation rate: K_rate = K_exchange/(rho_w*c_p*depth)
                                            ! Bracketed against the USGS gauges: 150 relaxed the Cannonsville
                                            ! release to equilibrium within 10 km (tailwaters +8 °C), while 25
                                            ! held the cold anomaly too long and left the far-field mainstem
                                            ! 3-7 °C too cold.  No single K satisfies both ends, since
                                            ! K_rate scales as 1/depth and the modeled West Branch depth
                                            ! (0.22-0.32 m) is biased shallow.
   real, parameter :: K_sw       = 400.0    ! SW equilibrium coefficient [W/m2/°C]
   integer, parameter :: surf_temp_opt = 2  ! temperature of non-bucket lateral inflow (QLateral - bucket):
                                            !  1 = reach instantaneous air temperature
                                            !  2 = Noah-MP infxsrt_t (0 C under snow, else ground temp TG)
                                            !      at the cells where the water enters the reach, weighted
                                            !      by delivered volume (falls back to 1 where unavailable)
                                            ! Controls solar contribution to T_e = T_air + SW_net/K_sw
                                            ! Higher value -> smaller diurnal solar swing in T_e
   real, parameter :: T_gw       =  12.0    ! fallback groundwater inflow temperature [°C],
                                            ! used only where no bucket temperature exists
                                            ! ~mean annual air temperature for DRB; buffers diurnal swing

   ! --- Physical bounds ---
   real, parameter :: T_min      = -15.0    ! lower clamp on AIR temperature [°C]
   real, parameter :: T_max      =  40.0    ! upper clamp on air and water temperature [°C]
   real, parameter :: T_water_min =  0.0    ! lower clamp on liquid WATER temperature [°C]
                                            ! (no ice physics: water is held at freezing)
   real, parameter :: min_depth  =  0.02    ! minimum depth for relaxation [m]
   real, parameter :: min_flow   =  1.0e-4  ! minimum flow for mixing    [m3/s]

   ! --- 24-hour EMA smoothing of T_air for equilibrium temperature ---
   ! Removes the sub-daily diurnal cycle from T_e while preserving synoptic
   ! weather variations.  EMA_alpha = DT / tau is computed at runtime from the
   ! namelist timestep so the e-folding time stays 24 hours at any DT.
   ! State lives in rt_domain(did)%T_air_smooth (per-domain, restart-persisted);
   ! entries still at the -9999 sentinel are seeded from current forcing.
   real, parameter :: T_smooth_tau = 86400.0              ! e-folding time [s]

contains

   !---------------------------------------------------------------------------
   ! Subroutine: run_stream_temp_eq
   !
   ! Computes equilibrium stream temperature for every reach in domain `did`
   ! and stores results in rt_domain(did)%T_QSUM [°C].
   !---------------------------------------------------------------------------
   subroutine run_stream_temp_eq(did)
      integer, intent(in) :: did

      integer :: i, k, n, m
      integer :: nlinksl, nlakes, ixrt, jxrt

      real    :: T_air_i, SW_net_i, T_e_i, T_in_i, T_out_i
      real    :: T_gw_i                       ! groundwater inflow temperature at this reach [°C]
      real    :: T_land_i                     ! non-bucket lateral inflow temperature at this reach [°C]
      real,    allocatable :: T_lat(:)        ! non-bucket lateral inflow temperature per reach [K]
      real    :: K_rate, travel_time, depth_i, vel_i, Q_this
      real    :: Q_total_i, QT_total_i, albedo_i, shade_i
      real    :: Tair_fallback, SW_fallback   ! domain-mean forcing for unmapped reaches
      real    :: T_air_raw_i                  ! per-reach instantaneous T_air before smoothing [°C]
      real    :: EMA_alpha                    ! EMA weight = DT / T_smooth_tau
      real*8  :: sumT, sumSW, ncnt            ! global sums for the fallback mean
      real,    allocatable :: T2_rt(:,:), SW_rt(:,:)  ! forcing on the routing grid

      real,    allocatable :: Q_acc(:)       ! upstream flow accumulated at each reach [m3/s]
      real,    allocatable :: QT_acc(:)      ! upstream Q*T accumulated at each reach  [m3/s * °C]
      real,    allocatable :: lk_Q(:)        ! lake inflow from upstream reaches [m3/s]
      real,    allocatable :: lk_QT(:)       ! lake inflow Q*T from upstream reaches [m3/s * °C]
#ifdef MPP_LAND
      real,    allocatable :: gQup(:)        ! upstream outflow in global to-node index space [m3/s]
      real,    allocatable :: gTup(:)        ! previous-step temperature, same index space    [°C]
#endif

      nlinksl  = rt_domain(did)%NLINKSL      ! local reach count (UDMP link space)
      nlakes   = rt_domain(did)%NLAKES

      if (nlst(did)%UDMP_OPT .ne. 1) then
         call hydro_stop("run_stream_temp_eq: water_temperature=2 requires UDMP_OPT=1")
      end if

      ! -----------------------------------------------------------------------
      ! Step 1: Reach-scale met forcing.
      !   Air temperature [K] and downward shortwave [W/m2] for each reach are
      !   the area-weighted means over the reach's catchment, using the same
      !   spatialweights (UDMP) mapping that aggregates land runoff to reaches.
      !   The LSM-grid fields are first copied onto the routing grid with the
      !   same index/halo convention as RunOffDisag, so they line up with the
      !   mapping's cell_i/cell_j.  Results go into Tair_fine_1D and
      !   SolarRad_fine_1D (1:nlinksl); reaches with no catchment polygon are
      !   set to 0 so the validity checks below route them to the fallback.
      !
      !   (Previously this used runWeatherDisagg with CHANXI/CHANYJ, which in
      !   UDMP mode index channel *pixels* in raster order, not reaches, so
      !   each reach received forcing from an unrelated pixel in its rank's
      !   subdomain.)
      ! -----------------------------------------------------------------------
      ixrt = rt_domain(did)%IXRT
      jxrt = rt_domain(did)%JXRT
      allocate(T2_rt(ixrt,jxrt), SW_rt(ixrt,jxrt))
      call lsm_to_routing_grid(weather2d%T2,    T2_rt, nlst(did)%AGGFACTRT, &
                               rt_domain(did)%IX, rt_domain(did)%JX)
      call lsm_to_routing_grid(weather2d%short, SW_rt, nlst(did)%AGGFACTRT, &
                               rt_domain(did)%IX, rt_domain(did)%JX)
      ! SW can legitimately be 0 (night); only negative/fill values are excluded,
      ! via the fill threshold.  Cells never written by the copy keep -9999.
      call UDMP_grid2reach_mean(T2_rt, rt_domain(did)%Tair_fine_1D(1:nlinksl),     &
                                fill_thresh=150.0, missing_val=0.0)
      call UDMP_grid2reach_mean(SW_rt, rt_domain(did)%SolarRad_fine_1D(1:nlinksl), &
                                fill_thresh=0.0,   missing_val=0.0)
      deallocate(T2_rt, SW_rt)

      ! Reaches with no catchment polygon (no mapped cells) take the forcing of
      ! their upstream reach(es); see fill_forcing_from_upstream.
      call fill_forcing_from_upstream(did, nlinksl)

      ! -----------------------------------------------------------------------
      ! Step 1b: Domain-mean Tair and SW as last-resort fallback for reaches
      !          with no catchment and no upstream reach with forcing (e.g.
      !          zero-flow headwater stubs; Tair_fine_1D still 0).  Summed
      !          across all ranks so every rank uses the same global mean.
      ! -----------------------------------------------------------------------
      sumT  = dble(sum(weather2d%T2))
      sumSW = dble(sum(max(0.0, weather2d%short)))
      ncnt  = dble(size(weather2d%T2))
#ifdef MPP_LAND
      call sum_double(sumT)
      call sum_double(sumSW)
      call sum_double(ncnt)
#endif
      Tair_fallback = real(sumT  / ncnt)
      SW_fallback   = real(sumSW / ncnt)

      ! -----------------------------------------------------------------------
      ! Step 1c: Update 24-hour EMA of T_air per reach.
      !   EMA_alpha is derived from the namelist timestep so the 24-h e-folding
      !   time is independent of DT.  Reaches still at the -9999 allocation
      !   sentinel (cold start, or restart file without T_air_smooth) are
      !   seeded from current forcing; reaches restored from a restart file
      !   continue their EMA seamlessly.
      ! -----------------------------------------------------------------------

      do i = 1, nlinksl
         if (rt_domain(did)%Tair_fine_1D(i) .gt. 200.0 .and. &
             rt_domain(did)%Tair_fine_1D(i) .lt. 350.0) then
            T_air_raw_i = rt_domain(did)%Tair_fine_1D(i) - 273.15
         else
            T_air_raw_i = Tair_fallback - 273.15
         end if
         T_air_raw_i = min(T_max, max(T_min, T_air_raw_i))
         EMA_alpha = min(1.0, nlst(did)%DT / rt_domain(did)%T_TauTair(i))
         if (rt_domain(did)%T_air_smooth(i) .lt. -100.0) then
            rt_domain(did)%T_air_smooth(i) = T_air_raw_i
         else
            rt_domain(did)%T_air_smooth(i) = EMA_alpha * T_air_raw_i &
                                           + (1.0 - EMA_alpha) * rt_domain(did)%T_air_smooth(i)
         end if
      end do

      ! -----------------------------------------------------------------------
      ! Upstream inflow pull (see module header).  For each local reach,
      ! accumulate upstream outflow Q and Q*T_prev the same way the flow
      ! routing does: via gtoNode into the globally gathered arrays (MPP), or
      ! via the LINKID==TO_NODE ComID match (serial).  TO_NODE holds ComIDs in
      ! UDMP mode and must never be used as an array index.
      ! -----------------------------------------------------------------------
      ! Non-bucket lateral inflow temperature per reach [K] (surf_temp_opt=2).
      ! Collective in MPP mode: every rank calls it (or none do).
      allocate(T_lat(max(1,nlinksl)), source=-9999.0)
      if (surf_temp_opt .eq. 2 .and. allocated(rt_domain(did)%INFXSRT_T_RT) .and. &
          allocated(rt_domain(did)%QSTRMVOLRT_TS) .and. allocated(rt_domain(did)%landRunOff)) then
         call UDMP_lateral_temp(rt_domain(did)%INFXSRT_T_RT, rt_domain(did)%QSTRMVOLRT_TS,   &
                                rt_domain(did)%landRunOff, nlst(did)%OVRTSWCRT,             &
                                nlst(did)%SUBRTSWCRT, 200.0, 350.0, T_lat(1:nlinksl), -9999.0)
      end if

      allocate(Q_acc(max(1,nlinksl)),  source=0.0)
      allocate(QT_acc(max(1,nlinksl)), source=0.0)

#ifdef MPP_LAND
      allocate(gQup(max(1, rt_domain(did)%nToInd)), source=0.0)
      allocate(gTup(max(1, rt_domain(did)%nToInd)), source=0.0)
      call gbcastReal2(rt_domain(did)%toNodeInd, rt_domain(did)%nToInd, &
                       rt_domain(did)%QLINK(1:nlinksl,1), nlinksl, gQup)
      call gbcastReal2(rt_domain(did)%toNodeInd, rt_domain(did)%nToInd, &
                       rt_domain(did)%T_QSUM_prev(1:nlinksl), nlinksl, gTup)
      do i = 1, nlinksl
         do n = 1, int(rt_domain(did)%gtoNode(i,1))
            m = int(rt_domain(did)%gtoNode(i,n+1))
            if (gQup(m) .gt. 0.0) then
               Q_acc(i)  = Q_acc(i)  + gQup(m)
               QT_acc(i) = QT_acc(i) + gQup(m) * gTup(m)
            end if
         end do
      end do
      deallocate(gQup, gTup)
#else
      do m = 1, nlinksl
         Q_this = max(0.0, rt_domain(did)%QLINK(m,1))
         if (Q_this .le. 0.0) cycle
         do i = 1, nlinksl
            if (rt_domain(did)%LINKID(i) .eq. rt_domain(did)%TO_NODE(m)) then
               Q_acc(i)  = Q_acc(i)  + Q_this
               QT_acc(i) = QT_acc(i) + Q_this * rt_domain(did)%T_QSUM_prev(m)
            end if
         end do
      end do
#endif

      ! -----------------------------------------------------------------------
      ! Step 2: Assign lake/reservoir temperatures from the 1-D CLM lake model.
      !         TYPEL=1 marks the lake pour-point reach; LAKEIDX maps that reach
      !         to the lake index in LakeOutflowTemp (units: K, top-layer).
      !
      !         Fallback: if LakeOutflowTemp is NaN or out of the physical range
      !         200–400 K (lake model not active / not yet spun up), compute the
      !         local equilibrium temperature at the pour-point reach location
      !         using the same T_e formula as for regular reaches.  This avoids
      !         leaving the outlet stuck at the 10 °C initialisation value.
      ! -----------------------------------------------------------------------
      if (nlakes .gt. 0) then
         ! Lake inflow temperature for the 1-D lake model's inflow heat: the
         ! upstream inflow pulled into each lake pour-point reach above is the
         ! lake inflow (flow-weighted, previous-step reach temperatures, the
         ! same lag as for regular reaches).  Summed over pour-point reaches
         ! and ranks so every rank holds the whole-lake value.
         allocate(lk_Q(nlakes), source=0.0)
         allocate(lk_QT(nlakes), source=0.0)
         do i = 1, nlinksl
            if (rt_domain(did)%TYPEL(i) .eq. 1) then
               k = rt_domain(did)%LAKEIDX(i)
               if (k .ge. 1 .and. k .le. nlakes) then
                  lk_Q(k)  = lk_Q(k)  + Q_acc(i)
                  lk_QT(k) = lk_QT(k) + QT_acc(i)
               end if
            end if
         end do
#ifdef MPP_LAND
         call sum_real1d(lk_Q, nlakes)
         call sum_real1d(lk_QT, nlakes)
#endif
         ! Lakes with no inflow this step keep their previous inflow temperature
         ! (unused while the inflow is zero).
         where (lk_Q .gt. 0.0) rt_domain(did)%LakeInflowTemp = lk_QT / lk_Q + 273.15
         deallocate(lk_Q, lk_QT)

         do i = 1, nlinksl
            if (rt_domain(did)%TYPEL(i) .eq. 1) then
               k = rt_domain(did)%LAKEIDX(i)
               if (k .ge. 1 .and. k .le. nlakes) then
                  associate(lake_t => rt_domain(did)%LakeOutflowTemp(k))
                     if (lake_t .ge. 200.0 .and. lake_t .le. 400.0) then
                        ! CLM 1-D lake model has a valid temperature — use it
                        rt_domain(did)%T_QSUM(i) = lake_t - 273.15
                     else
                        ! Fallback: equilibrium temperature at the pour-point location
                        ! (24-h smoothed T_air, same T_e formula as regular reaches)
                        if (rt_domain(did)%Tair_fine_1D(i) .gt. 200.0 .and. &
                            rt_domain(did)%Tair_fine_1D(i) .lt. 350.0) then
                           SW_net_i = max(0.0, rt_domain(did)%SolarRad_fine_1D(i))
                        else
                           SW_net_i = SW_fallback
                        end if
                        albedo_i = min(1.0, max(0.0, rt_domain(did)%S_Albedo(i)))
                        shade_i  = min(1.0, max(0.0, rt_domain(did)%F_Shade(i)))
                        SW_net_i = SW_net_i * (1.0 - albedo_i) * (1.0 - shade_i)
                        T_e_i    = rt_domain(did)%T_air_smooth(i) + SW_net_i / rt_domain(did)%T_Ksw(i)
                        rt_domain(did)%T_QSUM(i) = min(T_max, max(T_water_min, T_e_i))
                     end if
                  end associate
               end if
            end if
         end do
      end if

      ! -----------------------------------------------------------------------
      ! Step 3: Compute each reach's temperature from its (pulled) upstream
      !   inflow, lateral inflows, and exponential relaxation toward the local
      !   equilibrium temperature.  Upstream temperatures carry a 1-timestep
      !   lag (module header), so a single unordered pass suffices — no
      !   topological ordering and no downstream pushes.
      !
      !   Lateral inflows: surface runoff at T_air, groundwater at the bucket
      !   outflow temperature (fixed T_gw fallback).
      !   A future improvement would take these from Noah-MP.
      ! -----------------------------------------------------------------------
         ! Diagnostics (output only): reset, then filled with the values used below.
         rt_domain(did)%T_eq_diag(1:nlinksl)   = -9999.0
         rt_domain(did)%T_up_diag(1:nlinksl)   = -9999.0
         rt_domain(did)%T_lat_diag(1:nlinksl)  = -9999.0
         rt_domain(did)%T_gwin_diag(1:nlinksl) = -9999.0

         do i = 1, nlinksl

            ! Flow-weighted upstream inflow temperature (for lake pour-point
            ! reaches this is the lake inflow temperature)
            if (Q_acc(i) .gt. 0.0) rt_domain(did)%T_up_diag(i) = QT_acc(i) / Q_acc(i)

            ! --- Lake reaches: temperature already set in Step 2 ---
            if (rt_domain(did)%TYPEL(i) .eq. 1) cycle

            ! --- Met forcing at this reach (K → °C) ---
            ! Use domain-mean fallback for reaches with invalid disaggregation
            ! (Tair_fine_1D outside 200-350 K: either 0 K unfilled or fill value)
            if (rt_domain(did)%Tair_fine_1D(i) .gt. 200.0 .and. &
                rt_domain(did)%Tair_fine_1D(i) .lt. 350.0) then
               T_air_i  = rt_domain(did)%Tair_fine_1D(i) - 273.15
               SW_net_i = max(0.0, rt_domain(did)%SolarRad_fine_1D(i))
            else
               T_air_i  = Tair_fallback - 273.15
               SW_net_i = SW_fallback
            end if
            T_air_i = min(T_max, max(T_min, T_air_i))

            albedo_i = min(1.0, max(0.0, rt_domain(did)%S_Albedo(i)))
            shade_i  = min(1.0, max(0.0, rt_domain(did)%F_Shade(i)))
            SW_net_i = SW_net_i * (1.0 - albedo_i) * (1.0 - shade_i)

            ! --- Equilibrium temperature (24-hr smoothed T_air suppresses diurnal) ---
            ! T_air_i (instantaneous) is retained below for lateral inflow temperature.
            T_e_i = rt_domain(did)%T_air_smooth(i) + SW_net_i / rt_domain(did)%T_Ksw(i)
            T_e_i = min(T_max, max(T_water_min, T_e_i))

            ! --- Inflow temperature (volume-weighted mixing) ---
            ! Q_acc(i) holds upstream channel inflow (pulled above).  Surface
            ! runoff enters at instantaneous T_air, bucket flow at the bucket
            ! outflow temperature (K -> °C), or T_gw where none is available.
            T_gw_i = T_gw
            if (allocated(rt_domain(did)%qout_gwsubbas_t)) then
               if (i .le. size(rt_domain(did)%qout_gwsubbas_t)) then
                  if (rt_domain(did)%qout_gwsubbas_t(i) .gt. 200.0 .and. &
                      rt_domain(did)%qout_gwsubbas_t(i) .lt. 350.0) &
                     T_gw_i = rt_domain(did)%qout_gwsubbas_t(i) - 273.15
               end if
            end if
            T_gw_i = min(T_max, max(T_water_min, T_gw_i))
            Q_total_i  = Q_acc(i)                                               &
                       + max(0.0, rt_domain(did)%InflowLand(i))                 &
                       + max(0.0, rt_domain(did)%InflowGw(i))
            ! Non-bucket lateral inflow: Noah-MP runoff temperature at the cells
            ! where it enters the reach (surf_temp_opt=2), else air temperature.
            T_land_i = T_air_i
            if (T_lat(i) .gt. 200.0) T_land_i = T_lat(i) - 273.15
            T_land_i = min(T_max, max(T_water_min, T_land_i))
            QT_total_i = QT_acc(i)                                              &
                       + T_land_i * max(0.0, rt_domain(did)%InflowLand(i))     &
                       + T_gw_i  * max(0.0, rt_domain(did)%InflowGw(i))

            if (Q_total_i .gt. min_flow) then
               T_in_i = QT_total_i / Q_total_i
            else
               T_in_i = T_e_i   ! headwater: initialise at equilibrium
            end if
            T_in_i = min(T_max, max(T_water_min, T_in_i))

            ! --- Exponential relaxation toward equilibrium along the reach ---
            depth_i = rt_domain(did)%HLINK(i)
            vel_i   = max(0.0, rt_domain(did)%velocity(i))

            if (depth_i .gt. min_depth .and. vel_i .gt. 0.0) then
               K_rate      = rt_domain(did)%T_Kexch(i) / (rho_w * c_p_w * depth_i)
               travel_time = rt_domain(did)%CHANLEN(i) / vel_i
               T_out_i = T_e_i + (T_in_i - T_e_i) * exp(-K_rate * travel_time)
            else
               T_out_i = T_in_i   ! very low flow: pass inflow temperature through
            end if
            T_out_i = min(T_max, max(T_water_min, T_out_i))

            rt_domain(did)%T_QSUM(i) = T_out_i

            rt_domain(did)%T_eq_diag(i)   = T_e_i
            rt_domain(did)%T_lat_diag(i)  = T_land_i
            rt_domain(did)%T_gwin_diag(i) = T_gw_i

         end do   ! i

      deallocate(Q_acc, QT_acc, T_lat)

      ! Save this timestep's temperatures: next call gathers them as the
      ! (1-step lagged) upstream inflow temperatures.
      rt_domain(did)%T_QSUM_prev = rt_domain(did)%T_QSUM

   end subroutine run_stream_temp_eq

   !---------------------------------------------------------------------------
   ! Subroutine: init_stream_temp_params
   !
   ! Per-reach stream temperature parameters from optional Route_Link
   ! variables (dimension = reaches, same order as 'link'):
   !   Kexch   bulk surface heat exchange coefficient [W/m2/C]   default K_exchange
   !   Ksw     shortwave equilibrium coefficient [W/m2/C]        default K_sw
   !   TauTair air temperature smoothing e-folding time [s]      default T_smooth_tau
   !   Shade   shade fraction of incoming shortwave [0-1]        default 0
   !   Albedo  water surface albedo [0-1]                        default 0
   ! A variable that is absent uses its default for every reach; values outside
   ! the valid range (Kexch, Ksw, TauTair <= 0; Shade, Albedo outside 0-1;
   ! NaN or fill values) use the default for that reach and are counted in the
   ! log.  Read on the IO rank and distributed like the other Route_Link
   ! fields.  Shade and Albedo fill F_Shade / S_Albedo, which the reach-based
   ! configuration otherwise leaves at 0 (the CHANPARM.TBL stream-order values
   ! are used only by the gridded channel).
   !---------------------------------------------------------------------------
   subroutine init_stream_temp_params(did)
      integer, intent(in) :: did

      integer, parameter :: npar = 5
      character(len=8), parameter :: pname(npar) = &
           [character(len=8) :: 'Kexch', 'Ksw', 'TauTair', 'Shade', 'Albedo']
      real    :: pdef(npar), plo(npar), phi(npar)
      logical :: lo_open(npar)                ! lower bound exclusive (> lo) or inclusive (>= lo)
      integer :: nlinksl, gnlinksl, k, j, ncid, varid, iret, nbad
      logical :: found
      real, allocatable :: gval(:), lval(:)

      pdef = [K_exchange, K_sw, T_smooth_tau, 0.0, 0.0]
      plo  = [0.0, 0.0, 0.0, 0.0, 0.0]
      phi  = [huge(1.0), huge(1.0), huge(1.0), 1.0, 1.0]
      lo_open = [.true., .true., .true., .false., .false.]

      nlinksl  = rt_domain(did)%NLINKSL
#ifdef MPP_LAND
      gnlinksl = rt_domain(did)%GNLINKSL
#else
      gnlinksl = nlinksl
#endif
      if (.not. allocated(rt_domain(did)%T_Kexch))   allocate(rt_domain(did)%T_Kexch(max(1,nlinksl)))
      if (.not. allocated(rt_domain(did)%T_Ksw))     allocate(rt_domain(did)%T_Ksw(max(1,nlinksl)))
      if (.not. allocated(rt_domain(did)%T_TauTair)) allocate(rt_domain(did)%T_TauTair(max(1,nlinksl)))
      ! Diagnostics for CHRTOUT (written whenever water_temperature = 2)
      if (.not. allocated(rt_domain(did)%T_eq_diag))   allocate(rt_domain(did)%T_eq_diag(max(1,nlinksl)))
      if (.not. allocated(rt_domain(did)%T_up_diag))   allocate(rt_domain(did)%T_up_diag(max(1,nlinksl)))
      if (.not. allocated(rt_domain(did)%T_lat_diag))  allocate(rt_domain(did)%T_lat_diag(max(1,nlinksl)))
      if (.not. allocated(rt_domain(did)%T_gwin_diag)) allocate(rt_domain(did)%T_gwin_diag(max(1,nlinksl)))
      rt_domain(did)%T_eq_diag   = -9999.0
      rt_domain(did)%T_up_diag   = -9999.0
      rt_domain(did)%T_lat_diag  = -9999.0
      rt_domain(did)%T_gwin_diag = -9999.0
      allocate(lval(max(1,nlinksl)))

      ncid = -1
#ifdef MPP_LAND
      if (my_id .eq. io_id) then
#endif
         allocate(gval(max(1,gnlinksl)))
         iret = nf90_open(trim(nlst(did)%route_link_f), nf90_nowrite, ncid)
         if (iret /= nf90_noerr) call hydro_stop('init_stream_temp_params: cannot open Route_Link file ' // &
                                                 trim(nlst(did)%route_link_f))
#ifdef MPP_LAND
      else
         allocate(gval(1))
      end if
#endif

      do k = 1, npar
#ifdef MPP_LAND
         if (my_id .eq. io_id) then
#endif
            gval = pdef(k)
            found = nf90_inq_varid(ncid, trim(pname(k)), varid) == nf90_noerr
            nbad = 0
            if (found) then
               iret = nf90_get_var(ncid, varid, gval(1:gnlinksl))
               if (iret /= nf90_noerr) call hydro_stop('init_stream_temp_params: cannot read Route_Link ' // &
                                                       trim(pname(k)))
               do j = 1, gnlinksl
                  if (gval(j) /= gval(j) .or. gval(j) > phi(k) .or. &
                      (lo_open(k) .and. gval(j) <= plo(k)) .or. (.not. lo_open(k) .and. gval(j) < plo(k))) then
                     gval(j) = pdef(k)
                     nbad = nbad + 1
                  end if
               end do
               write(6,'(a,a,a,i0,a,es11.4,a,es11.4,a,i0,a)') ' init_stream_temp_params: ', trim(pname(k)), &
                    ' read from Route_Link for ', gnlinksl, ' reaches (min ', minval(gval(1:gnlinksl)), &
                    ', max ', maxval(gval(1:gnlinksl)), '); ', nbad, ' invalid values replaced by the default'
            else
               write(6,'(a,a,a,es11.4)') ' init_stream_temp_params: ', trim(pname(k)), &
                    ' not in Route_Link; default for all reaches = ', pdef(k)
            end if
#ifdef MPP_LAND
         end if
         call ReachLS_decomp(gval, lval)
#else
         lval(1:nlinksl) = gval(1:nlinksl)
#endif
         select case (k)
         case (1); rt_domain(did)%T_Kexch(1:nlinksl)   = lval(1:nlinksl)
         case (2); rt_domain(did)%T_Ksw(1:nlinksl)     = lval(1:nlinksl)
         case (3); rt_domain(did)%T_TauTair(1:nlinksl) = lval(1:nlinksl)
         case (4); rt_domain(did)%F_Shade(1:nlinksl)   = lval(1:nlinksl)
         case (5); rt_domain(did)%S_Albedo(1:nlinksl)  = lval(1:nlinksl)
         end select
      end do

#ifdef MPP_LAND
      if (my_id .eq. io_id) then
#endif
         iret = nf90_close(ncid)
#ifdef MPP_LAND
      end if
#endif
      deallocate(gval, lval)
      call flush(6)

   end subroutine init_stream_temp_params

   !---------------------------------------------------------------------------
   ! fill_forcing_from_upstream: reaches with no catchment in the spatialweights
   !   mapping (Tair_fine_1D outside 200-350 K after the catchment average) take
   !   the air temperature and shortwave of their upstream reach(es): an average
   !   weighted by upstream outflow QLINK(:,1), or an equal-weight average when
   !   all upstream outflows are zero.  Repeated passes fill chains of reaches
   !   without catchments.  New values are staged and applied after each pass,
   !   so the result is independent of reach order and MPI rank count.  Reaches
   !   with no upstream reach carrying forcing are left unset (the caller then
   !   uses the domain-mean fallback).  Collective in MPP mode.
   !---------------------------------------------------------------------------
   subroutine fill_forcing_from_upstream(did, nlinksl)
      integer, intent(in) :: did, nlinksl
      integer, parameter  :: max_pass = 5
      real,    parameter  :: T_lo = 200.0, T_hi = 350.0
      real,    allocatable :: newT(:), newS(:)
#ifdef MPP_LAND
      real,    allocatable :: gQ(:), gT(:), gS(:)
#endif
      real*8  :: nmiss, nfill
      real    :: wsum, wT, wS, esum, eT, eS, w, Tu, Su
      integer :: pass, i, n, m

      allocate(newT(max(1,nlinksl)), newS(max(1,nlinksl)))
#ifdef MPP_LAND
      allocate(gQ(max(1,rt_domain(did)%nToInd)), gT(max(1,rt_domain(did)%nToInd)), &
               gS(max(1,rt_domain(did)%nToInd)))
#endif
      do pass = 1, max_pass
         nmiss = 0.d0
         do i = 1, nlinksl
            if (.not. valid(rt_domain(did)%Tair_fine_1D(i))) nmiss = nmiss + 1.d0
         end do
#ifdef MPP_LAND
         call sum_double(nmiss)
#endif
         if (nmiss .le. 0.d0) exit

#ifdef MPP_LAND
         if (pass .eq. 1) &
            call gbcastReal2(rt_domain(did)%toNodeInd, rt_domain(did)%nToInd, &
                             rt_domain(did)%QLINK(1:nlinksl,1), nlinksl, gQ)
         call gbcastReal2(rt_domain(did)%toNodeInd, rt_domain(did)%nToInd, &
                          rt_domain(did)%Tair_fine_1D(1:nlinksl), nlinksl, gT)
         call gbcastReal2(rt_domain(did)%toNodeInd, rt_domain(did)%nToInd, &
                          rt_domain(did)%SolarRad_fine_1D(1:nlinksl), nlinksl, gS)
#endif
         newT(1:max(1,nlinksl)) = -1.0
         nfill = 0.d0
         do i = 1, nlinksl
            if (valid(rt_domain(did)%Tair_fine_1D(i))) cycle
            wsum = 0.0; wT = 0.0; wS = 0.0; esum = 0.0; eT = 0.0; eS = 0.0
#ifdef MPP_LAND
            do n = 1, int(rt_domain(did)%gtoNode(i,1))
               m  = int(rt_domain(did)%gtoNode(i,n+1))
               Tu = gT(m); Su = gS(m); w = max(0.0, gQ(m))
#else
            do m = 1, nlinksl
               if (rt_domain(did)%TO_NODE(m) .ne. rt_domain(did)%LINKID(i)) cycle
               Tu = rt_domain(did)%Tair_fine_1D(m); Su = rt_domain(did)%SolarRad_fine_1D(m)
               w  = max(0.0, rt_domain(did)%QLINK(m,1))
#endif
               if (.not. valid(Tu)) cycle
               wsum = wsum + w;   wT = wT + w*Tu;   wS = wS + w*max(0.0,Su)
               esum = esum + 1.0; eT = eT + Tu;     eS = eS + max(0.0,Su)
            end do
            if (esum .gt. 0.0) then
               if (wsum .gt. 0.0) then
                  newT(i) = wT / wsum;  newS(i) = wS / wsum
               else
                  newT(i) = eT / esum;  newS(i) = eS / esum
               end if
               nfill = nfill + 1.d0
            end if
         end do
         do i = 1, nlinksl
            if (newT(i) .gt. 0.0) then
               rt_domain(did)%Tair_fine_1D(i)     = newT(i)
               rt_domain(did)%SolarRad_fine_1D(i) = newS(i)
            end if
         end do
#ifdef MPP_LAND
         call sum_double(nfill)
#endif
         if (nfill .le. 0.d0) exit
      end do

      deallocate(newT, newS)
#ifdef MPP_LAND
      deallocate(gQ, gT, gS)
#endif
   contains
      logical function valid(t)
         real, intent(in) :: t
         valid = (t .gt. T_lo .and. t .lt. T_hi)
      end function valid
   end subroutine fill_forcing_from_upstream

   !---------------------------------------------------------------------------
   ! lsm_to_routing_grid: copy an LSM-grid field onto the routing grid
   !   (nearest/replicate), using the same fine-grid index and MPI halo
   !   convention as RunOffDisag in module_HYDRO_drv, so that the result lines
   !   up with the UDMP mapping indices.  Cells not covered keep -9999.
   !---------------------------------------------------------------------------
   subroutine lsm_to_routing_grid(fld_lsm, fld_rt, aggfactrt, ix, jx)
      real,    intent(in)  :: fld_lsm(:,:)
      real,    intent(out) :: fld_rt(:,:)
      integer, intent(in)  :: aggfactrt, ix, jx
      integer :: i, j, ax, ay, ii, jj

      fld_rt = -9999.0
      do j = 1, jx
         do i = 1, ix
            do ay = aggfactrt-1, 0, -1
               do ax = aggfactrt-1, 0, -1
                  ii = i*aggfactrt - ax
                  jj = j*aggfactrt - ay
#ifdef MPP_LAND
                  if (left_id .ge. 0) ii = ii + 1
                  if (down_id .ge. 0) jj = jj + 1
#endif
                  if (ii .le. size(fld_rt,1) .and. jj .le. size(fld_rt,2)) &
                     fld_rt(ii,jj) = fld_lsm(i,j)
               end do
            end do
         end do
      end do
   end subroutine lsm_to_routing_grid

end module module_stream_temperature_eq
