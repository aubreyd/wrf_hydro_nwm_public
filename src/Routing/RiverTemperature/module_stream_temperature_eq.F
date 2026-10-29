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
!     T_e = T_air + SW_net / K_e
!     SW_net = SW_in * (1 - albedo) * (1 - shade_factor)
!     K_e    = bulk surface heat exchange coefficient [W/m2/C]  (parameter K_exchange below)
!
!   Temperature routing (exponential relaxation along each reach):
!     T_out = T_e + (T_in - T_e) * exp(-K_rate * travel_time)
!     K_rate       = K_e / (rho_w * c_p * depth)              [1/s]
!     travel_time  = reach_length / velocity                   [s]
!
!   Inflow temperature (volume-weighted mixing, push-accumulation pattern):
!     T_in = [sum(Q_j * T_j, upstream j) + T_air*(Q_land + Q_gw)] /
!            [sum(Q_j, upstream j)        + Q_land + Q_gw]
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
!   - Calls fineWeathData%runWeatherDisagg() to obtain reach-scale met forcing from the
!     existing disaggregation infrastructure (same as water_temperature=1 path).
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
   use config_base,               only: noah_lsm, nlst
   use module_RT_data,            only: rt_domain
   use disaggregateWeatherModule, only: fineWeathData
#ifdef MPP_LAND
   use module_mpp_reachls,        only: gbcastReal2
#endif

   implicit none
   private
   public :: run_stream_temp_eq

   ! --- Physical constants ---
   real, parameter :: rho_w      = 1000.0   ! water density [kg/m3]
   real, parameter :: c_p_w      = 4182.0   ! water specific heat [J/kg/°C]

   ! --- Default model parameters (can be promoted to namelist later) ---
   real, parameter :: K_exchange =  50.0    ! bulk surface heat exchange coeff [W/m2/°C]
                                            ! Controls relaxation rate: K_rate = K_exchange/(rho_w*c_p*depth)
                                            ! Bracketed against the USGS gauges: 150 relaxed the Cannonsville
                                            ! release to equilibrium within 10 km (tailwaters +8 °C), while 25
                                            ! held the cold anomaly too long and left the far-field mainstem
                                            ! 3-7 °C too cold.  No single K satisfies both ends, since
                                            ! K_rate scales as 1/depth and the modeled West Branch depth
                                            ! (0.22-0.32 m) is biased shallow.
   real, parameter :: K_sw       = 400.0    ! SW equilibrium coefficient [W/m2/°C]
                                            ! Controls solar contribution to T_e = T_air + SW_net/K_sw
                                            ! Higher value -> smaller diurnal solar swing in T_e
   real, parameter :: T_gw       =  12.0    ! groundwater inflow temperature [°C]
                                            ! ~mean annual air temperature for DRB; buffers diurnal swing

   ! --- Physical bounds ---
   real, parameter :: T_min      = -15.0    ! lower temperature clamp [°C]
   real, parameter :: T_max      =  40.0    ! upper temperature clamp [°C]
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
      integer :: nlinks, nlinksl, nlakes, nsoil

      real    :: T_air_i, SW_net_i, T_e_i, T_in_i, T_out_i
      real    :: K_rate, travel_time, depth_i, vel_i, Q_this
      real    :: Q_total_i, QT_total_i, albedo_i, shade_i
      real    :: Tair_fallback, SW_fallback   ! domain-mean forcing for unmapped reaches
      real    :: T_air_raw_i                  ! per-reach instantaneous T_air before smoothing [°C]
      real    :: EMA_alpha                    ! EMA weight = DT / T_smooth_tau

      real,    allocatable :: Q_acc(:)       ! upstream flow accumulated at each reach [m3/s]
      real,    allocatable :: QT_acc(:)      ! upstream Q*T accumulated at each reach  [m3/s * °C]
#ifdef MPP_LAND
      real,    allocatable :: gQup(:)        ! upstream outflow in global to-node index space [m3/s]
      real,    allocatable :: gTup(:)        ! previous-step temperature, same index space    [°C]
#endif

      nlinks   = rt_domain(did)%NLINKS
      nlinksl  = rt_domain(did)%NLINKSL      ! local reach count (UDMP link space)
      nlakes   = rt_domain(did)%NLAKES
      nsoil    = noah_lsm%nsoil

      ! -----------------------------------------------------------------------
      ! Step 1: Disaggregate gridded met forcing to reach scale.
      !         Populates Tair_fine_1D, SolarRad_fine_1D, etc. in rt_domain.
      !         TMN_temp/TG_temp/TSLB_temp are currently zero (not yet coupled
      !         from NoahMP state), so TSLB_1D_fine will also be zero; we use
      !         T_air as a first-order groundwater temperature proxy instead.
      ! -----------------------------------------------------------------------
      call fineWeathData%runWeatherDisagg(                                       &
               rt_domain(did)%IX,                                                &
               rt_domain(did)%JX,                                                &
               rt_domain(did)%IXRT,                                              &
               rt_domain(did)%JXRT,                                              &
               weather2d%T2,                                                     &
               weather2d%short,                                                  &
               weather2d%q2x,                                                    &
               weather2d%wind,                                                   &
               weather2d%sfcp,                                                   &
               nlinks,                                                           &
               rt_domain(did)%overland%streams_and_lakes%CH_NETRT,               &
               rt_domain(did)%CHANXI,                                            &
               rt_domain(did)%CHANYJ,                                            &
               rt_domain(did)%overland%streams_and_lakes%lake_mask,              &
               rt_domain(did)%TMN_temp,                                          &
               rt_domain(did)%TG_temp,                                           &
               rt_domain(did)%TSLB_temp,                                         &
               nsoil,                                                            &
               nlst(did)%AGGFACTRT,                                              &
               rt_domain(did)%Tair_fine_1D,                                      &
               rt_domain(did)%SolarRad_fine_1D,                                  &
               rt_domain(did)%RelH_fine_1D,                                      &
               rt_domain(did)%windSpd_fine_1D,                                   &
               rt_domain(did)%TMN_1D_fine,                                       &
               rt_domain(did)%TG_1D_fine,                                        &
               rt_domain(did)%TSLB_1D_fine,                                      &
               rt_domain(did)%SFCP_1D                                            &
               )

      ! -----------------------------------------------------------------------
      ! Step 1b: Compute domain-mean Tair and SW as fallback for reaches whose
      !          disaggregation mapping is invalid (Tair_fine_1D outside the
      !          physically plausible 200-350 K range: either 0 K unfilled or
      !          a large fill value from an unmapped grid cell).
      !          Use the gridded weather2d fields directly — always correctly
      !          populated from HRLDAS forcing regardless of routing grid mapping.
      ! -----------------------------------------------------------------------
      Tair_fallback = sum(weather2d%T2)                     / real(size(weather2d%T2))
      SW_fallback   = sum(max(0.0, weather2d%short))        / real(size(weather2d%short))

      ! -----------------------------------------------------------------------
      ! Step 1c: Update 24-hour EMA of T_air per reach.
      !   EMA_alpha is derived from the namelist timestep so the 24-h e-folding
      !   time is independent of DT.  Reaches still at the -9999 allocation
      !   sentinel (cold start, or restart file without T_air_smooth) are
      !   seeded from current forcing; reaches restored from a restart file
      !   continue their EMA seamlessly.
      ! -----------------------------------------------------------------------
      EMA_alpha = min(1.0, nlst(did)%DT / T_smooth_tau)

      do i = 1, nlinks
         if (rt_domain(did)%Tair_fine_1D(i) .gt. 200.0 .and. &
             rt_domain(did)%Tair_fine_1D(i) .lt. 350.0) then
            T_air_raw_i = rt_domain(did)%Tair_fine_1D(i) - 273.15
         else
            T_air_raw_i = Tair_fallback - 273.15
         end if
         T_air_raw_i = min(T_max, max(T_min, T_air_raw_i))
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
      allocate(Q_acc(nlinks),  source=0.0)
      allocate(QT_acc(nlinks), source=0.0)

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
         do i = 1, nlinks
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
                        T_e_i    = rt_domain(did)%T_air_smooth(i) + SW_net_i / K_sw
                        rt_domain(did)%T_QSUM(i) = min(T_max, max(T_min, T_e_i))
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
      !   Lateral inflows (surface runoff + groundwater) enter at T_air as a
      !   first-order approximation.  A future improvement would use soil
      !   temperature (TSLB) for groundwater once NoahMP state is coupled in.
      ! -----------------------------------------------------------------------
         do i = 1, nlinksl

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
            T_e_i = rt_domain(did)%T_air_smooth(i) + SW_net_i / K_sw
            T_e_i = min(T_max, max(T_min, T_e_i))

            ! --- Inflow temperature (volume-weighted mixing) ---
            ! Q_acc(i) holds upstream channel inflow (from pre-seed in MPP mode,
            ! or from the ORDER-loop push in serial mode).  Both lateral
            ! inflows enter at local air/groundwater temperature.
            Q_total_i  = Q_acc(i)                                               &
                       + max(0.0, rt_domain(did)%InflowLand(i))                 &
                       + max(0.0, rt_domain(did)%InflowGw(i))
            QT_total_i = QT_acc(i)                                              &
                       + T_air_i * max(0.0, rt_domain(did)%InflowLand(i))      &
                       + T_gw    * max(0.0, rt_domain(did)%InflowGw(i))

            if (Q_total_i .gt. min_flow) then
               T_in_i = QT_total_i / Q_total_i
            else
               T_in_i = T_e_i   ! headwater: initialise at equilibrium
            end if
            T_in_i = min(T_max, max(T_min, T_in_i))

            ! --- Exponential relaxation toward equilibrium along the reach ---
            depth_i = rt_domain(did)%HLINK(i)
            vel_i   = max(0.0, rt_domain(did)%velocity(i))

            if (depth_i .gt. min_depth .and. vel_i .gt. 0.0) then
               K_rate      = K_exchange / (rho_w * c_p_w * depth_i)
               travel_time = rt_domain(did)%CHANLEN(i) / vel_i
               T_out_i = T_e_i + (T_in_i - T_e_i) * exp(-K_rate * travel_time)
            else
               T_out_i = T_in_i   ! very low flow: pass inflow temperature through
            end if
            T_out_i = min(T_max, max(T_min, T_out_i))

            rt_domain(did)%T_QSUM(i) = T_out_i

         end do   ! i

      deallocate(Q_acc, QT_acc)

      ! Save this timestep's temperatures: next call gathers them as the
      ! (1-step lagged) upstream inflow temperatures.
      rt_domain(did)%T_QSUM_prev = rt_domain(did)%T_QSUM

   end subroutine run_stream_temp_eq

end module module_stream_temperature_eq
