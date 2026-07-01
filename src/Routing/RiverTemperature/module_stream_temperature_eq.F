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
! MPI / cross-rank temperature (Approach B):
!   drive_CHANNEL_RSL decomposes the domain spatially: each rank owns a rectangular tile.
!   The push pattern (accumulate Q*T into the downstream reach) cannot cross rank boundaries
!   within a single timestep.  The fix used here is a 1-timestep lag:
!
!     1. At the start of each call, T_QSUM_prev is saved (these are the values that were
!        communicated across ranks by MPP_CHANNEL_COM_REAL at the end of the previous step).
!     2. Q_acc / QT_acc are pre-seeded for ALL reaches (local + remote) using T_QSUM_prev.
!        Remote upstream reaches will never have their contribution updated this timestep,
!        so they retain T_QSUM_prev — a 1-step lag which is negligible at hourly dt.
!     3. The ORDER loop iterates only over local reaches (yw_mpp_nlinks / nlinks_index in
!        MPP mode); each local reach pushes the DELTA Q*(T_out - T_QSUM_prev) to replace
!        its pre-seeded contribution with the freshly computed T_out.
!     4. After the loop, MPP_CHANNEL_COM_REAL synchronizes T_QSUM across rank boundaries
!        so that the pre-seed for the next timestep carries correct cross-rank temperatures.
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
   use module_mpp_land,           only: MPP_CHANNEL_COM_REAL
#endif

   implicit none
   private
   public :: run_stream_temp_eq

   ! --- Physical constants ---
   real, parameter :: rho_w      = 1000.0   ! water density [kg/m3]
   real, parameter :: c_p_w      = 4182.0   ! water specific heat [J/kg/°C]

   ! --- Default model parameters (can be promoted to namelist later) ---
   real, parameter :: K_exchange = 150.0    ! bulk surface heat exchange coeff [W/m2/°C]
                                            ! Controls relaxation rate: K_rate = K_exchange/(rho_w*c_p*depth)
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

contains

   !---------------------------------------------------------------------------
   ! Subroutine: run_stream_temp_eq
   !
   ! Computes equilibrium stream temperature for every reach in domain `did`
   ! and stores results in rt_domain(did)%T_QSUM [°C].
   !---------------------------------------------------------------------------
   subroutine run_stream_temp_eq(did)
      integer, intent(in) :: did

      integer :: i, iorder, k, j_down
      integer :: nlinks, nlakes, maxorder, nsoil
#ifdef MPP_LAND
      integer :: iyw
#endif

      real    :: T_air_i, SW_net_i, T_e_i, T_in_i, T_out_i
      real    :: K_rate, travel_time, depth_i, vel_i, Q_this
      real    :: Q_total_i, QT_total_i, albedo_i, shade_i
      real    :: Tair_fallback, SW_fallback   ! domain-mean forcing for unmapped reaches

      real, allocatable :: Q_acc(:)         ! upstream flow accumulated at each reach [m3/s]
      real, allocatable :: QT_acc(:)        ! upstream Q*T accumulated at each reach  [m3/s * °C]

      nlinks   = rt_domain(did)%NLINKS
      nlakes   = rt_domain(did)%NLAKES
      nsoil    = noah_lsm%nsoil

      ! In MPP mode rt_domain%MAXORDER is a domain-wide fill value (-9999), not
      ! the local-partition maximum.  Compute maxorder from the local ORDER array.
      maxorder = 0
      do i = 1, nlinks
         if (rt_domain(did)%ORDER(i) .gt. maxorder) maxorder = rt_domain(did)%ORDER(i)
      end do

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
      ! Allocate accumulators.  In MPP mode, also save T_QSUM from the previous
      ! timestep (already communicated across ranks) and pre-seed Q_acc/QT_acc
      ! for ALL reaches so that cross-rank upstream contributions enter as a
      ! 1-timestep lag rather than being dropped entirely (see module header).
      ! -----------------------------------------------------------------------
      allocate(Q_acc(nlinks),  source=0.0)
      allocate(QT_acc(nlinks), source=0.0)

#ifdef MPP_LAND
      ! Pre-seed: walk every reach (including remote cross-rank reaches whose
      ! QLINK is valid after drive_CHANNEL_RSL's MPP_CHANNEL_COM_REAL) and
      ! push its previous-timestep temperature to its downstream accumulator.
      ! rt_domain(did)%T_QSUM_prev holds T_QSUM as communicated at the end of
      ! the previous timestep, so cross-rank upstream temperatures are included.
      do i = 1, nlinks
         j_down = int(rt_domain(did)%TO_NODE(i))
         if (j_down .ge. 1 .and. j_down .le. nlinks) then
            Q_this         = max(0.0, rt_domain(did)%QLINK(i,1))
            Q_acc(j_down)  = Q_acc(j_down)  + Q_this
            QT_acc(j_down) = QT_acc(j_down) + Q_this * rt_domain(did)%T_QSUM_prev(i)
         end if
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
                        if (rt_domain(did)%Tair_fine_1D(i) .gt. 200.0 .and. &
                            rt_domain(did)%Tair_fine_1D(i) .lt. 350.0) then
                           T_air_i  = rt_domain(did)%Tair_fine_1D(i) - 273.15
                           SW_net_i = max(0.0, rt_domain(did)%SolarRad_fine_1D(i))
                        else
                           T_air_i  = Tair_fallback - 273.15
                           SW_net_i = SW_fallback
                        end if
                        T_air_i  = min(T_max, max(T_min, T_air_i))
                        albedo_i = min(1.0, max(0.0, rt_domain(did)%S_Albedo(i)))
                        shade_i  = min(1.0, max(0.0, rt_domain(did)%F_Shade(i)))
                        SW_net_i = SW_net_i * (1.0 - albedo_i) * (1.0 - shade_i)
                        T_e_i    = T_air_i + SW_net_i / K_sw
                        rt_domain(did)%T_QSUM(i) = min(T_max, max(T_min, T_e_i))
                     end if
                  end associate
               end if
            end if
         end do
      end if

      ! -----------------------------------------------------------------------
      ! Step 3: Route temperature downstream reach by reach (topological order).
      !
      !   "Push" pattern: after computing T_out for reach i, we update the
      !   downstream reach's accumulator.
      !
      !   MPP mode (see module header for full explanation):
      !     Q_acc/QT_acc are pre-seeded from T_QSUM_prev for all reaches.
      !     The ORDER loop iterates only over local reaches (nlinks_index).
      !     Each local reach pushes the DELTA Q*(T_out - T_prev) rather than
      !     Q*T_out, so the pre-seeded cross-rank contribution is preserved
      !     for reaches whose upstream is on a different rank.
      !
      !   Serial mode: standard Q*T push (identical to pre-MPP behavior).
      !
      !   Lateral inflows (surface runoff + groundwater) enter at T_air as a
      !   first-order approximation.  A future improvement would use soil
      !   temperature (TSLB) for groundwater once NoahMP state is coupled in.
      ! -----------------------------------------------------------------------
      do iorder = 1, maxorder

#ifdef MPP_LAND
         do iyw = 1, rt_domain(did)%yw_mpp_nlinks
            i = rt_domain(did)%nlinks_index(iyw)
#else
         do i = 1, nlinks
#endif

            if (rt_domain(did)%ORDER(i) .ne. iorder) cycle

            ! --- Lake reaches: temperature already set in Step 2; just push ---
            if (rt_domain(did)%TYPEL(i) .eq. 1) then
               j_down = int(rt_domain(did)%TO_NODE(i))
               if (j_down .ge. 1 .and. j_down .le. nlinks) then
                  Q_this = max(0.0, rt_domain(did)%QLINK(i,1))
#ifdef MPP_LAND
                  ! Replace pre-seeded T_QSUM_prev contribution with Step 2 lake temp
                  QT_acc(j_down) = QT_acc(j_down) &
                                 + Q_this * (rt_domain(did)%T_QSUM(i) - rt_domain(did)%T_QSUM_prev(i))
#else
                  Q_acc(j_down)  = Q_acc(j_down)  + Q_this
                  QT_acc(j_down) = QT_acc(j_down) + Q_this * rt_domain(did)%T_QSUM(i)
#endif
               end if
               cycle
            end if

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

            ! --- Equilibrium temperature ---
            T_e_i = T_air_i + SW_net_i / K_sw
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

            ! --- Push outflow temperature to downstream reach accumulator ---
            j_down = int(rt_domain(did)%TO_NODE(i))
            if (j_down .ge. 1 .and. j_down .le. nlinks) then
               Q_this = max(0.0, rt_domain(did)%QLINK(i,1))
#ifdef MPP_LAND
               ! Replace pre-seeded T_QSUM_prev contribution with current T_out_i.
               ! Q_acc(j_down) is already correct from the pre-seed; only QT needs updating.
               QT_acc(j_down) = QT_acc(j_down) + Q_this * (T_out_i - rt_domain(did)%T_QSUM_prev(i))
#else
               Q_acc(j_down)  = Q_acc(j_down)  + Q_this
               QT_acc(j_down) = QT_acc(j_down) + Q_this * T_out_i
#endif
            end if

#ifdef MPP_LAND
         end do   ! iyw
#else
         end do   ! i
#endif
      end do      ! iorder

      deallocate(Q_acc, QT_acc)
#ifdef MPP_LAND
      ! Synchronize T_QSUM across rank boundaries so the next timestep's pre-seed
      ! carries correct cross-rank temperatures.
      call MPP_CHANNEL_COM_REAL(rt_domain(did)%Link_location,  &
                                 rt_domain(did)%IXRT,           &
                                 rt_domain(did)%JXRT,           &
                                 rt_domain(did)%T_QSUM,         &
                                 nlinks, 99)
      ! Store communicated T_QSUM for use as pre-seed on the next timestep.
      rt_domain(did)%T_QSUM_prev = rt_domain(did)%T_QSUM
#endif

   end subroutine run_stream_temp_eq

end module module_stream_temperature_eq
