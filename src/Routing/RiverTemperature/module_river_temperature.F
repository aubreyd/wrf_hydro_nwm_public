!  Program Name: module_river_temperature
!  Author(s)/Contact(s): Reza Abdi (rabdi@ucar.edu)
!  Abstract:
!  History Log: First attept, Mar. 21
!
!  Usage:
!  Parameters:
!  Input Files:
!
!  Output Files:
!
!  Condition codes:
!
!  User controllable options: <if applicable>

MODULE module_river_temperature


	!use module_channel_routing, only: temperature_sim_H_inputs
	!use disaggregateWeatherModule,  only: fineWeatherData
	!use module_RT_data, only: rt_domain
	!use config_base, only: nlst
	IMPLICIT NONE


type  :: river_temperature_type
		!type (temperature_sim_H_inputs), pointer :: hydraulics_sim => null()
		!type (fineWeatherData), pointer :: weather_sim => null()

     logical :: pointer_allocation_guard = .false.

     contains
     procedure :: init => river_temperature_init
     procedure :: destroy => river_temperature_destroy
     !procedure :: run => river_temperature_run
     procedure :: runRiverTempModel => river_TEMPERATURE_SOLVER

end type river_temperature_type

type(river_temperature_type)  :: riverTemp

REAL, parameter ::  rho_water = 1000.
REAL, parameter ::  c_water = 4182.

contains

!subroutine river_temperature_init(this, NLINKS, IXRT, JXRT, reach_number, ChanCondConst, QSTRMVOLRT,  &
!					  QINFLOWBASE, QSUM,  XS_Peremeter, Bw, Tw, velocity, HLINK, ChSSlp, &
!					  CHANXI, CHANYJ, ZELEV, CHANLEN, CH_NETLNK, QLINK, DT_STEPS, DTCT, node_area, TO_NODE, TYPEL, &
!           IX, JX,  T2, q2x, u, short, CH_NETRT, LAKE_MSKRT, nsoil, ims,ime, jms,jme, TMN, TSK, TSLB, &
!            Tair_fine_1D, RelH_fine_1D, windSpd_fine_1D, SolarRad_fine_1D, TMN_1D_fine, TSK_1D_fine,TSLB_1D_fine, T_QSUM)


subroutine river_temperature_init(this, T_QSUM, T_QSUM_prev, NetRad, NetRad_prev, water_temp_counter)

	implicit none
	class(river_temperature_type), INTENT(inout)       :: this              ! the type object being initialize
        REAL, INTENT(inout), DIMENSION(:)             :: T_QSUM
				REAL, INTENT(inout), DIMENSION(:)             :: T_QSUM_prev
				REAL, INTENT(inout), DIMENSION(:)             :: NetRad
				REAL, INTENT(inout), DIMENSION(:)             :: NetRad_prev
				INTEGER, INTENT(inout)                        :: water_temp_counter

end subroutine river_temperature_init
! --------------------------------------------------------------------------
! river_temperature Destructor
subroutine river_temperature_destroy(this)
        implicit none
        class(river_temperature_type), intent(inout) :: this ! object being destroyed
 end subroutine river_temperature_destroy

! --------------------------------------------------------------------------
!subroutine river_temperature_run(this, did)

!	implicit none
!  class(river_temperature), INTENT(inout)             :: this

  ! Embeddedness and Substrate particle size (mm) respectively
!	REAL,    DIMENSION(:), allocatable                  :: C_Em, P_size
!	REAL,    DIMENSION(:), allocatable                  :: shadeF, sur_albedo, viewToSky, cloud   ! Shading factor and surface albedo for the reach
	!REAL,    INTENT(INOUT), DIMENSION(:)                :: T_QSUM

! INTEGER                                             ::  did  ! for the rt_domain


!  call river_TEMPERATURE_SOLVER(nlst(did)%dxrt0,      	        &
!				rt_domain(did)%DT_STEPS,        &
!                                rt_domain(did)%IXRT,            &
!                                rt_domain(did)%JXRT,            &
!                                rt_domain(did)%Tair_fine_1D,    &
!				rt_domain(did)%RelH_fine_1D,    &
!				rt_domain(did)%windSpd_fine_1D, &
!				rt_domain(did)%C_Cloud,         &
!				rt_domain(did)%TSLB_1D_fine,    &
!				rt_domain(did)%SolarRad_fine_1D,&
!				rt_domain(did)%XS_Peremeter,    &
!				rt_domain(did)%HLINK,	        &
!				rt_domain(did)%ChSSlp,	        &
!				rt_domain(did)%ZELEV,	        &
!                                rt_domain(did)%velocity,        &
!				rt_domain(did)%Tw,	        &
!				rt_domain(did)%QSUM,	        &
!				rt_domain(did)%QINFLOWBASE,     &
!				rt_domain(did)%QSTRMVOLRT,      &
!				rt_domain(did)%TMN_1D_fine,     &
!				rt_domain(did)%TSLB_1D_fine,    &
!				rt_domain(did)%TSK_1D_fine,     &
!                                rt_domain(did)%S_Albedo,        &
!                		rt_domain(did)%F_Shade,         &
!				rt_domain(did)%ChanCondConst,   &
!                                rt_domain(did)%C_Emb,	        &
!	                	rt_domain(did)%P_Size,          &
!				rt_domain(did)%CHANXI,          &
!				rt_domain(did)%CHANYJ,	        &
!				rt_domain(did)%CH_NETLNK,       &
!                                rt_domain(did)%CH_NETRT,        &
!                                rt_domain(did)%QLINK,	        &
!				rt_domain(did)%DT_STEPS,        &
!				rt_domain(did)%DTCT,	        &
!				rt_domain(did)%node_area,       &
!				rt_domain(did)%TO_NODE,	        &
!				rt_domain(did)%TYPEL,	        &
!				rt_domain(did)%T_QSUM,	        &
!                               rt_domain(did)%NLINKS           &
!                 		)

!end subroutine river_temperature_run

! ------------------------------------------------
!   SUBROUTINE RIVER_TEMPERATURE_SOLVER
! ------------------------------------------------
Subroutine river_TEMPERATURE_SOLVER(this, dx, DT, DTCT, IXRT, JXRT,  &
	      T_a, H, sfcp, V_wind, C_L, T_bed, s_in, Pw_river, y_w, S_0, z_elv, &
				V_w, W_w, Q_str, Q_GW_cms, Q_Trib_cms, T_GW, T_hyp, T_Trib, T_upstream, albedo, &
				SF, k_s_0, C_E, P_s, CHANXI, CHANYJ, CH_NETLNK, CH_NETRT, QLINK, &
				node_area, TO_NODE, TYPEL, T_QSUM_local, T_QSUM_prev, NetRad, NetRad_prev, water_temp_counter, totLength)

	implicit none
        class(river_temperature_type), INTENT(inout)             :: this
	! -------- DECLARATIONS ------------------------
        REAL, INTENT(IN)  :: dx		! Space intervals which should be based on the HYDRO grid size  -===================> CHANLEN
	REAL, INTENT(IN)  :: DT		! time intervals which in most of cases should be 1 hour           - ??? should be based on the hydraulic routing module
	REAL, INTENT(IN)  :: DTCT
	!INTEGER, INTENT(IN) :: t 			! the current and previous time steps for the simulation, t=2
	INTEGER, INTENT(IN)  :: totLength          	! Total length (# of cells) of the domain for simulation -===========================> NLINKS
	INTEGER, INTENT(IN)  :: IXRT, JXRT
	REAL, INTENT(IN), DIMENSION(:)    :: T_a		! Air temperature 						               - As an input
	REAL, INTENT(IN), DIMENSION(:)    :: H		        ! Humidity (%) - Specific humidity			                       - As an input
	REAL, INTENT(IN), DIMENSION(:)    :: sfcp         ! surface pressure (Pa)
	REAL, INTENT(IN), DIMENSION(:)    :: V_wind	        ! Wind speed (m/s)							       - As an input
	REAL, INTENT(IN), DIMENSION(:)	  :: C_L	        ! Cloudiness (0 to 1)					                       - As an input - should find it from somewhere!
	REAL, INTENT(IN), DIMENSION(:)	  :: T_bed	        ! Conduction Layer Temperature (C)			                       - As an input
	REAL, INTENT(IN), DIMENSION(:)    :: s_in		! Incoming sortwave radiation, as an input
	REAL, INTENT(IN), DIMENSION(:)    :: Pw_river	        ! Wetted perimeter in the river channle for the desird grid (m)                 =====> XS_Peremeter
	REAL, INTENT(IN), DIMENSION(:)    :: y_w		! Water Column Depth (m)                                    	   =====> HLINK
	REAL, INTENT(IN), DIMENSION(:)	  :: S_0		! Stream gradient (unitless)							=====> ChSSlp !TLE: changed to So. So always seems to be 0.01? seems high but there must be a reason
        REAL, INTENT(IN), DIMENSION(:)	  :: z_elv		! Elevation of station where met. data was obtained				=====> ZELEV
	REAL, INTENT(IN), DIMENSION(:)	  :: V_w		! Average Flow Velocity (m/s)							=====> velocity
	REAL, INTENT(IN), DIMENSION(:)	  :: W_w		! Wetted width (Top; m)								=====> Tw
	! Discharge data, instream, groundwatr flow, triburaty
	REAL, INTENT(IN), DIMENSION(:)	  :: Q_str		! Instream flow in the river channel - from sll sources (csm)                   =====> QSUM
!	REAL, INTENT(IN), DIMENSION(:,:)  :: Q_GW		! Groundwater flow in the gridded based (cms)					=====> QINFLOWBASE
!	REAL, INTENT(IN), DIMENSION(:,:)  :: Q_Trib		! Tributarty surface inflow in the gridded based (cms)			        =====> QSTRMVOLRT
	REAL, INTENT(IN), DIMENSION(:)  :: Q_GW_cms		! Groundwater flow in the gridded based (cms)					=====> InflowGw
	REAL, INTENT(IN), DIMENSION(:)  :: Q_Trib_cms		! Tributarty surface inflow in the gridded based (cms)			        =====> InflowLand

	REAL, INTENT(IN), DIMENSION(:)	  :: T_GW		! Groundwater flow temperature in the gridded based (C)		         	- As an input
	REAL, INTENT(IN), DIMENSION(:)	  :: T_hyp		! Hyporheic Flow temperature in the gridded based (C)		         	- As an input
	REAL, INTENT(IN), DIMENSION(:)	  :: T_Trib		! Tributary surface inflow temperature (C)					- As an input
	REAL, INTENT(IN), DIMENSION(:)	  :: T_upstream		! Temperature of flow from upstream node (C)					- As an input
	! Shading related parameters
	REAL, INTENT(IN), DIMENSION(:)    :: albedo		! albedo at the stream surface (0 to 1)						- ???
	REAL, INTENT(IN), DIMENSION(:)    :: SF			! Shading factor for each grid of the reach					- ???
	! Fluxes
	REAL                              :: SW_rad		! Shortwave radiation flux (W/m2)
	REAL                              :: LW_rad		! Longwave radiation flux (W/m2)
	REAL                              :: LH		        ! Latent heat energy flux (W/m2)
	REAL                              :: SH		        ! Sensible heat energy flux (W/m2)
	REAL                              :: CH		        ! Conduction flux (W/m2)
	! For hyporheic exchange incflow calculation & Porosity of the conduction layer (unitless)
	REAL, DIMENSION(IXRT,JXRT)        :: Q_hyp		! Darcy’s Law for Calculating for Hyporheic Flow (Domenico and Schwartz 1990)
	REAL, DIMENSION(IXRT,JXRT)        :: Disp_coef
        REAL		          	  :: A_s		! Cross-Sectional Across Seepage Face (m2)
	REAL			          :: k_s		! Composite Substrate Hydraulic Conductivity (m/s)
	REAL, INTENT(IN), DIMENSION(:)    :: C_E		! Embeddedness                                                          - As an input
	REAL, INTENT(IN), DIMENSION(:)    :: P_s		! Substrate particle size (mm)
	REAL, INTENT(IN), DIMENSION(:)    :: k_s_0              ! Dominant Substrate Hydraulic Conductivity (m/s)                       ===========> ChanCondConst

	REAL                              :: h_D		! Hydraulic Head for Darcy Calculation at location i (m)
	REAL                              :: h_D_dn		! Hydraulic Head for Darcy Calculation at location i+1 (m)
	REAL                              :: Dummy1		! For calcluating the hydraulic conductiveity
	REAL                              :: Dummy2		! For calcluating the hydraulic conductiveity
	REAL                              :: d_h		! Change in hydraulic head per unit distance
	! Two variables for the calculations
	!REAL, INTENT(OUT)                :: rho_water          ! Density of water (kg/m3)
	!REAL, INTENT(OUT)                :: c_water	        ! Specific heat capacity of water (J/kgC)
	! Component of the main equaion of the simulated water temperature
	REAL                              :: Advection  	!
	REAL                              :: Dispersion 	!
	REAL                              :: self
	REAL                              :: FluxTemp	        !
	REAL                              :: FluxDelta        ! TLE: added this as a tracker variable to see T update in this timestep
	!REAL                              :: Disp_coef      	! Dispersion Coefficient (m2/s)
	REAL                              :: Shear_Velocity	! Shear Velocity (m/s)
	REAL                              :: SlopeWS		! Water surface slope
	!REAL, INTENT(IN), DIMENSION(totLength)	  :: Z			                    ! Side slope (fraction)
	REAL,                   DIMENSION(totLength)          :: TLateral
	INTEGER, INTENT(IN),    DIMENSION(:)                  :: CHANXI, CHANYJ
	INTEGER, INTENT(IN),    DIMENSION(:,:)                :: CH_NETLNK
        INTEGER, INTENT(IN),    DIMENSION(:,:)                :: CH_NETRT
	REAL,                   DIMENSION(totLength)          :: QLateral
	REAL,    INTENT(IN),    DIMENSION(totLength,2)        :: QLINK
	REAL,    INTENT(in),    DIMENSION(:)                   :: node_area
	INTEGER, INTENT(in),    DIMENSION(:)                  :: TO_NODE
	INTEGER, INTENT(IN),    DIMENSION(:)                  :: TYPEL
	REAL,    INTENT(INOUT), DIMENSION(:)                  :: T_QSUM_local	            ! -- FINAL PRODUCT
	REAL,    INTENT(INOUT), DIMENSION(:)                  :: T_QSUM_prev  ! from previous timestep
	REAL,    INTENT(INOUT), DIMENSION(totLength)          :: NetRad, NetRad_prev	! Net heat flux from 5 sources (W/m2)
	INTEGER, INTENT(INOUT)                                :: water_temp_counter   !!!TLE: counter to keep track of how many routing time steps have elapsed
	INTEGER,                DIMENSION(totLength)          :: BC_Grid		    ! -- This is the new way I am thinking of to detec the upstream grids
	INTEGER,                DIMENSION(totLength)          :: FROM_NODE, N_UPSTR, node_order, node_headwater, node_branch
	INTEGER                                               :: neighbours, i, j, k, loop_break, c_index
	INTEGER                                               :: node_current, c_hw, c_br, denom
        ! Simulated water temeprature (C)
	REAL,                   DIMENSION(totLength)          :: T_w_sim	! Simulated water temperature for the reach in space & time (C)
  REAL                                                  :: Q_strTot, T_strTot, Q_strTmp, T_up, T_down, branch1_temp
	REAL,                   DIMENSION(totLength)          :: viewSky 		! View to sky factor (1-SF)
	REAL,                   DIMENSION(totLength)          :: T_w_sim_temp, s1t, s2t, branch_sum

        !REAL                                                  :: DT_STEPS               !-- number of timestep in routing
	!REAL,                   DIMENSION(totLength)                :: Q_GW_cms
	!REAL,                   DIMENSION(totLength)                :: Q_Trib_cms
	REAL,                   DIMENSION(totLength)                :: percent_Q_LINK, percent_Q_strTot, percent_Q_Lateral !TLE: temporary to look at flow fractions

	logical :: first_time_step_DT = .true.
	logical :: first_time_step_DTCT = .true.

	viewSky = 1. - SF
        !DT_STEPS = INT(DT/DTCT)   !-- fix the timestep
        ! ++++++++++++++++++++++++++++update this part
	! Assigning the upstream water temeprature boundary condition to the main matrix
	!-- Setting everything to zero first.
	Disp_coef = 0.
  Q_hyp = 0.
	s1t = 0.
	s2t = 0.

	!print *, "DTCT = ", DTCT

	! Calculating hyporheic fluxes
	Do i = 1, totLength
	    ! Find hyporheic fluxes
		  ! Ratio Conductivity of dominant sunstrate
		  Dummy1 = k_s_0(CH_NETLNK(CHANXI(i),CHANYJ(i))) * (1.0 - C_E(i))
	    !Dummy1 = k_s_0 * (1.0 - C_E(i))

	    ! Ratio Conductivity of sand - low range
	    ! 0.00002 is the embedded Substrate Hydraulic Conductivity for Silt/Sand (m/s)
	    Dummy2 = 0.00002 * C_E(i)
	    ! True cond. (m/s); Composite Substrate Hydraulic Conductivity
	    k_s = Dummy1 + Dummy2

	    h_D = y_w(CH_NETLNK(CHANXI(i),CHANYJ(i)))
	    h_D_dn = y_w(TO_NODE(i)) - (dx * S_0(CH_NETLNK(CHANXI(i),CHANYJ(i))))  !!!!!TLE: I'm pretty sure we want channel slope here, but Reza is using ChSSlp

	    !Change in hydraulic head per unit distance
	    d_h = (h_D - h_D_dn) / dx
			!	  end if


			! Cross-Sectional Area of the Seepage Face (m2)
			A_s = dx * S_0(CH_NETLNK(CHANXI(i),CHANYJ(i))) * Pw_river(CH_NETLNK(CHANXI(i),CHANYJ(i))) !!!! TLE: Again, I think we want channel slope
			! Darcy’s Law for Calculating for Hyporheic Flow (cms
			!		Q_hyp(CHANXI(i),CHANYJ(i)) = A_s * k_s * d_h
			Q_hyp(CHANXI(i),CHANYJ(i)) = 0

			!		print *, "S_0 = ", S_0(i)
			!		print *, "d_h = ", d_h
			!		print *, "A_s = ", A_s
			!		print *, "Q_hyp(i) = ", Q_hyp(CHANXI(i),CHANYJ(i))

      ! Calculating dispersion coefficient
			! Calculate dispersion using equations adapted from SubMcCormick1 - HeatSourceModel source code
			! (Martin and McCutcheon 1999)
			! To enable simulating the cooling/warming effect of riffle-pool
			SlopeWS = d_h
			if (SlopeWS < 0.0) then
			    ! For adverse slopes
					Shear_Velocity = V_w(CH_NETLNK(CHANXI(i),CHANYJ(i)))
			else
			    Shear_Velocity = sqrt(9.8 * y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) * SlopeWS)
			end if

			! Calculating wetted width for the river channel (top with - m)
			!W_w = W_b(node) + 2.0 * Z(node) * y_w(node) => Will grab it from the routing module

			! Physical Dispersion Coefficient (Fischer et. al. 1979)
			! Reza: Changed 0.011 to 0.06 in below eqn to account for natural (not straight) streams
    	Disp_coef(CHANXI(i),CHANYJ(i)) = (0.06 * (V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) ** 2.0) * (W_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) ** 2.0))/ &
							(y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) * Shear_Velocity)

			! A condition based on the adaptd methods from Martin & McCutcheon (1999) also considered with Boyd & Kasper (2003)
			! 3600 is the assumed dt in sectonds
			if ((Disp_coef(CHANXI(i),CHANYJ(i)) * ((3600.0 / dx) ** 2.0)) > 0.5) then
		      Disp_coef(CHANXI(i),CHANYJ(i)) = (0.45 * (dx ** 2.0)) / 3600.0
		  end if

		  if (y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) .EQ. 0.) then
		      Disp_coef(CHANXI(i),CHANYJ(i)) = 0.
		  end if
	end Do

			!!!TLE: trying an alternative to Reza's monster 'if' block. We need to establish FROM_NODE
			!      in order to use the MacCormick solver Reza presumably originally intended to use. I can't
			!      figure out where he got the alternative to the MacCormick solver, but whatever it is it
			!      doesn't seem to be working that well.

		  FROM_NODE = 0.
	  	N_UPSTR = 0.
  		node_order = 0.

			!loop over TO_NODE, use it to construct FROM_NODE, N_UPSTR (number of upstream nodes)
			! this should be enough info to construct a up->down order...?
	    DO i = 1,totLength
			    FROM_NODE(TO_NODE(i)) = FROM_NODE(TO_NODE(i)) + i
					N_UPSTR(TO_NODE(i)) = N_UPSTR(TO_NODE(i)) + 1
			end DO

#ifdef HYDRO_D
			print *, "FROM_NODE", FROM_NODE
			print *, "TO_NODE", TO_NODE
			print *, "N_UPSTR", N_UPSTR
#endif

			!second loop to construct list of headwater nodes, list of branch nodes
			c_hw = 1.
			c_br = 1.
			node_headwater = -1.
			node_branch = -1.
			DO i = 1, totLength

					! check if headwater
			    if (FROM_NODE(i) .eq. 0) then
					    node_headwater(c_hw) = i
					    c_hw = c_hw + 1
					end if

					!check if branch
					if (N_UPSTR(i) .gt. 1) then
					    node_branch(c_br) = i
							c_br = c_br + 1
					end if
			end DO

#ifdef HYDRO_D
			print *, "node_branch = ", node_branch
			print *, "node_headwater = ", node_headwater
#endif

			!construct the up->downstream node order
			loop_break = 0.
			c_index = 1.
			c_hw = 1.
			c_br = 1.
			node_current = node_headwater(c_hw)
			DO WHILE (loop_break .eq. 0)
					!print *, "c_index =", c_index
					!print *, "node_current =", node_current
					!print *, "N_UPSTR(TO_NODE(node_current)) = ", N_UPSTR(TO_NODE(node_current))
			    if (c_index .eq. totLength) then
						  node_order(c_index) = node_current
					    loop_break = 1
					end if
					if (TO_NODE(node_current) .lt. 0) then
						  c_hw = c_hw + 1
							node_order(c_index) = node_current
						  node_current = node_headwater(c_hw)
						  c_index = c_index + 1
							cycle
					end if
					if (N_UPSTR(TO_NODE(node_current)) .eq. 1) then
!						  print *, "node_current = ", node_current
!							print *, "TO_NODE(node_current) = ", TO_NODE(node_current)
					    node_order(c_index) = node_current
						  c_index = c_index + 1
							node_current = TO_NODE(node_current)
					else if (N_UPSTR(TO_NODE(node_current)) .gt. 1) then
						  if (node_branch(c_br) .gt. 0) then                          !move to the next headwater
!								  print *, "node_current = ", node_current
!									print *, "TO_NODE(node_current) = ", TO_NODE(node_current)
							    node_order(c_index) = node_current
									c_index = c_index + 1
							    c_hw = c_hw + 1
									c_br = c_br + 1
							    node_current = node_headwater(c_hw)
							else                                                        !continue to mainstem
							    node_order(c_index) = node_current
								  c_index = c_index + 1
								  node_current = TO_NODE(node_current)
							end if
					end if
			end DO
#ifdef HYDRO_D
			print *, "node_order = ", node_order
#endif

	if (first_time_step_DT) then	!-- To do this only one time!

			!initialize homogeneous stream channel temperature
			T_w_sim_temp = T_Trib - 273.2
			T_QSUM_prev = T_Trib - 273.2
			!T_w_sim_temp = T_a - 273.2
			!T_QSUM_prev = T_a - 273.2

			!starting off with no rad in previous timestep (not sure if that's a big deal or not)
			NetRad_prev = 0.

			first_time_step_DT = .false.

  else
	    !-- It's supposed to be in Deg C
      T_w_sim_temp = T_QSUM_local

  end if

  ! Calculate net radiation in this time step for all nodes
  Do k = 1 , totLength
			i = node_order(k)
	    ! Calculate the heat fluxes

			SW_rad = SW_RAD_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), s_in, SF, albedo)
			LW_rad = LW_RAD_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), T_w_sim_temp, T_a, H, sfcp, C_L, viewSky)
			LH     = LatentHeat_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), SW_rad, LW_rad, T_w_sim_temp, T_a, H, sfcp, V_wind, z_elv)
			SH     = SensibleHeat_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), T_w_sim_temp, LH, T_a, H, sfcp, z_elv)
			!!!TLE: turning off channel bed sediment conduction for now.
			!CH     = BedSediment_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)),  dx, T_w_sim_temp, P_s, C_E, Pw_river, T_bed)
			CH     = 0.0

			! Net radiation -> Heat Flux
			NetRad(CH_NETLNK(CHANXI(i),CHANYJ(i))) = (SW_rad + LW_rad + LH + SH + CH) / &
			    (rho_water * c_water * y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))

			! Calculate mixture temperature for each node and time step
#ifdef HYDRO_D
			print *, "k = ", k
			print *, "i = ", i
			print *, "-----------ADVECTION-----------"
			print *, "CH_NETLNK(CHANXI(i),CHANYJ(i)) = ", CH_NETLNK(CHANXI(i),CHANYJ(i))
			print *, "-----SW_rad-----"
			print *, "SW_rad = ", SW_rad
			print *, "-----LW_rad-----"
			print *, "LW_rad = ", LW_rad
			print *, "-----LH-----"
			print *, "LH = ", LH
			print *, "-----SH-----"
			print *, "SH = ", SH
			print *, "-----CH-----"
			print *, "CH = ", CH
			print *, "------------"
			print *, "NetRad(CH_NETLNK(CHANXI(i),CHANYJ(i))) = ", NetRad(CH_NETLNK(CHANXI(i),CHANYJ(i)))
#endif

  end do

#ifdef HYDRO_D
	print *, "QLINK = ", QLINK(:,1)
	print *, "T_w_sim_temp = ", T_w_sim_temp
	print *, "shortwave = ", s_in
	print *, "NetRad = ", NetRad
#endif
 	!!!TLE: Incorporating the MacCormick (1969) finite difference approximation of
	!      the one-dimensional heat transfer equation (as documented in Boyd and Kasper (2003), p.68, Eq. 2-119 -- 2-122)
	!MACCORMICK STEP 1
	branch_sum = 0.
	branch1_temp = 0.
	denom = 0.
  Do k = 1,totLength
			!converter to keep node order correct
      i = node_order(k)

			if (TO_NODE(i) .ne. node_order(k+1)) then
				 branch_sum(TO_NODE(i)) = branch_sum(TO_NODE(i)) + T_QSUM_prev(i)
				! branch1_temp = branch1_temp + T_QSUM_prev(i)
				! denom = denom + 1
			end if

			!set up and downstream T's
			if (N_UPSTR(i) .eq. 0) then
				  T_up = (T_trib(i)-273.2)
			else if (N_UPSTR(i) .gt. 1) then
				  T_up = (branch_sum(i) + T_QSUM_prev(node_order(k-1))) / N_UPSTR(i)
	    else
			    T_up = T_QSUM_prev(node_order(k-1))
			end if

			if (TO_NODE(i) .lt. 0) then
			    T_down = (T_trib(i)-273.2)
			else
				  if (TO_NODE(i) .ne. node_order(k+1)) then
						  T_down = T_QSUM_prev(TO_NODE(i))
				else
			    	  T_down = T_QSUM_prev(node_order(k+1))
				end if
			end if

			!calculate s1t
			s1t(i) = (-1. * V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) * ((T_down - T_QSUM_prev(i))/dx)) + &
			      (Disp_coef(CHANXI(i),CHANYJ(i)) * ((T_down-(2.*T_QSUM_prev(i))+T_up)/(dx**2.))) + &
						(NetRad_prev(CH_NETLNK(CHANXI(i),CHANYJ(i))))

#ifdef HYDRO_D
				print *, "---MacCormick Step 1---"
				print *, "term1 = ", (-1. * V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) * ((T_down - T_QSUM_prev(i))/dx))
				print *, "V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) =", V_w(CH_NETLNK(CHANXI(i),CHANYJ(i)))
				print *, "T_QSUM_prev(i) =", T_QSUM_prev(i)
				print *, "T_down =", T_down
				print *, "term2 = ", (Disp_coef(CHANXI(i),CHANYJ(i)) * ((T_down-(2.*T_QSUM_prev(i))+T_up)/(dx**2.)))
				print *, "term3 = ", (NetRad_prev(CH_NETLNK(CHANXI(i),CHANYJ(i))))
				print *, "s1t = ", s1t(i)
#endif

			!first approximation of T_w
			T_w_sim_temp(i) = T_QSUM_prev(i) + (s1t(i) * DT)

	end Do

	!MACCORMICK STEP 2
	branch_sum = 0.
	branch1_temp = 0.
	denom = 0.
	Do k = 1,totLength
			!converter to keep node order correct
			i = node_order(k)

			if (TO_NODE(i) .ne. node_order(k+1)) then
				 branch_sum(TO_NODE(i)) = branch_sum(TO_NODE(i)) + T_w_sim_temp(i)
			end if

			!set up and downstream T's
			if (N_UPSTR(i) .eq. 0) then
				  T_up = (T_trib(i)-273.2)
			else if (N_UPSTR(i) .gt. 1) then
				  T_up = (branch_sum(i) + T_w_sim_temp(node_order(k-1))) / N_UPSTR(i)
			else
					T_up = T_w_sim_temp(node_order(k-1))
			end if
			if (TO_NODE(i) .lt. 0) then
				  T_down = (T_trib(i)-273.2)
			else
					if (TO_NODE(i) .ne. node_order(k+1)) then
							T_down = T_w_sim_temp(TO_NODE(i))
					else
							T_down = T_w_sim_temp(node_order(k+1))
					end if
			end if

			!calculate s2t
			s2t(i) = (-1. * V_w(i) * ((T_w_sim_temp(i) - T_up)/dx)) + &
			      (Disp_coef(CHANXI(i),CHANYJ(i)) * ((T_down - (2.*T_w_sim_temp(i)) + T_up)/(dx**2.))) + &
						(NetRad(i))

#ifdef HYDRO_D
				print *, "---MacCormick Step 2---"
				print *, "term1"
				print *, "v_w(i) = ", V_w(i)
				print *, "(T_w_sim_temp(i) - T_up) = ", (T_w_sim_temp(i) - T_up)
				print *, "term2 = ", (Disp_coef(CHANXI(i),CHANYJ(i)) * ((T_down - (2.*T_w_sim_temp(i)) + T_up)/(dx**2.)))
				print *, "term3 = ", (NetRad(i))
				print *, "s2t = ", s2t(i)
#endif

			!final T_w estimate
			T_w_sim_temp(i) = T_QSUM_prev(i) + (((s1t(i) + s2t(i))/2.)*DT)

			!convert back to Kelvin
			T_w_sim(CH_NETLNK(CHANXI(i),CHANYJ(i))) = T_w_sim_temp(i) + 273.2
 end Do


  Do i = 1,totLength


	   !-- getting a weighted average
	   if (CH_NETLNK(CHANXI(i),CHANYJ(i)) .gt. 0) then
		     if ((Q_Trib_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))+Q_GW_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))+Q_hyp(CHANXI(i),CHANYJ(i))) .EQ. 0.) then
		         TLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = 0.
		     else
		           TLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = ((Q_Trib_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))*T_Trib(i)) + &
                                                            (Q_GW_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))*T_GW(i)) + &
							    (Q_hyp(CHANXI(i),CHANYJ(i))*T_hyp(i))) / &
							    (Q_Trib_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))+Q_GW_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))+ &
							    Q_hyp(CHANXI(i),CHANYJ(i)))

		     end if

	       QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = Q_Trib_cms(CH_NETLNK(CHANXI(i),CHANYJ(i))) + Q_GW_cms(CH_NETLNK(CHANXI(i),CHANYJ(i))) + Q_hyp(CHANXI(i),CHANYJ(i))
				 !!!TLE Hack: we need non-zero flow in the outlet node for the first time step or it goes to NaN. set Qlateral to an arbitrary non-zero value
				 if (first_time_step_DTCT) then
				     if (TO_NODE(i) .le. 0) then
						     if (QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) .le. 0.) then
								     QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = 0.01
										 Tlateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = T_Trib(CH_NETLNK(CHANXI(i),CHANYJ(i)))
	!									 print *, "TLE hack"
								 end if
						 end if
				 end if
	   end if


	  ! TLE: Final round of re-writing (hopefully). We now have an explicit 'inflow from upstream'
		!      variable courtesy of Aubrey, so no more finagling with negatives. Can also get
		!      rid of the loops since that is now taken care of in module_channel_routing
		!      TODO: check in about the backwatering thing

		Q_strTot = Q_str(i)
		T_strTot = T_upstream(i)


    !-- to Deg K
	  T_strTot = T_strTot + 273.2

	  T_QSUM_local(i) = (((Q_strTot) * T_strTot) + & !(((Q_strTot*-1.) * T_strTot) + &
		 (QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1)*T_w_sim(i)) + &
		 (QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i)))*TLateral(i))) / &
		 ((Q_strTot)+QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1)+QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))))

#ifdef HYDRO_D
		!TLE TEMPORARY: calculate the relative contribution of each Q source (percentage of total Q)
		percent_Q_LINK(i) = (ABS(QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1)) / &
		  (ABS(Q_strTot)+ABS(QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1))+ABS(QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))))))*100
		percent_Q_strTot(i) = (ABS(Q_strTot) / &
 		  (ABS(Q_strTot)+ABS(QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1))+ABS(QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))))))*100
		percent_Q_Lateral(i) = (ABS(QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i)))) / &
 		  (ABS(Q_strTot)+ABS(QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1))+ABS(QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))))))*100

		   print *,"---------------------------"
				 print *,"i=", i
				 print *,"CHANXI(i)=",CHANXI(i)
				 print *,"CHANYJ(i)=",CHANYJ(i)
				 print *,"CH_NETLNK(CHANXI(i),CHANYJ(i))=", CH_NETLNK(CHANXI(i),CHANYJ(i))
				 print *,"Q_strTot=", Q_strTot
			   print *,"T_strTot=", T_strTot
			   print *,"QLINK(i,1)=", QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1)
!				 print *,"T_w_sim_temp=", T_w_sim_temp(i)
			   print *,"T_w_sim(i,2)=", T_w_sim(i)
			   print *,"T_QSUM_local(i)", T_QSUM_local(i)
				 print *,"Q_hyp(i)", Q_hyp(CHANXI(i),CHANYJ(i))
	!			 print *,"Q_GW(i)", Q_GW(CHANXI(i),CHANYJ(i))
				 print *,"Q_GW_cms(i)", Q_GW_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))
!				 print *,"Q_trib(i)", Q_trib(CHANXI(i),CHANYJ(i))
				 print *,"Q_Trib_cms(i)", Q_Trib_cms(CH_NETLNK(CHANXI(i),CHANYJ(i)))
				 print *,"QLateral=", QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i)))
				 print *,"TLateral=", TLateral(i)
				 print *,"% QLINK=", percent_Q_LINK(i)
				 print *,"% Qstr=", percent_Q_strTot(i)
				 print *,"% Qlat=", percent_Q_Lateral(i)
			   print *,"---------------------------"
!				 print *, "QLINK:", QLINK
!			  end if
              ! Deg K to Deg C for the output and energy balance
	      T_QSUM_local(i) = T_QSUM_local(i) - 273.2
!				print *, 'i = ', i
#endif


  end DO

#ifdef HYDRO_D
	print *, "T_QSUM_local = ", T_QSUM_local
	print *, "QLINK = ", QLINK(:,1)
#endif
    !print *,"---------------------------"
    !print *,"QLINK(:,1)="
    !print *,QLINK(:,1)
    !print *,"T_w_sim="
    !print *,T_w_sim
    !print *,"QLateral="
    !print *,QLateral
    !print *,"TLateral="
    !print *,TLateral
    !print *,"T_QSUM_local:"
    !print *,T_QSUM_local

	!increment stuff
	NetRad_prev = NetRad
	T_QSUM_prev = T_QSUM_local
	water_temp_counter = water_temp_counter + 1
  first_time_step_DTCT = .false.


!#endif


end subroutine river_TEMPERATURE_SOLVER



! ------------------------------------------------
!   FUNCTION Shortwave radiation flux
! ------------------------------------------------
real function SW_RAD_FLUX(space, s_in, SF, albedo)
	! Shortwave radiation flux for the desired time step and grid
	implicit none
        ! The node we are calculating the flux for that
	INTEGER, INTENT(IN)                         :: space
        REAL,    INTENT(IN),  DIMENSION(:)  :: s_in, SF, albedo

	! For more information, see Magnusson et al., (2012) and Maidment, (1993)
	SW_RAD_FLUX = s_in(space) * (1.0 - SF(space)) * (1.0 -  albedo(space))

end function SW_RAD_FLUX

! ------------------------------------------------
!   FUNCTION Longwave radiation flux
! ------------------------------------------------

real function LW_RAD_FLUX(space, T_w_sim, T_a, H, sfcp, C_L, viewSky)
	! Longwav radiation flux for the desired time step and grid
	implicit none

	REAL    :: SBConst	! Stefan-Boltzman constant, 5.6696*10-8 (kg·s-3·K-4)
	REAL    :: StrTemp	! simulated water temperature for the back to atmosphere part of the LW flux
	REAL    :: satVapPre	! Saturation Vapor Pressure (mbar)
	REAL    :: VapPreReal  	! Real vapor pressure (mbar)
	REAL    :: emis		! Emissivity of the Atmosphere (unitless)
	REAL    :: LW_Atm	! Atmospheric longwave radiation
	REAL    :: LW_LC	! Landcover radiation
	REAL    :: LW_Back	! Back radiation
	REAL    :: H_tmp		! Temporary value for assigning the humidity from the input
      INTEGER :: space, spaceTemp        ! The node we are calculating the flux for that
        !REAL,    INTENT(IN),  DIMENSION(:,:)  ::  T_w_sim
	REAL,    INTENT(IN),  DIMENSION(:)    ::  T_w_sim
        REAL,    INTENT(IN),  DIMENSION(:)    ::  T_a, H, sfcp, C_L, viewSky
	!INTEGER, INTENT(IN),  DIMENSION(:)    ::  BC_Grid, FROM_NODE

	SBConst = 0.000000056696;	! Initializing the Stefan-Boltzman constant
	StrTemp = T_w_sim(space)

	! Saturation Vapor Pressure (modified from Chapra 1997 and Raudkivi 1979)
	!satVapPre = 0.611 * exp((17.27 * (T_a(space)-273.2)) / (T_a(space)-273.2))
	!!!TLE: See my comment below in latent heat section. This is different in a significant way
	!       from Raudkivi 1979. Changing to that equation (for kPa) until we can get Reza to weigh
	!       in or some other form of closure
	satVapPre = 0.611 * exp((17.27 * (T_a(space)-273.2)) / (237.3 + (T_a(space)-273.2)))

	! Vapor Pressure (Bedient and Huber 1992)
	H_tmp = H(space)
	if (H_tmp .EQ. 100.0) then
			H_tmp = 99.0
	endif
  VapPreReal = (H_tmp / 100.0) * satVapPre
	!!!TLE: Replacing Reza's relative-humidity based calculation with the Noah-MP one.
	!       See latent heat subroutine for additional details
	!VapPreReal = (H(space)*sfcp(space)) / (0.622+(0.378*H(space)))

	! Atmospheric Emissivity (Dingman 2002, Kustas et al. 1994)
	!emis = 1.72 * (((VapPreReal) / (T_a(space)-273.2))**(1.0 / 7.0)) * &
	!	(1.0 + ((0.22 * C_L(space)) ** 2.0))		! See Dingman (2002) p 282
	!!!TLE: Dingman (2002) is a textbook I don't have access to, so can't check that precisely.
	!       Kustas et al. (1994) emissivity equation is as follows:
	!       emis = 0.642 * ((VapPreReal / T_a(space))**(1.0/7.0))
	!       Trying atmospheric emissivity from Boyd and Kasper
	emis = 1.72 * ( ( (0.1*VapPreReal)/(T_a(space)) )**(1.0/7.0) ) * ( 1.22 + (C_L(space)**2.0) )

        ! Put a control here
	if (emis > 0.96) then
    emis = 0.96
  end if
	! Atmospheric Longwave Radiation Flux Attenuated in Water Column
	! (Benyahya et al.,2010; Boyd and Kasper, 2003; Westhoff et al., 2007, 2011)
	LW_Atm = 0.96 * viewSky(space) * emis * SBConst * ((T_a(space))**4.0)

	! Land cover longwave radiation (Boyd and Kasper, 2003; Westhoff et al., 2007)
	LW_LC = 0.96 * (1.0 - viewSky(space)) * 0.96 * SBConst * ((T_a(space)) ** 4.0)

	! Back radiation
	! (Boyd and Kasper, 2003; Evans et al., 1998; Kim and Chapra, 1997; Ouellet et al., 2014; Westhoff et al., 2007)
	LW_Back = -0.96 * SBConst * ((StrTemp + 273.2) ** 4.0)

	! Atmospheric, back, and landcover longwave radiation
	LW_RAD_FLUX = LW_Atm + LW_LC + LW_Back

!print *, "emis = ", emis
!print *, "StrTemp = ", (StrTemp+273.2)
!print *, "T_a(node) = ", T_a(space)
!print *, "LW_Atm = ", LW_Atm
!print *, "LW_LC = ", LW_LC
!print *, "LW_Back = ", LW_Back

end function LW_RAD_FLUX


! ------------------------------------------------
!   FUNCTION Latent heat flux
! ------------------------------------------------
real function LatentHeat_FLUX(space, SW_rad, LW_rad, T_w_sim, T_a, H, sfcp, V_wind, z_elv)

	implicit none
	REAL    :: c_air = 1004.6	! Heat capacity of air (J/kg deg C)
	REAL    :: rho_air = 1.2041	! Density of air at 20 deg C (kg/m^3)
	REAL    :: StrTemp		! simulated water temperature at the previous location
	REAL    :: L_e			! Latent Heat of Vaporization
	REAL    :: satVapPreLat		! Saturation Vapor Pressure (mbar)
	REAL    :: VapPreReal		! Real vapor pressure (mbar)
	REAL    :: s			! Slope of the saturation vapor pressure curve (kPa/C)
	REAL    :: r_a			! Aerodynamic resistance (s/m)
	REAL    :: Pa_air		! Air pressure (kPa)
	REAL    :: psy			! psychrometric constant(kPa/C)
	REAL    :: E			! Penman open water evaporation
	REAL    :: H_tmp		! Temporary value for assigning the humidity from the input
  INTEGER :: space         	! The node we are calculating the flux for that
        REAL    :: SW_rad, LW_rad
        !REAL,    INTENT(IN),  DIMENSION(:,:)  ::  T_w_sim
	REAL,    INTENT(IN),  DIMENSION(:)    ::  T_w_sim
        REAL,    INTENT(IN),  DIMENSION(:)    ::  T_a, H, sfcp, V_wind, z_elv
	!INTEGER, INTENT(IN),  DIMENSION(:)    ::  BC_Grid, FROM_NODE
	StrTemp = T_w_sim(space)

	! <<<<<<<<< Evaporation section >>>>>>>>
	! Calculate the latent heat of vaporization J/kg ((Maidment, 1993)
	!L_e = 1000000 * (2.501 - (0.002361 * StrTemp))
	!L_e = 1000.0 * (2501.4 + StrTemp)  ! this version is from Westhoff (2007) in J/kg
!	print *, "L_e = ", L_e
	!!!TLE: trying Boyd and Kasper version
	L_e = 1000.0 * ( 2501.4 + (1.83 + StrTemp) )

	! Calculate the slope of the saturation vapor pressure curve at a given air temperature
	! (modified from Chapra 1997 and Raudkivi 1979)
	!satVapPreLat = 0.611 * exp((17.27 * T_a(space)-273.2) / (T_a(space)-273.2))

	!!!TLE: As written, the equation above spits out infinity. I have not been able to find this exact equation
	!       (couldn't read Chapra 1997, and not sure where the 'modifications' come in). Trying the
	!       Raudkivi (1979) formulation, which seems to be 611*exp((17.27*(T_a(space)-273.2))/(237.2+(T_a(space)-273.2)))--answer in Pa
	satVapPreLat = 0.611 * exp((17.27 * (T_a(space)-273.2)) / (237.3 + (T_a(space)-273.2))) !0.611 to get answer in kPa

	! Vapor Pressure Based on Saturation Vapor Pressure Using Stream Temperature
	! (Bedient and Huber 1992)
	H_tmp = H(space)
	if (H(space) .EQ. 100.0) then
			H_tmp = 99.0
	endif
  VapPreReal = (H_tmp / 100.0) * satVapPreLat
  !!!TLE: This method for calculating vapor pressure doesn't work, because H is specific humidity,
	!       not relative humidity as Reza assumed. Replacing the above calculation with the method used
	!       in Noah-MP, which calculates vapor pressure based on surface pressure and spec. humidity.
	!VapPreReal = (H(space)*sfcp(space)) / (0.622+(0.378*H(space)))

	! Slope of the saturation vapor pressure curve (kPa/C) is calculated as
	! (Maidment, 1993)
	!s = (4100.0 * satVapPreLat) / ((T_a(space)-273.2) ** 2.0)
	!!!TLE: Finding a slightly different equation...we'll try the one I found.
	!s = (4100.0 * satVapPreLat) / ((237.3 + (T_a(space)-273.2)) ** 2.0)
	!!!TLE: re-writing to try the pure formulation from the Oregon manual
	s = satVapPreLat - 0.611 * exp((17.27 * (T_a(space)-273.2-1.0)) / (237.3 + (T_a(space)-273.2-1.0)))

	! Calculate the aerodynamic resistance (s/m)
	! (Westhoff et al., 2007)
	!r_a = 245.0 / ((0.54 * V_wind(space)) + 0.5)
	!!!TLE: repurposing this variable to Boyd and Kasper's aerodynamic evaporation
	!!!    For their wind function (fW = a+b*V_wind), assuming a = 3e-9 and b = 1e-9
  r_a = (0.000000003 + 0.000000001 * V_wind(space)) * (satVapPreLat - VapPreReal)

	! Calculate air pressure
	!Pa_air = 101.3 - (0.01055 * z_elv(space))
	!!!TLE: this also comes from Westhoff et al. (2007), but now we're passing in surface pressure so just use that
	Pa_air = sfcp(space)

	! Calculate psychrometric constant(kPa / deg C) (based on air pressure(Pa),
	! (value should be adjusted for different site elevations)
	! Ratio of water to dry air = .622, and the latent heat of water vaporization = 2.45E6 (J/kgC)
	! (Cuenca 1989)
	!psy = (c_air * Pa_air) / (0.622 * 2450000.0)
	!!!TLE: trying Boyd and Kasper (basically the same, but with our own calculated latent heat)
	psy = (c_air * Pa_air) / (0.622 * L_e)

	! Calculate the Penman open water evaporation
	! (Maidment, 1993; Westhoff et al., 2007)
	!E = ((s * (SW_rad + LW_rad)) / (1000.0 * L_e * (s + psy))) &
	!	+ ((c_air * rho_air * psy * (satVapPreLat - VapPreReal)) / (1000.0 * L_e * r_a * (s + psy)))
	!!!TLE: This is slightly different from what's in Westhoff et al. (2007). There's an extra psychrometric constant in here.
	!E = ((s * (SW_rad + LW_rad)) / (1000.0 * L_e * (s + psy))) + &
	!((c_air * rho_air * (satVapPreLat - VapPreReal)) / (rho_water * L_e * r_a * (s + psy)))
	!!!TLE: trying the Boyd and Kasper version
	E = (((s * (SW_rad + LW_rad)) / (rho_water * L_e))+( r_a * psy )) / (s + psy)

	LatentHeat_FLUX = -rho_water * L_e * E

!  print *, "-----LATENTHEAT_FLUX-----"
  !print *, "psy = ", psy
	!print *, "L_e = ", L_e
	!print *, "E = ", E
!  print *, "s = ", s
!	print *, "SW_rad = ", SW_rad
!	print *, "LW_rad = ", LW_rad
!	print *, "L_e = ", L_e
!	print *, "psy = ", psy
!	print *, "c_air = ", c_air
!	print *, "rho_air = ", rho_air
!	print *, "satVapPreLat = ", satVapPreLat
!	print *, "VapPreReal = ", VapPreReal
!	print *, "r_a = ", r_a
!	print *, "-----end-----"

end function LatentHeat_FLUX


! ------------------------------------------------
!   FUNCTION Sensible heat flux
! ------------------------------------------------
real function SensibleHeat_FLUX(space, T_w_sim, LH, T_a, H, sfcp, z_elv)

	implicit none
	REAL    :: c_air = 1004.6	! Heat capacity of air (J/kg deg C)
	REAL    :: rho_air = 1.2041	! Density of air at 20 deg C (kg/m3)
	REAL    :: ews			! Saturation vapor pressure of the evaporating surface (KPa)
	REAL    :: ewa			! Actual vapor pressure of the evaporation surface (KPa)
	REAL    :: StrTemp		! Simulated water temperature at the previous location
	REAL    :: Pa_air		! Air pressure (mbar)
	REAL    :: H_tmp		! Temporary value for assigning the humidity from the input
	REAL    :: B_r			! Bowen's ratio
	REAL    :: psy      ! psychrometric constant
	INTEGER :: space
        REAL    :: LH
        !REAL,    INTENT(IN),  DIMENSION(:,:)  ::  T_w_sim
	REAL,    INTENT(IN),  DIMENSION(:)    ::  T_w_sim
        REAL,    INTENT(IN),  DIMENSION(:)    ::  T_a, H, sfcp, z_elv
	!INTEGER, INTENT(IN),  DIMENSION(:)    ::  BC_Grid, FROM_NODE

	StrTemp = T_w_sim(space)

	! Saturation vapor pressure of the evaporating surface (KPa)
	! (Dingman, 1994; Maidment, 1993)
	!ews = 0.61275 * exp((17.27 * StrTemp) / (StrTemp))
	!!!TLE: The above is a slightly different equation from what is elsewhere.
	!      Putting in the Raudkivi formulation for now while I work on aligning
	!      the pressure calculations with Noah-MP's version
	ews = 0.611 * exp((17.27 * StrTemp)/(237.3 + StrTemp))

	! Actual vapor pressure of the evaporation surface (KPa)
	! To handle the infinte problem in B_r equation, consider the umidity of 100 as 99
	H_tmp = H(space)
        if (H(space) .EQ. 100.0) then
            H_tmp = 99.0
        endif
	ewa = (H_tmp / 100.0) * ews

	!!!TLE: replacing vapor pressure formulation. See latent heat subroutine for more info.
	!ewa = (H(space)*sfcp(space)) / (0.622+(0.378*H(space)))

	! Adiabatic atmospheric pressure adjusted for elevation above sea level
	! (Westhoff et al., 2007)
	!Pa_air = 101.3 * (((293 - (0.0065 * z_elv(space))) / 293)**5.256)
	!!!TLE: The above is NOT the equation in Westhoff et al. (2007). I don't know where it
	!       came from, but going to change it to be consistent with the air pressure
	!       calculation in Westhoff, which is also what is used in the psychrometric
	!       constant calculation in the latent heat module
	!Pa_air = 101.3 - (0.01055 * z_elv(space))
  !!!TLE: Update 2. We are now reading in surface pressure, so don't need to calculate
	Pa_air = sfcp(space)

	!psy = (c_air * Pa_air) / (0.622 * 2450000.0)

	!!! For now, switching back to Reza's method (relative humidity based) with RH hard coded at 30%
	!!! TLE temporary hack
	! as formulated, sometimes ewa is greater than ews, which messes up the sign on the Bowen ratio
	! adding an if statement to ensure ewa is always less than ews
	!if (ewa .ge. ews) then
	!    ewa = ews * 0.95
	!end if

	! Bowen ratio
	! (Evans et al.,1998; Magnusson et al., 2012; Webb and Zhang, 1997; Westhoff et al., 2007)
	! (modified for units from Bedient and Huber 1992)
	B_r = 0.00061 * Pa_air * ((StrTemp - (T_a(space)-273.2)) / (ews - ewa))

	!!!!TLE temporary hack
	!! seeing bowen ratio values as low as -25. This seems very high, so just testing a temporary lower threshold
	!if (B_r .lt. -1.) then
	!    B_r = -1.
	!end if

	!!!TLE: modified B_r to use the psychrometric constant we calculate
	!B_r = psy * ((StrTemp - (T_a(space)-273.2)) / (ews - ewa))

	SensibleHeat_FLUX = B_r * LH
!	print *, "B_r = ", B_r
!  print *, "Pa_air = ", Pa_air
!	print *, "ews = ", ews
!	print *, "ewa = ", ewa

end function SensibleHeat_FLUX


! ------------------------------------------------
!   FUNCTION Conduction Flux
! ------------------------------------------------
real function BedSediment_FLUX(space, dx, T_w_sim, P_s, C_E, Pw_river, T_bed)

	implicit none

	REAL  :: StrTemp		! Simulated water temperature at the previous location
	REAL  :: Dummy1			!
	REAL  :: Dummy2			!
	REAL  :: P			! Estimated porosity
	REAL  :: ThermalDiffuse		! Volumetric Weighted Thermal Diffusivity (m2/s)
	REAL  :: Sed_ThermalDiffuse	! Substrate Thermal Diffusivity (m2/s)
	REAL  :: H2O_ThermalDiffuse	! Water Thermal Diffusivity (m2/s)
	REAL  :: Ratio_Sediment		! To save the value of Volume_Sediment / Volume_cl
	REAL  :: Sed_Depth		! Sediment depth used in bed conduction (m)
	REAL  :: Volume_Sediment	! Substrate Volume in Conduction Layer (m3)
	REAL  :: Volume_cl		! Total Volume in Conduction Layer (m3)
	REAL  :: Volume_H2O		! Water Volume in Conduction Layer (m3)
	REAL  :: Density		! Volumetric Weighted Density (kg/m3)
	REAL  :: Sed_Density		! Substrate Density (kg/m3)
	REAL  :: HeatCapacity		! Volumetric Weighted Specific Heat Capacity (J/kgC)
	REAL  :: Sed_HeatCapacity	! Substrate Specific Heat Capacity (J/kgC)
        REAL  :: Ratio_H2O, dx
        INTEGER :: space
        !REAL,    INTENT(IN),  DIMENSION(:,:)  ::  T_w_sim
	REAL,    INTENT(IN),  DIMENSION(:)    ::  T_w_sim
        REAL,    INTENT(IN),  DIMENSION(:)    ::  P_s, C_E, Pw_river, T_bed
	!INTEGER, INTENT(IN),  DIMENSION(:)    ::  BC_Grid, FROM_NODE
	StrTemp = T_w_sim(space)
	! Assigning the numbers
	!!!TLE: from where?!?!
	Sed_ThermalDiffuse = 0.0000045
	H2O_ThermalDiffuse = 0.00000014331 !this looks right (in m2/s)
	Sed_Density = 1600.0      ! google agrees with this
	Sed_HeatCapacity = 2219.0 !from Google, spec. heat capacity of sand is ~2900 J/kgC

	! ======================================================
	! Calculate Volumetric Ratio of Water and Substrate
	! Code uses this ratio to estimate conduction constants
	! First, calculating the porosity
	! Ratio Size of dominant substrate
	Dummy1 = P_s(space) * (1.0 - C_E(space))
	!!!TLE: from CHANPARM: P_s = 50 mm; C_E = 0.5
	!      So, Dummy1 always = 25 mm

	! Ratio Conductivity of sand - low range
	Dummy2 = 0.062 * C_E(space)
	!!!TLE: Dummy2 always = 0.031 [units?]

	! Estimated Porosity
	! (adopted from Bedient and Huber 1992)
	P = (0.3683 * ((Dummy1 + Dummy2)**(-0.0641)))
	print *, "Porosity = ", P
	!!!TLE: haven't been able to locate this equation...
	!      According to this and the above,
	!      P always = 0.2994
	!      This seems ~reasonable, depending on substrate

	! ======================================================
	! Variables used in bed conduction
	! Calculate the sediment depth (conduction layer)
	! (Assumes 100 Particles Thickness)
	! 0.1m < particleSize < 0.3m
	!!!TLE: confused by the above, since particle size is set to 50 mm = 0.05 m
	Sed_Depth = 10.0 * P_s(space) / 1000.0;
	!!!TLE: so, sed_depth always = 0.5
	if (Sed_Depth > 1.0) then
		Sed_Depth = 1.0;
	else if (Sed_Depth < 0.1) then
		Sed_Depth = 0.1
	end if
	! ======================================================

	Volume_Sediment = (1.0 - P) * Pw_river(space) * Sed_Depth * dx
	Volume_cl = Pw_river(space) * Sed_Depth * dx

	if (Volume_cl .EQ. 0.) then
	   Ratio_Sediment = 0.
	else
	   Ratio_Sediment = Volume_Sediment / Volume_cl
	end if

	Volume_H2O = P * Pw_river(space) * Sed_Depth * dx

	if (Volume_cl .EQ. 0.) then
	   Ratio_H2O = 0.
	else
	   Ratio_H2O = Volume_H2O / Volume_cl
	end if

	!!!TLE: As written, porosity is constant at ~0.3, so we should have a
	!      ~70:30 ratio of sediment to water
	print *, "Ratio_H2O = ", Ratio_H2O
	print *, "Ratio_Sediment = ", Ratio_Sediment

	ThermalDiffuse = (Sed_ThermalDiffuse * Ratio_Sediment) + (H2O_ThermalDiffuse * Ratio_H2O)
	Density = (Sed_Density * Ratio_Sediment) + (rho_water * Ratio_H2O)
	HeatCapacity = (Sed_HeatCapacity * Ratio_Sediment) + (c_water * Ratio_H2O)

	BedSediment_FLUX = ThermalDiffuse * Density * HeatCapacity * ((T_bed(space)-273.2) - StrTemp) / (Sed_Depth / 2.0)

	print *, "-----BedSediment_FLUX-----"
	print *, "c_water = ", c_water
	print *, "ThermalDiffuse = ", ThermalDiffuse
	print *, "Density = ", Density
	print *, "HeatCapacity = ", HeatCapacity
	print *, "Sed_Depth = ", Sed_Depth
	print *, "(T_bed(space)-273.2) - StrTemp = ", (T_bed(space)-273.2) - StrTemp
	print *, "----------"

end function BedSediment_FLUX

end module module_river_temperature
