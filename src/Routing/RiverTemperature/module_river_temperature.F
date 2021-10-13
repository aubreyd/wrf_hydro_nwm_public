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
	    
	    
subroutine river_temperature_init(this, T_QSUM)

	implicit none
	class(river_temperature_type), INTENT(inout)       :: this              ! the type object being initialize
        REAL, INTENT(inout), DIMENSION(:)             :: T_QSUM        
		
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
Subroutine river_TEMPERATURE_SOLVER(this, dx, DT, DTCT, IXRT, JXRT,  T_a, H, V_wind, C_L, T_bed, s_in, Pw_river, y_w, S_0, z_elv, V_w, W_w, Q_str, Q_GW, &
				Q_Trib, T_GW, T_hyp, T_Trib, albedo, SF, k_s_0, C_E, P_s, CHANXI, CHANYJ, CH_NETLNK, &
				CH_NETRT, QLINK, node_area, TO_NODE, TYPEL, T_QSUM_local, totLength)
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
	REAL, INTENT(IN), DIMENSION(:)    :: H		        ! Humidity (%) - relative humidity			                       - As an input
	REAL, INTENT(IN), DIMENSION(:)    :: V_wind	        ! Wind speed (m/s)							       - As an input
	REAL, INTENT(IN), DIMENSION(:)	  :: C_L	        ! Cloudiness (0 to 1)					                       - As an input - should find it from somewhere!
	REAL, INTENT(IN), DIMENSION(:)	  :: T_bed	        ! Conduction Layer Temperature (C)			                       - As an input
	REAL, INTENT(IN), DIMENSION(:)    :: s_in		! Incoming sortwave radiation, as an input
	REAL, INTENT(IN), DIMENSION(:)    :: Pw_river	        ! Wetted perimeter in the river channle for the desird grid (m)                 =====> XS_Peremeter
	REAL, INTENT(IN), DIMENSION(:)    :: y_w		! Water Column Depth (m)                                    	                =====> HLINK
	REAL, INTENT(IN), DIMENSION(:)	  :: S_0		! Stream gradient (unitless)							=====> ChSSlp
        REAL, INTENT(IN), DIMENSION(:)	  :: z_elv		! Elevation of station where met. data was obtained				=====> ZELEV
	REAL, INTENT(IN), DIMENSION(:)	  :: V_w		! Average Flow Velocity (m/s)							=====> velocity
	REAL, INTENT(IN), DIMENSION(:)	  :: W_w		! Wetted width (Top; m)								=====> Tw											
	! Discharge data, instream, groundwatr flow, triburaty
	REAL, INTENT(IN), DIMENSION(:)	  :: Q_str		! Instream flow in the river channel - from sll sources (csm)                   =====> QSUM
	REAL, INTENT(IN), DIMENSION(:,:)  :: Q_GW		! Groundwater flow in the gridded based (cms)					=====> QINFLOWBASE
	REAL, INTENT(IN), DIMENSION(:,:)  :: Q_Trib		! Tributarty surface inflow in the gridded based (cms)			        =====> QSTRMVOLRT
	REAL, INTENT(IN), DIMENSION(:)	  :: T_GW		! Groundwater flow temperature in the gridded based (C)		         	- As an input
	REAL, INTENT(IN), DIMENSION(:)	  :: T_hyp		! Hyporheic Flow temperature in the gridded based (C)		         	- As an input
	REAL, INTENT(IN), DIMENSION(:)	  :: T_Trib		! Tributarty surface inflow temperature (C)					- As an input
	! Shading related parameters
	REAL, INTENT(IN), DIMENSION(:)    :: albedo		! albedo at the stream surface (0 to 1)						- ???
	REAL, INTENT(IN), DIMENSION(:)    :: SF			! Shading factor for each grid of the reach					- ???
	! Fluxes
	REAL                              :: SW_rad		! Shortwave radiation flux (W/m2)
	REAL                              :: LW_rad		! Longwave radiation flux (W/m2)
	REAL                              :: LH		        ! Latent heat energy flux (W/m2)
	REAL                              :: SH		        ! Sensible heat energy flux (W/m2)
	REAL                              :: CH		        ! Conduction flux (W/m2)
	REAL,  DIMENSION(totLength)       :: Net_HeatFlux	! Net heat flux from 5 sources (W/m2)			
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
	INTEGER,                DIMENSION(totLength)          :: BC_Grid		    ! -- This is the new way I am thinking of to detec the upstream grids
	INTEGER                                               :: neighbours, i, j
        ! Simulated water temeprature (C)
	REAL,                   DIMENSION(totLength)          :: T_w_sim	! Simulatd water temperature for the reach in space & time (C)
        REAL                                                  :: Q_strTot, T_strTot 
	REAL,                   DIMENSION(totLength)          :: viewSky 		! View to sky factor (1-SF)
	REAL,                   DIMENSION(totLength)          :: T_w_sim_temp	
						      
        !REAL                                                  :: DT_STEPS               !-- number of timestep in routing
	REAL,    allocatable,   DIMENSION(:,:)                :: Q_GW_csm, Q_Trib_csm
	logical :: first_time_step_DT = .true.
	logical :: first_time_step_DTCT = .true.
	
	viewSky = 1. - SF
        !DT_STEPS = INT(DT/DTCT)   !-- fix the timestep
        ! ++++++++++++++++++++++++++++update this part
	! Assigning the upstream water temeprature boundary condition to the main matrix


	Disp_coef = 0.
	!-- Setting everything to zero first.
        Q_hyp = 0.
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

	       ! For the last grid where there is no water surface in the downstream
	       if (CH_NETLNK(CHANXI(i),CHANYJ(i)) .EQ. totLength) then
			!h_D = y_w(CH_NETLNK(CHANXI(i),CHANYJ(i)))
			d_h = 0
	       else
	       h_D = y_w(CH_NETLNK(CHANXI(i),CHANYJ(i)))
	       !-- Here also I have the issue with the directiono of the cell. I am using TO_NODE here to 
	       !-- detect the cell in the downstream.
	       !h_D_dn = y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))+1) - (dx * S_0(CH_NETLNK(CHANXI(i),CHANYJ(i))))
	       h_D_dn = y_w(TO_NODE(i)) - (dx * S_0(CH_NETLNK(CHANXI(i),CHANYJ(i))))

	       !Change in hydraulic head per unit distance
	       d_h = (h_D - h_D_dn) / dx
	       end if

		!h_D = y_w(CH_NETLNK(CHANXI(i),CHANYJ(i)))

		! Cross-Sectional Area of the Seepage Face (m2)
		A_s = dx * S_0(CH_NETLNK(CHANXI(i),CHANYJ(i))) * Pw_river(CH_NETLNK(CHANXI(i),CHANYJ(i)))
		! Darcy’s Law for Calculating for Hyporheic Flow (cms 
		Q_hyp(CHANXI(i),CHANYJ(i)) = A_s * k_s * d_h 

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
		! I Changed 0.011 to 0.06 in below eqn to account for natural (not straight) streams
			
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
	        
  
	!-- I have this long if statement for determining the headwaters within the grid. 
	BC_Grid = 0
	!-- WITH NEW UPDATE AND LOOKING DOWNSTREAM, WE REALLY DON'T NEED THIS IF AND UPDATING BC_GRID,
	!-- BUT I KEEP THAT HERE FOR NOW! 
	if (first_time_step_DT) then	!-- To do this only one time!
		neighbours = 0
		DO i = 1,IXRT  !rows
			DO j = 1,JXRT   !colsumns
				if ( CH_NETRT(i,j) .EQ. 0) then
				if ((i .gt. 1) .and. (i .lt. IXRT) .and. (j .gt. 1) .and. (j .lt. JXRT)) then
					if (CH_NETRT(i-1,j-1) .NE. -9999) then
                  				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i-1,j) .NE. -9999) then 
                  				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i-1,j+1) .NE. -9999) then 
                  				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j-1) .NE. -9999) then
                  				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j+1) .NE. -9999) then
                  				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j-1) .NE. -9999) then
                  				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i+1,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1	
             	 			end if
				else if ((i .ne. 1) .and. (i .ne. IXRT) .and. (j .eq. 1)) then 
					if (CH_NETRT(i-1,j) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i-1,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i+1,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if	
				else if ((i .eq. IXRT) .and. (j .ne. 1) .and. (j .ne. JXRT)) then
					if (CH_NETRT(i,j-1) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i-1,j-1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i-1,j) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i-1,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
				else if ((j .ne. 1) .and. (j.ne. JXRT) .and. (i .eq. 1)) then
					if (CH_NETRT(i,j-1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j+1) .NE. -9999) then 
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
				else if ((i .ne. 1) .and. (i .ne. IXRT) .and. (j .eq. JXRT)) then
					if (CH_NETRT(i-1,j) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i-1,j-1) .NE. -9999) then
				                 neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j-1) .NE. -9999) then
				                 neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i+1,j) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j) .NE. -9999) then 
                 				neighbours = neighbours + 1
              				end if
				else if ((i .eq. 1) .and. (j .eq. 1)) then
					if (CH_NETRT(i,j+1) .NE. -9999) then
						neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i+1,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if			
				else if ((i .eq. IXRT) .and. (j .eq. 1)) then
					if (CH_NETRT(i-1,j) .NE. -9999) then
                 				neighbours = neighbours + 1	
              				end if
					if (CH_NETRT(i-1,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j+1) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
				else if ((i .eq. 1) .and. (j .eq. JXRT)) then
					if (CH_NETRT(i,j-1) .NE. -9999) then
                 				neighbours = neighbours + 1
				        end if
	      				if (CH_NETRT(i+1,j-1) .NE. -9999) then 
                 				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i+1,j) .NE. -9999) then
                 				neighbours = neighbours + 1
              				end if
				else if ((i .eq. IXRT) .and. (j .eq. JXRT)) then
					if (CH_NETRT(i-1,j-1) .NE. -9999) then
                				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i-1,j) .NE. -9999) then
                				neighbours = neighbours + 1
              				end if
					if (CH_NETRT(i,j-1) .NE. -9999) then
                				neighbours = neighbours + 1	
              				end if
				end if

				if ((neighbours .eq. 1) .and. (CH_NETLNK(i,j) .ne. totLength)) then
					BC_Grid(CH_NETLNK(i,j)) = 1
				end if

				neighbours = 0
          			end if
			end do
		end do
		!-- to convert to Deg C 
		T_w_sim_temp = T_Trib - 273.2
		first_time_step_DT = .false.
      
      !-- Preparing the T_w_sim to be used in the first time step directly from the inputs from NoahMP
      !T_w_sim = 0.
      !print *,"T_trib:"
      !print *,T_Trib
      !T_w_sim(:,1) =  T_Trib - 273.2
      !T_w_sim_temp = T_Trib - 273.2
      
      
      else
	 !-- It's supposed to be in Deg C
         T_w_sim_temp = T_QSUM_local
 
      end if
     
      !-- Converting two inflows to cms
      Q_Trib_csm = Q_Trib  / 1000 * node_area(CH_NETLNK(CHANXI(i),CHANYJ(i)))/DT
      Q_GW_csm = Q_GW  / 1000 * node_area(CH_NETLNK(CHANXI(i),CHANYJ(i)))/DT

      	
      ! Caclulating the temperature
      Do i = 1 , totLength
		! Calculate the heat fluxes
		SW_rad = SW_RAD_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), s_in, SF, albedo)
		LW_rad = LW_RAD_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), T_w_sim_temp, T_a, H, C_L, viewSky)	
		LH     = LatentHeat_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), SW_rad, LW_rad, T_w_sim_temp, T_a, H, V_wind, z_elv)
		SH     = SensibleHeat_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)), T_w_sim_temp, LH, T_a, H, z_elv)
		CH     = BedSediment_FLUX(CH_NETLNK(CHANXI(i),CHANYJ(i)),  dx, T_w_sim_temp, P_s, C_E, Pw_river, T_bed)
		! Net heat flux
		Net_HeatFlux(CH_NETLNK(CHANXI(i),CHANYJ(i))) = SW_rad + LW_rad + LH + SH + CH
		! Calculate mixture temperature for each node and time step 
		if (y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) .EQ. 0.) then
		   Advection = 0.
		else
		   Advection = Net_HeatFlux(CH_NETLNK(CHANXI(i),CHANYJ(i))) / (rho_water * c_water * &
			    y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))
		end if

		! For the first time step, where the temperature is not calculated for that time step yet:
		! This is impprtant though! As, basically, here I am definig the upstream boundary condition.
		! -- I think there is no need for the next "if" after my my new edits right after thif "if"
		! -- I am gonna turn this off
		! -- ##### FROM HERE #####
		!if (BC_Grid(CH_NETLNK(CHANXI(i),CHANYJ(i))) .eq. 1) then ! -- somethoing like this if statement.
		!   if ((Q_Trib_csm(CHANXI(i),CHANYJ(i))+Q_GW_csm(CHANXI(i),CHANYJ(i))+ Q_hyp(CHANXI(i),CHANYJ(i))) .EQ. 0.) then 
		!      T_w_sim(CH_NETLNK(CHANXI(i),CHANYJ(i)),2) = 0
		!   else
		!      T_w_sim(CH_NETLNK(CHANXI(i),CHANYJ(i)),2) =(Q_Trib_csm(CHANXI(i),CHANYJ(i))*(T_Trib(i)-273.1)+ &                 
		!					 Q_GW_csm(CHANXI(i),CHANYJ(i))*(T_GW(i)-273.1)+ &
		!					 Q_hyp(CHANXI(i),CHANYJ(i))*(T_hyp(i)-273.1))/ &
		!					(Q_Trib_csm(CHANXI(i),CHANYJ(i))+Q_GW_csm(CHANXI(i),CHANYJ(i))+ &
		!					 Q_hyp(CHANXI(i),CHANYJ(i)))
		!   end if
		!end if
		! -- ##### TO HERE #####

		! For the very upstream grid where there is no information in the prior grid:
		! Since I am looking downstream inplace of upstream now!
		!if (BC_Grid(CH_NETLNK(CHANXI(i),CHANYJ(i))) .eq. 1) then
		!	self = 0.
		!	if ((y_w(CH_NETLNK(CHANXI(i),CHANYJ(i)))) .EQ. 0.) .OR. (V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) .EQ. 0.)) then
		!	   Dispersion = 0.
 		!	   !FluxTemp = T_w_sim(CH_NETLNK(CHANXI(i),CHANYJ(i)),1)
		!	   FluxTemp = T_w_sim_temp(CH_NETLNK(CHANXI(i),CHANYJ(i)))
		!	else
		!	   Dispersion = (-1.0 * Disp_coef(CHANXI(i),CHANYJ(i)) * ((Net_HeatFlux(CH_NETLNK(CHANXI(i),CHANYJ(i))) / &
		!				 (rho_water * c_water * y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) * &
		!				  V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))) - self)) / dx
		!	   !FluxTemp = T_w_sim(CH_NETLNK(CHANXI(i),CHANYJ(i)),1) + (Advection + Dispersion) * (dx/V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))
		!	   FluxTemp = T_w_sim_temp(CH_NETLNK(CHANXI(i),CHANYJ(i))) + (Advection + Dispersion) * (dx/V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))
		!	end if

		!-- Sooo, I had to update this code here because I don't know how to get access to the upper cell which is discharing to our target cell.
		!-- I though the "FROM_NODE" would be the key but after printing the array, it turns out that the array is actually starts form 1 to the NLINKS.
		!-- nothing speccific about the FROM_NODE
		!-- To me, that means the model doesn't care about that FROM_NODE. 
		!-- Therefore, I updated the model's formulation to wotk with TO_NODE instead 
		!else
		!-- Self
                 if ((y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) .EQ. 0.) .OR. (V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) .EQ. 0.)) then
		     self = 0.
		 else
		     self = Net_HeatFlux(CH_NETLNK(CHANXI(i),CHANYJ(i))) / &
				(rho_water * c_water * y_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) * V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))
		 end if
		 !-- looking downstream
		 if ((y_w(TO_NODE(i)) .EQ. 0.) .OR. (V_w(TO_NODE(i)) .EQ. 0.)) then
		     Dispersion = (-1.0 * Disp_coef(CHANXI(i),CHANYJ(i)) * -1. * self) / dx  
		 else
		     Dispersion = (-1.0 * Disp_coef(CHANXI(i),CHANYJ(i)) * ((Net_HeatFlux(TO_NODE(i)) / &
				 (rho_water * c_water * y_w(TO_NODE(i)) * V_w(TO_NODE(i)))) - self)) / dx
		 end if
		 !-- Putting everything together, for the advection and dispersion
		 if (V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))) .EQ. 0.) then 
			   !FluxTemp = T_w_sim(CH_NETLNK(CHANXI(i),CHANYJ(i))-1,1) + (Advection + Dispersion) * (dx/V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))
			   FluxTemp = T_w_sim_temp(CH_NETLNK(CHANXI(i),CHANYJ(i)))
		 else
			   FluxTemp = T_w_sim_temp(CH_NETLNK(CHANXI(i),CHANYJ(i))) + (Advection + Dispersion) * (dx/V_w(CH_NETLNK(CHANXI(i),CHANYJ(i))))
		 end if

		!end if

		! -- We don't apply the lateral inflows to T_w_sim; They will be applied to the final porduct, T_QSUM later, 
                !-- Deg C to Deg K
		T_w_sim(CH_NETLNK(CHANXI(i),CHANYJ(i))) = FluxTemp + 273.2

		!-- getting a weighted average
		if (CH_NETLNK(CHANXI(i),CHANYJ(i)) .gt. 0) then
		   if ((Q_Trib_csm(CHANXI(i),CHANYJ(i))+Q_GW_csm(CHANXI(i),CHANYJ(i))+ Q_hyp(CHANXI(i),CHANYJ(i))) .EQ. 0.) then
		      TLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = 0.
		   else
		      !-- in Deg K
		      TLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = ((Q_Trib_csm(CHANXI(i),CHANYJ(i))*T_Trib(i)) + &                                                     
                                                            (Q_GW_csm(CHANXI(i),CHANYJ(i))*T_GW(i)) + &
							    (Q_hyp(CHANXI(i),CHANYJ(i))*T_hyp(i))) / &
							   (Q_Trib_csm(CHANXI(i),CHANYJ(i))+Q_GW_csm(CHANXI(i),CHANYJ(i))+ &
							    Q_hyp(CHANXI(i),CHANYJ(i)))
		   end if

		!-- For now, I realized that we're dealing with 4 different discharge values: QSUM, QINFLOWBASE, QSTORMVOLRT, QLINK
		!-- I am getting these values from module_channel_routing. Overall, it seems they all are in m3/s
		!-- Okay, let's look at the units of these variables based on the code we have in module_channel_routing:
		!-- L779:  QINFLOWBASE, QSTORMVOLRT: Based on the conversion here (Q/1000*node_area/Time), in code mm is going to be m3/s (I converted them before the DO)
		!-- L1089: QLINK: Well, this one is based on the diffusion function and as we can see from L181 - L187,
		!--  the diffusion function uses Manning's equation for calculating the answer. Therefore, it should also be on m3/s					    
		!-- L1110: QSUM: Since it's a function of QLINK with no additional conversion, we can say, that one also is in m3/s			    
                !-- Okay, based on what I learned from module_channel_routing, I am gonna remove the DT_STEPS from caluating the QLateral
		
	            QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))) = Q_Trib_csm(CHANXI(i),CHANYJ(i)) + Q_GW_csm(CHANXI(i),CHANYJ(i)) + Q_hyp(CHANXI(i),CHANYJ(i))										  
	        end if

                              
	      !-- Here I wanna detect the cells/grids discharging to the current cell
	      !-- Not very efficient though
	      Q_strTot = 0.
	      T_strTot = 0.
	      Do j = 1, totLength
	      	 if (TO_NODE(j) .eq. i) then
	      	     Q_strTot = Q_strTot + Q_str(j)
	      	     if ((Q_strTot +  Q_str(j)) .EQ. 0.) then
	      		T_strTot = 0.
		     else
			T_strTot = ((Q_strTot * T_strTot)+(Q_str(j)*T_w_sim_temp(j))) / &
				(Q_strTot + Q_str(j))		
		     end if   
		  end if
		  
	      end do
	      
              !-- to Deg K
	      T_strTot = T_strTot + 273.2
			
	      !-- inflow to node across each face  
	      !if(TYPEL(i) .eq. 0) then                         ! -- If stream; only regular nodes have to attribute
	      !-- in Deg K
	      T_QSUM_local(i) = ((ABS(Q_strTot) * T_strTot) + &
			   (ABS(QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1))*T_w_sim(i)) + &
			   (ABS(QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i))))*TLateral(i))) / &
			   (ABS(Q_strTot)+ABS(QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1))+ABS(QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i)))))

			   !print *,"---------------------------"
			   !print *,"Q_strTot=", Q_strTot
			   !print *,"T_strTot=", T_strTot
			   !print *,"QLINK(i,1)=", QLINK(CH_NETLNK(CHANXI(i),CHANYJ(i)),1)
			   !print *,"T_w_sim(i,2)=", T_w_sim(i)  
			   !print *,"QLateral=", QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i)))
			   !print *,"TLateral=", TLateral(i)   
			   !print *,"T_QSUM_local(i)", T_QSUM_local(i)
			   !print *,"---------------------------"
              
              ! Deg K to Deg C for the output and energy balance
	      T_QSUM_local(i) = T_QSUM_local(i) - 273.2

              ! Update the vector to be used in the next time step.		 
	      !T_w_sim(i) = T_QSUM_local(i)						

	      ! -- Q_str = QSUM(TO_NODE(i))                   : The amount of the water is gonna be sent to the grid in the downstream?
	      ! -- 						QSUM's unit is volume, that means the QLateral also should be converted to the volume. 
	      ! -- T_QSUM(TO_NODE(i))                         : The temperature we are gonna attach to the QSUM to affect it's downstream
	      ! -- QLINK(i,1)                                 : Volume of the water that updates the QSUM and in the case of cells 9 & 13 go to 14 (my example),  
	      !                                                 They could both update the flux of the cell #14. It's unit also is volume
	      ! -- T_w_sim(i)                                 : Simulated water temperature based on the energy budget calculations, not affected by the lateral inlflows
	      ! -- QLateral(CH_NETLNK(CHANXI(i),CHANYJ(i)))   : Lateral surface and subsurface inflows
	      ! -- TLateral(i)                                : Temperature of the lateral surface and subsurface inflows,(combined) 
	      !						        based on the equation I showed earlier today

	      ! For now, the idea is to run the temperature module based on the time step we have for the hydraulics routing module but report the results only for every hour. 

     end DO

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

    first_time_step_DTCT = .false.
             
    !-- Updating the T_w_sim for the next timestep
    !T_w_sim(:,1) = T_w_sim(:,2)

		
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

real function LW_RAD_FLUX(space, T_w_sim, T_a, H, C_L, viewSky)
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
        INTEGER :: space, spaceTemp        ! The node we are calculating the flux for that
        !REAL,    INTENT(IN),  DIMENSION(:,:)  ::  T_w_sim
	REAL,    INTENT(IN),  DIMENSION(:)    ::  T_w_sim
        REAL,    INTENT(IN),  DIMENSION(:)    ::  T_a, H, C_L, viewSky
	!INTEGER, INTENT(IN),  DIMENSION(:)    ::  BC_Grid, FROM_NODE
        
	SBConst = 0.000000056696;	! Initializing the Stefan-Boltzman constant
	StrTemp = T_w_sim(space)
	
	! Saturation Vapor Pressure (modified from Chapra 1997 and Raudkivi 1979)
	satVapPre = 0.611 * exp((17.27 * (T_a(space)-273.2)) / (T_a(space)-273.2))
	! Vapor Pressure (Bedient and Huber 1992)
	VapPreReal = (H(space) / 100.0) * satVapPre
	! Atmospheric Emissivity (Dingman 2002, Kustas et al. 1994)
	emis = 1.72 * (((VapPreReal) / (T_a(space)-273.2))**(1.0 / 7.0)) * &
		(1.0 + ((0.22 * C_L(space)) ** 2.0))		! See Dingman (2002) p 282

        ! Put a control here		
	if (emis > 0.96) then
            emis = 0.96
        end if
	
	! Atmospheric Longwave Radiation Flux Attenuated in Water Column 
	! (Benyahya et al.,2010; Boyd and Kasper, 2003; Westhoff et al., 2007, 2011)
	LW_Atm = 0.96 * emis * SBConst * ((T_a(space))**4.0) * viewSky(space)
	! Land cover longwave radiation (Boyd and Kasper, 2003; Westhoff et al., 2007)
	LW_LC = 0.96 * (1.0 - viewSky(space)) * 0.96 * SBConst * ((T_a(space)) ** 4.0)
	! Back radiation
	! (Boyd and Kasper, 2003; Evans et al., 1998; Kim and Chapra, 1997; Ouellet et al., 2014; Westhoff et al., 2007)
	LW_Back = -0.96 * SBConst * ((StrTemp + 273.2) ** 4.0)
	
	! Atmospheric, back, and landcover longwave radiation
	LW_RAD_FLUX = LW_Atm + LW_LC + LW_Back

end function LW_RAD_FLUX
	
	
! ------------------------------------------------
!   FUNCTION Latent heat flux
! ------------------------------------------------
real function LatentHeat_FLUX(space, SW_rad, LW_rad, T_w_sim, T_a, H, V_wind, z_elv)
	
	implicit none
	REAL    :: c_air = 1004.	! Heat capacity of air (J/kg deg C)
	REAL    :: rho_air = 1.2041	! Density of air at 20 deg C (kg/m^3)
	REAL    :: StrTemp		! simulated water temperature at the previous location
	REAL    :: L_e			! Latent Heat of Vaporization
	REAL    :: satVapPreLat		! Saturation Vapor Pressure (mbar)
	REAL    :: VapPreReal		! Real vapor pressure (mbar)
	REAL    :: s			! Slope of the saturation vapor pressure curve (kPa/C)
	REAL    :: r_a			! Aerodynamic resistance (s/m)
	REAL    :: Pa_air		! Air pressure (mbar)
	REAL    :: psy			! psychrometric constant(kPa/C)
	REAL    :: E			! Penman open water evaporation
	INTEGER :: space         	! The node we are calculating the flux for that
        REAL    :: SW_rad, LW_rad
        !REAL,    INTENT(IN),  DIMENSION(:,:)  ::  T_w_sim
	REAL,    INTENT(IN),  DIMENSION(:)    ::  T_w_sim
        REAL,    INTENT(IN),  DIMENSION(:)    ::  T_a, H, V_wind, z_elv
	!INTEGER, INTENT(IN),  DIMENSION(:)    ::  BC_Grid, FROM_NODE 
	StrTemp = T_w_sim(space)
	
	! <<<<<<<<< Evaporation section >>>>>>>>
	! Calculate the latent heat of vaporization J/kg ((Maidment, 1993)
	L_e = 1000000 * (2.501 - (0.002361 * StrTemp))		
	! Calculate the slope of the saturation vapor pressure curve at a given air temperature
	! (modified from Chapra 1997 and Raudkivi 1979)
	satVapPreLat = 0.611 * exp((17.27 * T_a(space)-273.2) / (T_a(space)-273.2))
	! Vapor Pressure Based on Saturation Vapor Pressure Using Stream Temperature
	! (Bedient and Huber 1992)
	VapPreReal = (H(space) / 100.0) * satVapPreLat
	! Slope of the saturation vapor pressure curve (kPa/C) is calculated as 
	! (Maidment, 1993)
	s = (4100.0 * satVapPreLat) / ((T_a(space)-273.2) ** 2.0)
	! Calculate the aerodynamic resistance (s/m)
	! (Westhoff et al., 2007)
	r_a = 245.0 / ((0.54 * V_wind(space)) + 0.5)
	! Calculate air pressure
	Pa_air = 101.3 - (0.0105 * z_elv(space))
	! Calculate psychrometric constant(kPa / deg C) (based on air pressure(Pa), 
	! (value should be adjusted for different site elevations)
	! Ratio of water to dry air = .622, and the latent heat of water vaporization = 2.45E6 (J/kgC)
	! (Cuenca 1989)
	psy = (c_air * Pa_air) / (0.622 * 2450000.0)
	! Calculate the Penman open water evaporation
	! (Maidment, 1993; Westhoff et al., 2007)
	E = ((s * (SW_rad + LW_rad)) / (1000.0 * L_e * (s + psy))) &
		+ ((c_air * rho_air * psy * (satVapPreLat - VapPreReal)) / (1000.0 * L_e * r_a * (s + psy)))
	
	LatentHeat_FLUX = -rho_water * L_e * E
	

end function LatentHeat_FLUX


! ------------------------------------------------
!   FUNCTION Sensible heat flux
! ------------------------------------------------
real function SensibleHeat_FLUX(space, T_w_sim, LH, T_a, H, z_elv)
	
	implicit none
	REAL    :: c_air = 1004.	! Heat capacity of air (J/kgC)
	REAL    :: rho_air = 1.2041	! Density of air at 20 deg C (kg/m3)
	REAL    :: ews			! Saturation vapor pressure of the evaporating surface (KPa)
	REAL    :: ewa			! Actual vapor pressure of the evaporation surface (KPa)			
	REAL    :: StrTemp		! Simulated water temperature at the previous location
	REAL    :: Pa_air		! Air pressure (mbar)
	REAL    :: H_tmp		! Temporary value for assigning the humidity from the input
	REAL    :: B_r			! Bowen's ratio
	INTEGER :: space
        REAL    :: LH
        !REAL,    INTENT(IN),  DIMENSION(:,:)  ::  T_w_sim
	REAL,    INTENT(IN),  DIMENSION(:)    ::  T_w_sim
        REAL,    INTENT(IN),  DIMENSION(:)    ::  T_a, H, z_elv
	!INTEGER, INTENT(IN),  DIMENSION(:)    ::  BC_Grid, FROM_NODE
       
	StrTemp = T_w_sim(space)
	
	! Saturation vapor pressure of the evaporating surface (KPa)
	! (Dingman, 1994; Maidment, 1993)
	ews = 0.61275 * exp((17.27 * StrTemp) / (StrTemp))
	! Actual vapor pressure of the evaporation surface (KPa)
	! To handle the infinte problem in B_r equation, consider the umidity of 100 as 99
	H_tmp = H(space)
        if (H(space) .EQ. 100.0) then 
            H_tmp = 99.0 
        endif
	ewa = (H_tmp / 100.0) * ews
	! Adiabatic atmospheric pressure adjusted for elevation above sea level
	! (Westhoff et al., 2007)
	Pa_air = 101.3 * (((293 - (0.0065 * z_elv(space))) / 293)**5.256)
	! Bowen ratio
	! (Evans et al.,1998; Magnusson et al., 2012; Webb and Zhang, 1997; Westhoff et al., 2007)
	! (modified for units from Bedient and Huber 1992)
	B_r = 0.00061 * Pa_air * ((StrTemp - (T_a(space)-273.2)) / (ews - ewa))

	SensibleHeat_FLUX = B_r * LH
	
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
	Sed_ThermalDiffuse = 0.0000045
	H2O_ThermalDiffuse = 0.00000014331
	Sed_Density = 1600.0
	Sed_HeatCapacity = 2219.0

	! ======================================================
	! Calculate Volumetric Ratio of Water and Substrate
	! Code uses this ratio to estimate conduction constants
	! First, calculating the porosity
	! Ratio Size of dominant substrate
	Dummy1 = P_s(space) * (1.0 - C_E(space))
	! Ratio Conductivity of sand - low range
	Dummy2 = 0.062 * C_E(space)
	! Estimated Porosity
	! (adopted from Bedient and Huber 1992)
	P = (0.3683 * ((Dummy1 + Dummy2)**(-0.0641)))
	! ======================================================
	! Variables used in bed conduction
	! Calculate the sediment depth (conduction layer)
	! (Assumes 100 Particles Thickness)
	! 0.1m < particleSize < 0.3m
	Sed_Depth = 10.0 * P_s(space) / 1000.0;
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

	ThermalDiffuse = (Sed_ThermalDiffuse * Ratio_Sediment) + (H2O_ThermalDiffuse * Ratio_H2O)
	Density = (Sed_Density * Ratio_Sediment) + (rho_water * Ratio_H2O)
	HeatCapacity = (Sed_HeatCapacity * Ratio_Sediment) + (c_water * Ratio_H2O)

	BedSediment_FLUX = ThermalDiffuse * Density * HeatCapacity * ((T_bed(space)-273.) - StrTemp) / (Sed_Depth / 2.0)


end function BedSediment_FLUX

end module module_river_temperature
