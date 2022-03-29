!  Program Name:
!  Author(s)/Contact(s):
!  Abstract:
!  History Log:
! 
!  Usage:
!  Parameters: <Specify typical arguments passed>
!  Input Files:
!        <list file names and briefly describe the data they include>
!  Output Files:
!        <list file names and briefly describe the information they include>
! 
!  Condition codes:
!        <list exit condition or error codes returned >
!        If appropriate, descriptive troubleshooting instructions or
!        likely causes for failures could be mentioned here with the
!        appropriate error code
! 
!  User controllable options: <if applicable>

module module_NoahMP_hrldas_driver

  USE module_hrldas_netcdf_io
  USE module_sf_noahmp_groundwater
  USE module_sf_noahmpdrv, only: noahmp_init, noahmplsm
  USE module_date_utilities
!  USE module_sf_watertracer, only: tell_wt_flag, initialized_wt 
  use module_mpp_land, only: MPP_LAND_PAR_INI, mpp_land_init, getLocalXY, mpp_land_bcast_char, mpp_land_sync
  use module_mpp_land, only: check_land, node_info, fatal_error_stop, numprocs
  use module_cpl_land, only: cpl_land_init
  use module_NWM_io, only: output_NoahMP_NWM

  IMPLICIT NONE

  include "mpif.h"

   REAL,    allocatable, DIMENSION(:,:)   :: infxsrt,sfcheadrt, soldrain
   REAL,    allocatable, DIMENSION(:,:)   :: infxsrt_tr, sfcheadrt_tr, soldrain_tr    ! water tracer vars
   !LRK - Remove HRLDAS_ini_typ for WRF-Hydro
   integer :: forc_typ, snow_assim, HRLDAS_ini_typ
   !integer :: forc_typ, snow_assim
   REAL,    allocatable, DIMENSION(:,:)   :: etpnd, greenfrac, prcp0
   real :: etpnd1
   character(len=19) :: forcDate
   ! LRK - Remove GEO_STATIC_FLNM
   !character(len = 256):: GEO_STATIC_FLNM
   real, allocatable, dimension(:) :: zsoil    
   integer :: kk
   integer  :: finemesh, finemesh_factor
 ! INOUT (new accumulator variables for output water balance) 
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACCPRCP ! accumulated precip [mm]
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACCECAN  ! accumulated canopy evap [mm]
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACCETRAN ! accumulated transpiration [mm]
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACCEDIR ! accumulated direct soil evap [mm]
   integer :: io_config_outputs=0
   integer :: t0OutputFlag=0   
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SOILSAT_TOP ! top 2 layer soil saturation [fraction]
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SOILSAT ! column integrated soil saturation [fraction]
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SOILICE ! fraction of soil moisture that is ice [fraction]
   REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWT_AVG ! snowpack average temperature (by layer mass) [K]
   character(len=15), DIMENSION(100)        ::  IOCVARS

  character(len=9), parameter :: version = "v20150506"
  integer :: LDASIN_VERSION

!------------------------------------------------------------------------
! Begin exact copy of declaration section from driver (substitute allocatable, remove intent)
!------------------------------------------------------------------------

! IN only (as defined in WRF)

  INTEGER                                 ::  ITIMESTEP ! timestep number
  INTEGER                                 ::  YR        ! 4-digit year
  REAL                                    ::  JULIAN_IN ! Julian day
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  COSZEN    ! cosine zenith angle
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XLAT_URB2D! latitude [rad]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  DZ8W      ! thickness of atmo layers [m]
  REAL                                    ::  DTBL      ! timestep [s]
  REAL,    ALLOCATABLE, DIMENSION(:)      ::  DZS       ! thickness of soil layers [m]
  INTEGER                                 ::  NSOIL     ! number of soil layers
  INTEGER                                 ::  NUM_SOIL_LAYERS     ! number of soil layers
  REAL                                    ::  DX        ! horizontal grid spacing [m]
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    ::  IVGTYP    ! vegetation type
  INTEGER, ALLOCATABLE, DIMENSION(:,:)    ::  ISLTYP    ! soil type
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  VEGFRA    ! vegetation fraction []
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TMN       ! deep soil temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XLAND     ! =2 ocean; =1 land/seaice
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XICE      ! fraction of grid that is seaice
  REAL                                    ::  XICE_THRESHOLD! fraction of grid determining seaice
  INTEGER                                 ::  ISICE     ! land cover category for ice
  INTEGER                                 ::  ISURBAN   ! land cover category for urban
  INTEGER                                 ::  ISWATER   ! land cover category for water
  INTEGER                                 ::  IDVEG     ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1   
  INTEGER                                 ::  IOPT_CRS  ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  INTEGER                                 ::  IOPT_BTR  ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  INTEGER                                 ::  IOPT_RUN  ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  INTEGER                                 ::  IOPT_SFC  ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  INTEGER                                 ::  IOPT_FRZ  ! supercooled liquid water (1-> NY06; 2->Koren99)
  INTEGER                                 ::  IOPT_INF  ! frozen soil permeability (1-> NY06; 2->Koren99)
  INTEGER                                 ::  IOPT_RAD  ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  INTEGER                                 ::  IOPT_ALB  ! snow surface albedo (1->BATS; 2->CLASS)
  INTEGER                                 ::  IOPT_SNF  ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
  INTEGER                                 ::  IOPT_TBOT ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
  INTEGER                                 ::  IOPT_STC  ! snow/soil temperature time scheme
  INTEGER                                 ::  IOPT_GLA  ! glacier option (1->phase change; 2->simple)
  INTEGER                                 ::  IOPT_RSF  ! surface resistance option (1->Zeng; 2->simple)
  INTEGER                                 ::  IZ0TLND   ! option of Chen adjustment of Czil (not used)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  T_PHY     ! 3D atmospheric temperature valid at mid-levels [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  QV_CURR   ! 3D water vapor mixing ratio [kg/kg_dry]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  U_PHY     ! 3D U wind component [m/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  V_PHY     ! 3D V wind component [m/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SWDOWN    ! solar down at surface [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GLW       ! longwave down at surface [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  P8W       ! 3D pressure, valid at interface [Pa]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINBL, RAINBL_tmp    ! precipitation entering land model [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SR        ! frozen precip ratio entering land model [-]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINCV    ! convective precip forcing [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINNCV   ! non-convective precip forcing [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RAINSHV   ! shallow conv. precip forcing [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWNCV   ! non-covective snow forcing (subset of rainncv) [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GRAUPELNCV! non-convective graupel forcing (subset of rainncv) [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  HAILNCV   ! non-convective hail forcing (subset of rainncv) [mm]

! New spatially varying fields

  CHARACTER(LEN = 256)                    ::  spatial_filename 

! INOUT (with generic LSM equivalent) (as defined in WRF)

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TSK       ! surface radiative temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  HFX       ! sensible heat flux [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QFX       ! latent heat flux [kg s-1 m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LH        ! latent heat flux [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GRDFLX    ! ground/snow heat flux [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SMSTAV    ! soil moisture avail. [not used]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SMSTOT    ! total soil water [mm][not used]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SFCRUNOFF ! accumulated surface runoff [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  UDRUNOFF  ! accumulated sub-surface runoff [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ALBEDO    ! total grid albedo []
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWC     ! snow cover fraction []
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SMOISEQ   ! volumetric soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SMOIS     ! volumetric soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SH2O      ! volumetric liquid soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  TSLB      ! soil temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOW      ! snow water equivalent [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWH     ! physical snow depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANWAT    ! total canopy water + ice [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACSNOM    ! accumulated snow melt leaving pack
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACSNOW    ! accumulated snow on grid
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EMISS     ! surface bulk emissivity
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSFC      ! bulk surface specific humidity

! INOUT (with no Noah LSM equivalent) (as defined in WRF)

  INTEGER, ALLOCATABLE, DIMENSION(:,:)    ::  ISNOWXY   ! actual no. of snow layers
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TVXY      ! vegetation leaf temperature
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TGXY      ! bulk ground surface temperature
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANICEXY  ! canopy-intercepted ice (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANLIQXY  ! canopy-intercepted liquid water (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EAHXY     ! canopy air vapor pressure (pa)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TAHXY     ! canopy air temperature (k)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CMXY      ! bulk momentum drag coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHXY      ! bulk sensible heat exchange coefficient
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FWETXY    ! wetted or snowed fraction of the canopy (-)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNEQVOXY  ! snow mass at last time step(mm h2o)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ALBOLDXY  ! snow albedo at last time step (-)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSNOWXY   ! snowfall on the ground [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WSLAKEXY  ! lake water storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ZWTXY     ! water table depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WAXY      ! water in the "aquifer" [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WTXY      ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SMCWTDXY  ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  DEEPRECHXY! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RECHXY    ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  TSNOXY    ! snow temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  ZSNSOXY   ! snow layer depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SNICEXY   ! snow layer ice [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SNLIQXY   ! snow layer liquid water [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LFMASSXY  ! leaf mass [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RTMASSXY  ! mass of fine roots [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  STMASSXY  ! stem mass [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WOODXY    ! mass of wood (incl. woody roots) [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  STBLCPXY  ! stable carbon in deep soil [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FASTCPXY  ! short-lived carbon, shallow soil [g/m2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LAI       ! leaf area index
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XSAIXY    ! stem area index
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TAUSSXY   ! snow age factor

! OUT (with no Noah LSM equivalent) (as defined in WRF)
   
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  T2MVXY    ! 2m temperature of vegetation part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  T2MBXY    ! 2m temperature of bare ground part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  Q2MVXY    ! 2m mixing ratio of vegetation part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  Q2MBXY    ! 2m mixing ratio of bare ground part
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TRADXY    ! surface radiative temperature (k)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  NEEXY     ! net ecosys exchange (g/m2/s CO2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GPPXY     ! gross primary assimilation [g/m2/s C]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  NPPXY     ! net primary productivity [g/m2/s C]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FVEGXY    ! Noah-MP vegetation fraction [-]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RUNSFXY   ! surface runoff [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RUNSBXY   ! subsurface runoff [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ECANXY    ! evaporation of intercepted water (mm/s)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EDIRXY    ! soil surface evaporation rate (mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ETRANXY   ! transpiration rate (mm/s)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FSAXY     ! total absorbed solar radiation (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FIRAXY    ! total net longwave rad (w/m2) [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  APARXY    ! photosyn active energy by canopy (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PSNXY     ! total photosynthesis (umol co2/m2/s) [+]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SAVXY     ! solar rad absorbed by veg. (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SAGXY     ! solar rad absorbed by ground (w/m2)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RSSUNXY   ! sunlit leaf stomatal resistance (s/m)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RSSHAXY   ! shaded leaf stomatal resistance (s/m)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  ALBSNDXY  ! snow albedo (direct)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  ALBSNIXY  ! snow albedo (diffuse)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  BGAPXY    ! between gap fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WGAPXY    ! within gap fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TGVXY     ! under canopy ground temperature[K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TGBXY     ! bare ground temperature [K]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHVXY     ! sensible heat exchange coefficient vegetated
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHBXY     ! sensible heat exchange coefficient bare-ground
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SHGXY     ! veg ground sen. heat [w/m2]   [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SHCXY     ! canopy sen. heat [w/m2]   [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SHBXY     ! bare sensible heat [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EVGXY     ! veg ground evap. heat [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EVBXY     ! bare soil evaporation [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GHVXY     ! veg ground heat flux [w/m2]  [+ to soil]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GHBXY     ! bare ground heat flux [w/m2] [+ to soil]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  IRGXY     ! veg ground net LW rad. [w/m2] [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  IRCXY     ! canopy net LW rad. [w/m2] [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  IRBXY     ! bare net longwave rad. [w/m2] [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TRXY      ! transpiration [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EVCXY     ! canopy evaporation heat [w/m2]  [+ to atm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHLEAFXY  ! leaf exchange coefficient 
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHUCXY    ! under canopy exchange coefficient 
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHV2XY    ! veg 2m exchange coefficient 
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHB2XY    ! bare 2m exchange coefficient 
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  Z0        ! roughness length output to WRF
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ZNT       ! roughness length output to WRF
  INTEGER   ::  ids,ide, jds,jde, kds,kde,  &  ! d -> domain
   &            ims,ime, jms,jme, kms,kme,  &  ! m -> memory
   &            its,ite, jts,jte, kts,kte      ! t -> tile

!------------------------------------------------------------------------
! Needed for Water tracer capability
!------------------------------------------------------------------------

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QFX_TR       ! latent heat flux [kg s-1 m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  LH_TR        ! latent heat flux [W m-2]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SFCRUNOFF_TR ! accumulated surface runoff [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  UDRUNOFF_TR  ! accumulated sub-surface runoff [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SMOIS_TR     ! volumetric soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SH2OXY_TR      ! volumetric liquid soil moisture [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWXY_TR      ! snow water equivalent [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SNOWHXY_TR     ! physical snow depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANWAT_TR    ! total canopy water + ice [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACSNOM_TR    ! accumulated snow melt leaving pack
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ACSNOW_TR    ! accumulated snow on grid
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANICEXY_TR  ! canopy-intercepted ice (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CANLIQXY_TR  ! canopy-intercepted liquid water (mm)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSNOWXY_TR   ! snowfall on the ground [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WSLAKEXY_TR  ! lake water storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WAXY_TR      ! water in the "aquifer" [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  WTXY_TR      ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SMCWTDXY_TR  ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  DEEPRECHXY_TR! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RECHXY_TR    ! groundwater storage [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  ZSNSOXY_TR   ! snow layer depth [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SNICEXY_TR   ! snow layer ice [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SNLIQXY_TR   ! snow layer liquid water [mm]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RUNSFXY_TR   ! surface runoff [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RUNSBXY_TR   ! subsurface runoff [mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ECANXY_TR    ! evaporation of intercepted water (mm/s)
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EDIRXY_TR    ! soil surface evaporation rate (mm/s]
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  ETRANXY_TR   ! transpiration rate (mm/s)
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  FLUXXY_TR    ! tracer flux at soil bottom [m]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SMOIS_SUB    ! volumetric soil moisture in sublayers [m3/m3]
  REAL,    ALLOCATABLE, DIMENSION(:,:,:)  ::  SH2OXY_SUB     ! volumetric liquid soil moisture in sublayers [m3/m3]
  !------namelist variables----------------
  integer                                 :: water_tracer_option ! 0-water tracer off, 1-water tracer on
  character(len=256)                      :: track_start   ! YYYY-MM-DD
  character(len=256)                      :: track_end     ! YYYY-MM-DD
  integer                                 :: partial_mixing_option ! partial mixing option (0-complete mixing, 1-partial mixing)
  !-----variables read fro namelist--------
  integer                                 :: OPT_WT       ! 0-water tracer off, 1-water tracer on
  character(len=19)                       :: wvtstart     ! YYYY-MM-DD_hh:mm:ss
  character(len=19)                       :: wvtend       ! YYYY-MM-DD_hh:mm:ss
  integer                                 :: IOPT_PM      ! partial mixing option (0-complete mixing, 1-partial mixing)
  integer                                 :: NTRACER      ! sum of NSUB
  !------tracer flag and raio variables to "tag" water from precip------
  INTEGER                                 :: wvtflag      ! 0->before track_start; 1->at track_start
                                                          ! 2->during tracking period; 3->after track_end
  REAL                                    :: wvt_ratio    ! ratio of tracer precip to total precip
  INTEGER                                 :: K

!------------------------------------------------------------------------
! Needed for NoahMP init
!------------------------------------------------------------------------

  LOGICAL                                 ::  FNDSOILW    ! soil water present in input
  LOGICAL                                 ::  FNDSNOWH    ! snow depth present in input
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  CHSTARXY    ! for consistency with MP_init; delete later
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  SEAICE      ! seaice fraction

!------------------------------------------------------------------------
! Needed for MMF_RUNOFF (IOPT_RUN = 5); not part of MP driver in WRF
!------------------------------------------------------------------------

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  MSFTX
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  MSFTY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  EQZWT
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RIVERBEDXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  RIVERCONDXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  PEXPXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  FDEPTHXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  AREAXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QRFSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSPRINGSXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QRFXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSPRINGXY
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  QSLATXY
  REAL                                    ::  WTDDT  = 30.0    ! frequency of groundwater call [minutes]
  INTEGER                                 ::  STEPWTD          ! step of groundwater call

!------------------------------------------------------------------------
! 2D variables not used in WRF - should be removed?
!------------------------------------------------------------------------

  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  XLONIN      ! longitude
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  TERRAIN     ! terrain height
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GVFMIN      ! annual minimum in vegetation fraction
  REAL,    ALLOCATABLE, DIMENSION(:,:)    ::  GVFMAX      ! annual maximum in vegetation fraction

!------------------------------------------------------------------------
! End 2D variables not used in WRF
!------------------------------------------------------------------------

  CHARACTER(LEN=256) :: MMINSL  = 'STAS'  ! soil classification
  CHARACTER(LEN=256) :: LLANDUSE          ! (=USGS, using USGS landuse classification)

!------------------------------------------------------------------------
! Timing:
!------------------------------------------------------------------------

  INTEGER :: NTIME          ! total timesteps
  INTEGER :: ITIME          ! LSM time step
  integer :: clock_count_1 = 0
  integer :: clock_count_2 = 0
  integer :: clock_rate    = 0
  real    :: timing_sum    = 0.0

  integer :: sflx_count_sum
  integer :: count_before_sflx
  integer :: count_after_sflx

!---------------------------------------------------------------------
!  DECLARE/Initialize constants
!---------------------------------------------------------------------

    INTEGER                             :: I
    INTEGER                             :: J
    INTEGER                             :: SLOPETYP
    INTEGER                             :: YEARLEN
    INTEGER, PARAMETER                  :: NSNOW = 3    ! number of snow layers fixed to 3
    REAL, PARAMETER                     :: undefined_real = 9.9692099683868690E36 ! NetCDF float   FillValue
    INTEGER, PARAMETER                  :: undefined_int = -2147483647            ! NetCDF integer FillValue
    LOGICAL                             :: update_lai, update_veg

!---------------------------------------------------------------------
!  File naming, parallel
!---------------------------------------------------------------------

  character(len=19)  :: olddate, newdate, startdate
  character          :: hgrid, hgrid_hydro
  integer            :: igrid, igrid_hydro
  logical            :: lexist
  integer            :: imode, iimode
  integer            :: ixfull
  integer            :: jxfull
  integer            :: ixpar
  integer            :: jxpar
  integer            :: xstartpar
  integer            :: ystartpar
  integer            :: rank = 0
  CHARACTER(len=256) :: inflnm, outflnm, inflnm_template
  logical            :: restart_flag
  character(len=256) :: restart_flnm
  integer            :: ierr

!---------------------------------------------------------------------
! Attributes from LDASIN input file (or HRLDAS_SETUP_FILE, as the case may be)
!---------------------------------------------------------------------

  INTEGER           :: IX
  INTEGER           :: JX
  REAL              :: DY
  REAL              :: TRUELAT1
  REAL              :: TRUELAT2
  REAL              :: CEN_LON
  INTEGER           :: MAPPROJ
  REAL              :: LAT1
  REAL              :: LON1

  integer ix_tmp, jx_tmp

!---------------------------------------------------------------------
!  NAMELIST start
!---------------------------------------------------------------------

  character(len=256) :: indir
  ! nsoil defined above
  integer            :: forcing_timestep
  integer            :: noah_timestep
  integer            :: start_year
  integer            :: start_month
  integer            :: start_day
  integer            :: start_hour
  integer            :: start_min
  character(len=256) :: outdir = "."
  character(len=256) :: restart_filename_requested = " "
  integer            :: restart_frequency_hours
  integer            :: output_timestep

  integer            :: dynamic_veg_option
  integer            :: canopy_stomatal_resistance_option
  integer            :: btr_option
  integer            :: runoff_option
  integer            :: surface_drag_option
  integer            :: supercooled_water_option
  integer            :: frozen_soil_option
  integer            :: radiative_transfer_option
  integer            :: snow_albedo_option
  integer            :: pcp_partition_option
  integer            :: tbot_option
  integer            :: temp_time_scheme_option
  integer            :: glacier_option
  integer            :: surface_resistance_option

  integer            :: split_output_count = 1
  integer            :: khour
  integer            :: kday
  real               :: zlvl 
  character(len=256) :: hrldas_setup_file = " "
  character(len=256) :: mmf_runoff_file = " "
  character(len=256) :: external_veg_filename_template = " "
  character(len=256) :: external_lai_filename_template = " "
  integer            :: xstart = 1
  integer            :: ystart = 1
  integer            ::   xend = 0
  integer            ::   yend = 0
  integer, PARAMETER    :: MAX_SOIL_LEVELS = 20   ! maximum soil levels in namelist
  REAL, DIMENSION(MAX_SOIL_LEVELS) :: soil_thick_input       ! depth to soil interfaces from namelist [m]
  integer, dimension(MAX_SOIL_LEVELS)     :: tracer_sublayer ! number of sublayers for tracer
  integer, dimension(MAX_SOIL_LEVELS)     :: NSUB         ! number of sublayer for tracer
  integer :: rst_bi_out, rst_bi_in !0: default netcdf format. 1: binary write/read by each core.

  character(len=19)  :: tmpdate
  character(len=1000) :: VARLIST
  integer :: brkflag = 0
  integer :: varind = 1
  namelist /WRF_HYDRO_OFFLINE/ &
       !LRK - Remove HRLDAS_ini_typ and GEO_STATIC_FLNM for WRF-Hydro
       finemesh,finemesh_factor,forc_typ, snow_assim
       !finemesh,finemesh_factor,forc_typ, snow_assim , GEO_STATIC_FLNM, HRLDAS_ini_typ
  namelist / NOAHLSM_OFFLINE /    &
       indir, nsoil, soil_thick_input, forcing_timestep, noah_timestep, &
       start_year, start_month, start_day, start_hour, start_min, &
       outdir, &
       restart_filename_requested, restart_frequency_hours, output_timestep, &

       water_tracer_option, track_start, track_end, partial_mixing_option, tracer_sublayer, &

       dynamic_veg_option, canopy_stomatal_resistance_option, &
       btr_option, runoff_option, surface_drag_option, supercooled_water_option, &
       frozen_soil_option, radiative_transfer_option, snow_albedo_option, &
       pcp_partition_option, tbot_option, temp_time_scheme_option, &
       glacier_option, surface_resistance_option, &

       split_output_count, & 
       khour, kday, zlvl, hrldas_setup_file, mmf_runoff_file, &
       spatial_filename, &
       external_veg_filename_template, external_lai_filename_template, &
       xstart, xend, ystart, yend, rst_bi_out, rst_bi_in

  contains

  subroutine land_driver_ini(NTIME_out,wrfits,wrfite,wrfjts,wrfjte)
     use module_HYDRO_drv, only: HYDRO_ini
     use module_namelist, only: nlst_rt
     use module_rt_data, only: rt_domain
     use module_hrldas_HYDRO, only: open_print_mpp
     implicit  none
!---------------------------------------------
! Define interface for hrldas_drv_HYDRO_ini [Huancui Hu]
!---------------------------------------------
INTERFACE
  SUBROUTINE hrldas_drv_HYDRO_ini(STC_io,SMC_io,SH2OX_io,infxsrt,sfcheadrt,soldrain,&
                                  ii,jj,kk,kt,dt,olddate,zsoil   &
                                 ,NTRACER,OPT_WT, SMOIS_TR, SH2OXY_TR,SMOIS_SUB,SH2OXY_SUB &    !wt vars
                                 ,infxsrt_tr,sfcheadrt_tr,soldrain_tr)       !wt vars
          implicit none
          integer,            intent(in)                             :: ii,jj,kk, kt
          character(len=*),   intent(in)                             :: olddate
          real                                                       :: dt
          real,               intent(inout),dimension(ii,kk,jj)      :: STC_io,SMC_io,SH2OX_io
          real,               intent(inout),dimension(ii,jj)         :: infxsrt,sfcheadrt, soldrain
          real,               intent(in),   dimension(kk)            :: zsoil
!--------------------------------------------------------
! Needed for water tracer capability
!--------------------------------------------------------
          integer,            intent(in)                             :: NTRACER
          integer,optional,   intent(in)                             :: OPT_WT
          real,   optional,   intent(inout),dimension(ii,NTRACER,jj) :: SMOIS_TR, SH2OXY_TR, SMOIS_SUB,SH2OXY_SUB
          real,   optional,   intent(inout),dimension(ii,jj)         :: infxsrt_tr,sfcheadrt_tr, soldrain_tr
  END SUBROUTINE hrldas_drv_HYDRO_ini
END INTERFACE
!---------------------------------------------

     integer:: NTIME_out
     integer, parameter :: did=1

    ! initilization for stand alone parallel code.
    integer, optional, intent(in) :: wrfits,wrfite,wrfjts,wrfjte
    call  MPP_LAND_INIT()

! Initialize namelist variables to dummy values, so we can tell
! if they have not been set properly.

  nsoil                   = -999
  soil_thick_input        = -999
  dtbl                    = -999
  start_year              = -999
  start_month             = -999
  start_day               = -999
  start_hour              = -999
  start_min               = -999
  khour                   = -999
  kday                    = -999
  zlvl                    = -999
  forcing_timestep        = -999
  noah_timestep           = -999
  output_timestep         = -999
  restart_frequency_hours = -999
  tracer_sublayer         = -999       ! water tracer sublayers
  wvtflag                 = -999       ! flag for water tracer
  nlst_rt(did)%wvtflag    = -999
  nlst_rt(did)%NTRACER    = -999
  nlst_rt(did)%PMOPT      = -999

  open(30, file="namelist.hrldas", form="FORMATTED")
  read(30, NML=NOAHLSM_OFFLINE, iostat=ierr)
  if (ierr /= 0) then
     write(*,'(/," ***** ERROR: Problem reading namelist NOAHLSM_OFFLINE",/)')
     rewind(30)
     read(30, NOAHLSM_OFFLINE)
     stop "FATAL ERROR: Problem reading namelist NOAHLSM_OFFLINE"
  endif
  read(30, NML=WRF_HYDRO_OFFLINE, iostat=ierr)
  if (ierr /= 0) then
     write(*,'(/," ***** ERROR: Problem reading namelist WRF_HYDRO_OFFLINE",/)')
     rewind(30)
     read(30, NOAHLSM_OFFLINE)
     call hydro_stop (" FATAL ERROR: Problem reading namelist WRF_HYDRO_OFFLINE")
  endif


  close(30)

  dtbl = real(noah_timestep)
  num_soil_layers = nsoil      ! because surface driver uses the long form
  IDVEG = dynamic_veg_option ! transfer from namelist to driver format
  IOPT_CRS = canopy_stomatal_resistance_option
  IOPT_BTR = btr_option
  IOPT_RUN = runoff_option
  IOPT_SFC = surface_drag_option
  IOPT_FRZ = supercooled_water_option
  IOPT_INF = frozen_soil_option
  IOPT_RAD = radiative_transfer_option
  IOPT_ALB = snow_albedo_option
  IOPT_SNF = pcp_partition_option
  IOPT_TBOT = tbot_option
  IOPT_STC = temp_time_scheme_option
  IOPT_GLA = glacier_option
  IOPT_RSF = surface_resistance_option
  !-----water tracer variables--------
  OPT_WT   = water_tracer_option
  NTRACER  = 1
  if (OPT_WT .eq. 1) then 
    nlst_rt(did)%PMOPT   = partial_mixing_option
    if (partial_mixing_option .eq. 1 .or. partial_mixing_option .eq. 3) then
      IOPT_PM  = 1
    else
      IOPT_PM  = 0
    end if
    NSUB     = tracer_sublayer
    NTRACER  = sum(NSUB(1:nsoil))
    nlst_rt(did)%NTRACER = NTRACER
    allocate(nlst_rt(did)%NSUB(MAX_SOIL_LEVELS))
    nlst_rt(did)%NSUB    = NSUB
  end if
  print *, "FLAG1"
  print *, "OPT_WT,IOPT_PM:", OPT_WT, IOPT_PM
  print *, "NTRACER", NTRACER 

  !!----------------------------------------------------------------------------
  !! channel_only
  call updateNameList("channel_only",0)
  if(forc_typ .eq. 9)  call updateNameList("channel_only",1)
  call updateNameList("channelBucket_only",0)
  if(forc_typ .eq. 10) call updateNameList("channelBucket_only",1)

  if(forc_typ .eq. 9 .or. forc_typ .eq. 10) then
     write(olddate,'(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
          start_year, start_month, start_day, start_hour, start_min, 0
     forcdate = olddate
     if ((khour < 0) .and. (kday < 0)) then
        write(*, '("FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini() - "// &
             "Namelist error: Either KHOUR or KDAY must be defined.")')
        stop
     else if (( khour < 0 ) .and. (kday > 0)) then
        khour = kday * 24
     else if ((khour > 0) .and. (kday > 0)) then
        write(*, '("WARNING: In land_driver_ini() - Check Namelist: KHOUR and KDAY both defined.")')
        stop
     endif
     NTIME=(KHOUR)*3600./nint(dtbl)
     NTIME_out = NTIME 

     if(.not. RT_DOMAIN(did)%initialized) then  
        nlst_rt(did)%dt = real(noah_timestep)
        nlst_rt(did)%khour = khour ! Adding kHOUR to be used in the NWM output routines.
        nlst_rt(did)%olddate(1:19) = olddate(1:19)
        nlst_rt(did)%startdate(1:19) = olddate(1:19)
        nlst_rt(did)%nsoil = -999999
        call mpp_land_bcast_int1(nlst_rt(did)%nsoil)
        allocate(nlst_rt(did)%zsoil8(nlst_rt(did)%nsoil))
        nlst_rt(did)%zsoil8(1:nlst_rt(did)%nsoil) = zsoil(1:nlst_rt(did)%nsoil)       
        call HYDRO_ini(ntime,did,ix0=1,jx0=1)
        RT_DOMAIN(did)%initialized = .true.
     endif ! if .not. RT_DOMAIN(did)%initialized

     call open_print_mpp(6)
     
     return  !! no more init necessary if channel_only or channelBucket_only

  endif
  !!----------------------------------------------------------------------------
  !! channel_only


!---------------------------------------------------------------------
!  NAMELIST end
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  NAMELIST check begin
!---------------------------------------------------------------------

  update_lai = .true.   ! default: use LAI if present in forcing file
  if (dynamic_veg_option == 2 .or. dynamic_veg_option == 5 .or. dynamic_veg_option == 6) &
    update_lai = .false.

  update_veg = .false.  ! default: don't use VEGFRA if present in forcing file
  if (dynamic_veg_option == 1 .or. dynamic_veg_option == 6 .or. dynamic_veg_option == 7) &
    update_veg = .true.

  if (nsoil < 0) then
     stop "FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini()"// &
          " - NSOIL must be set in the namelist."
  endif

  if ((khour < 0) .and. (kday < 0)) then
     write(*, '("FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini() - "// &
                "Namelist error.")')
     write(*, '(" ***** ")')
     write(*, '(" *****      Either KHOUR or KDAY must be defined.")')
     write(*, '(" ***** ")')
     stop
  else if (( khour < 0 ) .and. (kday > 0)) then
     khour = kday * 24
  else if ((khour > 0) .and. (kday > 0)) then
     write(*, '("WARNING: In land_driver_ini() - Check Namelist: KHOUR and KDAY both defined.")')
  else
     ! all is well.  KHOUR defined
  endif

  if (forcing_timestep < 0) then
        write(*, *)
        write(*, '("FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini()- "// &
                   "Namelist error.")')
        write(*, '(" ***** ")')
        write(*, '(" *****       FORCING_TIMESTEP needs to be set greater than zero.")')
        write(*, '(" ***** ")')
        write(*, *)
        stop
  endif

  if (noah_timestep < 0) then
        write(*, *)
        write(*, '("FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini()"// &
                   " - Namelist error.")')
        write(*, '(" ***** ")')
        write(*, '(" *****       NOAH_TIMESTEP needs to be set greater than zero.")')
        write(*, '(" *****                     900 seconds is recommended.       ")')
        write(*, '(" ***** ")')
        write(*, *)
        stop
  endif

  !
  ! Check that OUTPUT_TIMESTEP fits into NOAH_TIMESTEP:
  !
  if (output_timestep /= 0) then
     if (mod(output_timestep, noah_timestep) > 0) then
        write(*, *)
        write(*, '("FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini() - "// &
                   "Namelist error.")')
        write(*, '(" ***** ")')
        write(*, '(" *****       OUTPUT_TIMESTEP should set to an integer multiple of NOAH_TIMESTEP.")')
        write(*, '(" *****            OUTPUT_TIMESTEP = ", I12, " seconds")') output_timestep
        write(*, '(" *****            NOAH_TIMESTEP   = ", I12, " seconds")') noah_timestep
        write(*, '(" ***** ")')
        write(*, *)
        stop
     endif
  endif

  !
  ! Check that RESTART_FREQUENCY_HOURS fits into NOAH_TIMESTEP:
  !
  if (restart_frequency_hours /= 0) then
     if (mod(restart_frequency_hours*3600, noah_timestep) > 0) then
        write(*, *)
        write(*, '("FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini() - "// &
                   "Namelist error.")')
        write(*, '(" *****       RESTART_FREQUENCY_HOURS (converted to seconds) should set to an ")')
        write(*, '(" *****       integer multiple of NOAH_TIMESTEP.")')
        write(*, '(" *****            RESTART_FREQUENCY_HOURS = ", I12, " hours:  ", I12, " seconds")') &
             restart_frequency_hours, restart_frequency_hours*3600
        write(*, '(" *****            NOAH_TIMESTEP           = ", I12, " seconds")') noah_timestep
        write(*, '(" ***** ")')
        write(*, *)
        stop
     endif
  endif

  if (dynamic_veg_option == 2 .or. dynamic_veg_option == 5 .or. dynamic_veg_option == 6) then
     if ( canopy_stomatal_resistance_option /= 1) then
        write(*, *)
        write(*, '("FATAL ERROR: In module_NoahMP_hrldas_driver.F land_driver_ini() - "// &
                   "Namelist error.")')
        write(*, '(" ***** ")')
        write(*, '(" *****       CANOPY_STOMATAL_RESISTANCE_OPTION must be 1 when DYNAMIC_VEG_OPTION == 2/5/6")')
        write(*, *)
        stop
     endif
  endif

!---------------------------------------------------------------------
!  NAMELIST check end
!---------------------------------------------------------------------

!----------------------------------------------------------------------
! Initialize gridded domain
!----------------------------------------------------------------------

    if(finemesh .ne. 0) then
         xstart = (wrfits-1)*finemesh_factor + 1
         xend = (wrfite-1)*finemesh_factor
         ystart = (wrfjts-1)*finemesh_factor + 1
         yend = (wrfjte-1)*finemesh_factor
         call CPL_LAND_INIT(xstart,xend, ystart,yend)
         ix_tmp = xend - xstart + 1
         jx_tmp = yend - ystart + 1
    else
print *, "flag new1"
       call read_dim(hrldas_setup_file,ix_tmp,jx_tmp)
print *, "flag new2"
       call MPP_LAND_PAR_INI(1,ix_tmp,jx_tmp,1)
       call getLocalXY(ix_tmp,jx_tmp,xstart,ystart,xend,yend)
    endif

print *, "flag new3"
  call read_hrldas_hdrinfo(hrldas_setup_file, ix, jx, xstart, xend, ystart, yend, &
       iswater, isurban, isice, llanduse, dx, dy, truelat1, truelat2, cen_lon, lat1, lon1, &
       igrid, mapproj)
print *, "flag new4"
  write(hgrid,'(I1)') igrid

  write(olddate,'(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2)') &
       start_year, start_month, start_day, start_hour, start_min, 0

  startdate = olddate

   ix = ix_tmp
   jx = jx_tmp

    forcdate = olddate

!-----water tracer tracking times-------
  if (OPT_WT .eq. 1) then 
    write(wvtstart,'(A4,"-",A2,"-",A2,"_",A2,":",A2,":",A2)') &
       track_start(1:4), track_start(5:6), track_start(7:8), track_start(9:10), track_start(11:12),0

    write(wvtend,'(A4,"-",A2,"-",A2,"_",A2,":",A2,":",A2)') &
       track_end(1:4), track_end(5:6), track_end(7:8), track_end(9:10), track_end(11:12),0
  end if
!---------------------------------------
  
  ids = xstart
  ide = xend
  jds = ystart
  jde = yend
  kds = 1
  kde = 2
  its = xstart
  ite = xend
  jts = ystart
  jte = yend
  kts = 1
  kte = 2
  ims = xstart
  ime = xend
  jms = ystart
  jme = yend
  kms = 1
  kme = 2

!---------------------------------------------------------------------
!  Allocate multi-dimension fields for subwindow calculation
!---------------------------------------------------------------------

  ixfull = xend-xstart+1
  jxfull = yend-ystart+1

  ixpar = ixfull
  jxpar = jxfull
  xstartpar = 1
  ystartpar = 1

  ALLOCATE ( COSZEN    (XSTART:XEND,YSTART:YEND) )    ! cosine zenith angle
  ALLOCATE ( XLAT_URB2D(XSTART:XEND,YSTART:YEND) )    ! latitude [rad]
  ALLOCATE ( DZ8W      (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! thickness of atmo layers [m]
  ALLOCATE ( DZS       (1:NSOIL)                   )  ! thickness of soil layers [m]
  ALLOCATE ( IVGTYP    (XSTART:XEND,YSTART:YEND) )    ! vegetation type
  ALLOCATE ( ISLTYP    (XSTART:XEND,YSTART:YEND) )    ! soil type
  ALLOCATE ( VEGFRA    (XSTART:XEND,YSTART:YEND) )    ! vegetation fraction []
  ALLOCATE ( TMN       (XSTART:XEND,YSTART:YEND) )    ! deep soil temperature [K]
  ALLOCATE ( XLAND     (XSTART:XEND,YSTART:YEND) )    ! =2 ocean; =1 land/seaice
  ALLOCATE ( XICE      (XSTART:XEND,YSTART:YEND) )    ! fraction of grid that is seaice
  ALLOCATE ( T_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D atmospheric temperature valid at mid-levels [K]
  ALLOCATE ( QV_CURR   (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D water vapor mixing ratio [kg/kg_dry]
  ALLOCATE ( U_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D U wind component [m/s]
  ALLOCATE ( V_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D V wind component [m/s]
  ALLOCATE ( SWDOWN    (XSTART:XEND,YSTART:YEND) )    ! solar down at surface [W m-2]
  ALLOCATE ( GLW       (XSTART:XEND,YSTART:YEND) )    ! longwave down at surface [W m-2]
  ALLOCATE ( P8W       (XSTART:XEND,KDS:KDE,YSTART:YEND) )  ! 3D pressure, valid at interface [Pa]
  ALLOCATE ( RAINBL    (XSTART:XEND,YSTART:YEND) )    ! total precipitation entering land model [mm]
  ALLOCATE ( RAINBL_tmp    (XSTART:XEND,YSTART:YEND) )    ! precipitation entering land model [mm]
  ALLOCATE ( SR        (XSTART:XEND,YSTART:YEND) )    ! frozen precip ratio entering land model [-]
  ALLOCATE ( RAINCV    (XSTART:XEND,YSTART:YEND) )    ! convective precip forcing [mm]
  ALLOCATE ( RAINNCV   (XSTART:XEND,YSTART:YEND) )    ! non-convective precip forcing [mm]
  ALLOCATE ( RAINSHV   (XSTART:XEND,YSTART:YEND) )    ! shallow conv. precip forcing [mm]
  ALLOCATE ( SNOWNCV   (XSTART:XEND,YSTART:YEND) )    ! non-covective snow forcing (subset of rainncv) [mm]
  ALLOCATE ( GRAUPELNCV(XSTART:XEND,YSTART:YEND) )    ! non-convective graupel forcing (subset of rainncv) [mm]
  ALLOCATE ( HAILNCV   (XSTART:XEND,YSTART:YEND) )    ! non-convective hail forcing (subset of rainncv) [mm]


! INOUT (with generic LSM equivalent) (as defined in WRF)

  ALLOCATE ( TSK       (XSTART:XEND,YSTART:YEND) )  ! surface radiative temperature [K]
  ALLOCATE ( HFX       (XSTART:XEND,YSTART:YEND) )  ! sensible heat flux [W m-2]
  ALLOCATE ( QFX       (XSTART:XEND,YSTART:YEND) )  ! latent heat flux [kg s-1 m-2]
  ALLOCATE ( LH        (XSTART:XEND,YSTART:YEND) )  ! latent heat flux [W m-2]
  ALLOCATE ( GRDFLX    (XSTART:XEND,YSTART:YEND) )  ! ground/snow heat flux [W m-2]
  ALLOCATE ( SMSTAV    (XSTART:XEND,YSTART:YEND) )  ! soil moisture avail. [not used]
  ALLOCATE ( SMSTOT    (XSTART:XEND,YSTART:YEND) )  ! total soil water [mm][not used]
  ALLOCATE ( SFCRUNOFF (XSTART:XEND,YSTART:YEND) )  ! accumulated surface runoff [m]
  ALLOCATE ( UDRUNOFF  (XSTART:XEND,YSTART:YEND) )  ! accumulated sub-surface runoff [m]
  ALLOCATE ( ALBEDO    (XSTART:XEND,YSTART:YEND) )  ! total grid albedo []
  ALLOCATE ( SNOWC     (XSTART:XEND,YSTART:YEND) )  ! snow cover fraction []
  ALLOCATE ( SMOISEQ   (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! eq volumetric soil moisture [m3/m3]
  ALLOCATE ( SMOIS     (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! volumetric soil moisture [m3/m3]
  ALLOCATE ( SH2O      (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! volumetric liquid soil moisture [m3/m3]
  ALLOCATE ( TSLB      (XSTART:XEND,1:NSOIL,YSTART:YEND) )     ! soil temperature [K]
  ALLOCATE ( SNOW      (XSTART:XEND,YSTART:YEND) )  ! snow water equivalent [mm]
  ALLOCATE ( SNOWH     (XSTART:XEND,YSTART:YEND) )  ! physical snow depth [m]
  ALLOCATE ( CANWAT    (XSTART:XEND,YSTART:YEND) )  ! total canopy water + ice [mm]
  ALLOCATE ( ACSNOM    (XSTART:XEND,YSTART:YEND) )  ! accumulated snow melt leaving pack
  ALLOCATE ( ACSNOW    (XSTART:XEND,YSTART:YEND) )  ! accumulated snow on grid
  ALLOCATE ( EMISS     (XSTART:XEND,YSTART:YEND) )  ! surface bulk emissivity
  ALLOCATE ( QSFC      (XSTART:XEND,YSTART:YEND) )  ! bulk surface specific humidity

! INOUT (with no Noah LSM equivalent) (as defined in WRF)

  ALLOCATE ( ISNOWXY   (XSTART:XEND,YSTART:YEND) )  ! actual no. of snow layers
  ALLOCATE ( TVXY      (XSTART:XEND,YSTART:YEND) )  ! vegetation leaf temperature
  ALLOCATE ( TGXY      (XSTART:XEND,YSTART:YEND) )  ! bulk ground surface temperature
  ALLOCATE ( CANICEXY  (XSTART:XEND,YSTART:YEND) )  ! canopy-intercepted ice (mm)
  ALLOCATE ( CANLIQXY  (XSTART:XEND,YSTART:YEND) )  ! canopy-intercepted liquid water (mm)
  ALLOCATE ( EAHXY     (XSTART:XEND,YSTART:YEND) )  ! canopy air vapor pressure (pa)
  ALLOCATE ( TAHXY     (XSTART:XEND,YSTART:YEND) )  ! canopy air temperature (k)
  ALLOCATE ( CMXY      (XSTART:XEND,YSTART:YEND) )  ! bulk momentum drag coefficient
  ALLOCATE ( CHXY      (XSTART:XEND,YSTART:YEND) )  ! bulk sensible heat exchange coefficient
  ALLOCATE ( FWETXY    (XSTART:XEND,YSTART:YEND) )  ! wetted or snowed fraction of the canopy (-)
  ALLOCATE ( SNEQVOXY  (XSTART:XEND,YSTART:YEND) )  ! snow mass at last time step(mm h2o)
  ALLOCATE ( ALBOLDXY  (XSTART:XEND,YSTART:YEND) )  ! snow albedo at last time step (-)
  ALLOCATE ( QSNOWXY   (XSTART:XEND,YSTART:YEND) )  ! snowfall on the ground [mm/s]
  ALLOCATE ( WSLAKEXY  (XSTART:XEND,YSTART:YEND) )  ! lake water storage [mm]
  ALLOCATE ( ZWTXY     (XSTART:XEND,YSTART:YEND) )  ! water table depth [m]
  ALLOCATE ( WAXY      (XSTART:XEND,YSTART:YEND) )  ! water in the "aquifer" [mm]
  ALLOCATE ( WTXY      (XSTART:XEND,YSTART:YEND) )  ! groundwater storage [mm]
  ALLOCATE ( SMCWTDXY  (XSTART:XEND,YSTART:YEND) )  ! soil moisture below the bottom of the column (m3m-3)
  ALLOCATE ( DEEPRECHXY(XSTART:XEND,YSTART:YEND) )  ! recharge to the water table when deep (m)
  ALLOCATE ( RECHXY    (XSTART:XEND,YSTART:YEND) )  ! recharge to the water table (diagnostic) (m)
  ALLOCATE ( TSNOXY    (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow temperature [K]
  ALLOCATE ( ZSNSOXY   (XSTART:XEND,-NSNOW+1:NSOIL,YSTART:YEND) )  ! snow layer depth [m]
  ALLOCATE ( SNICEXY   (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow layer ice [mm]
  ALLOCATE ( SNLIQXY   (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow layer liquid water [mm]
  ALLOCATE ( LFMASSXY  (XSTART:XEND,YSTART:YEND) )  ! leaf mass [g/m2]
  ALLOCATE ( RTMASSXY  (XSTART:XEND,YSTART:YEND) )  ! mass of fine roots [g/m2]
  ALLOCATE ( STMASSXY  (XSTART:XEND,YSTART:YEND) )  ! stem mass [g/m2]
  ALLOCATE ( WOODXY    (XSTART:XEND,YSTART:YEND) )  ! mass of wood (incl. woody roots) [g/m2]
  ALLOCATE ( STBLCPXY  (XSTART:XEND,YSTART:YEND) )  ! stable carbon in deep soil [g/m2]
  ALLOCATE ( FASTCPXY  (XSTART:XEND,YSTART:YEND) )  ! short-lived carbon, shallow soil [g/m2]
  ALLOCATE ( LAI       (XSTART:XEND,YSTART:YEND) )  ! leaf area index
  ALLOCATE ( XSAIXY    (XSTART:XEND,YSTART:YEND) )  ! stem area index
  ALLOCATE ( TAUSSXY   (XSTART:XEND,YSTART:YEND) )  ! snow age factor
  
! OUT (with no Noah LSM equivalent) (as defined in WRF)
   
  ALLOCATE ( T2MVXY    (XSTART:XEND,YSTART:YEND) )  ! 2m temperature of vegetation part
  ALLOCATE ( T2MBXY    (XSTART:XEND,YSTART:YEND) )  ! 2m temperature of bare ground part
  ALLOCATE ( Q2MVXY    (XSTART:XEND,YSTART:YEND) )  ! 2m mixing ratio of vegetation part
  ALLOCATE ( Q2MBXY    (XSTART:XEND,YSTART:YEND) )  ! 2m mixing ratio of bare ground part
  ALLOCATE ( TRADXY    (XSTART:XEND,YSTART:YEND) )  ! surface radiative temperature (k)
  ALLOCATE ( NEEXY     (XSTART:XEND,YSTART:YEND) )  ! net ecosys exchange (g/m2/s CO2)
  ALLOCATE ( GPPXY     (XSTART:XEND,YSTART:YEND) )  ! gross primary assimilation [g/m2/s C]
  ALLOCATE ( NPPXY     (XSTART:XEND,YSTART:YEND) )  ! net primary productivity [g/m2/s C]
  ALLOCATE ( FVEGXY    (XSTART:XEND,YSTART:YEND) )  ! Noah-MP vegetation fraction [-]
  ALLOCATE ( RUNSFXY   (XSTART:XEND,YSTART:YEND) )  ! surface runoff [mm/s]
  ALLOCATE ( RUNSBXY   (XSTART:XEND,YSTART:YEND) )  ! subsurface runoff [mm/s]
  ALLOCATE ( ECANXY    (XSTART:XEND,YSTART:YEND) )  ! evaporation of intercepted water (mm/s)
  ALLOCATE ( EDIRXY    (XSTART:XEND,YSTART:YEND) )  ! soil surface evaporation rate (mm/s]
  ALLOCATE ( ETRANXY   (XSTART:XEND,YSTART:YEND) )  ! transpiration rate (mm/s)
  ALLOCATE ( FSAXY     (XSTART:XEND,YSTART:YEND) )  ! total absorbed solar radiation (w/m2)
  ALLOCATE ( FIRAXY    (XSTART:XEND,YSTART:YEND) )  ! total net longwave rad (w/m2) [+ to atm]
  ALLOCATE ( APARXY    (XSTART:XEND,YSTART:YEND) )  ! photosyn active energy by canopy (w/m2)
  ALLOCATE ( PSNXY     (XSTART:XEND,YSTART:YEND) )  ! total photosynthesis (umol co2/m2/s) [+]
  ALLOCATE ( SAVXY     (XSTART:XEND,YSTART:YEND) )  ! solar rad absorbed by veg. (w/m2)
  ALLOCATE ( SAGXY     (XSTART:XEND,YSTART:YEND) )  ! solar rad absorbed by ground (w/m2)
  ALLOCATE ( RSSUNXY   (XSTART:XEND,YSTART:YEND) )  ! sunlit leaf stomatal resistance (s/m)
  ALLOCATE ( RSSHAXY   (XSTART:XEND,YSTART:YEND) )  ! shaded leaf stomatal resistance (s/m)
  ALLOCATE ( ALBSNDXY  (XSTART:XEND,1:2,YSTART:YEND) )  ! snow albedo (direct)
  ALLOCATE ( ALBSNIXY  (XSTART:XEND,1:2,YSTART:YEND) )  ! snow albedo (diffuse)
  ALLOCATE ( BGAPXY    (XSTART:XEND,YSTART:YEND) )  ! between gap fraction
  ALLOCATE ( WGAPXY    (XSTART:XEND,YSTART:YEND) )  ! within gap fraction
  ALLOCATE ( TGVXY     (XSTART:XEND,YSTART:YEND) )  ! under canopy ground temperature[K]
  ALLOCATE ( TGBXY     (XSTART:XEND,YSTART:YEND) )  ! bare ground temperature [K]
  ALLOCATE ( CHVXY     (XSTART:XEND,YSTART:YEND) )  ! sensible heat exchange coefficient vegetated
  ALLOCATE ( CHBXY     (XSTART:XEND,YSTART:YEND) )  ! sensible heat exchange coefficient bare-ground
  ALLOCATE ( SHGXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground sen. heat [w/m2]   [+ to atm]
  ALLOCATE ( SHCXY     (XSTART:XEND,YSTART:YEND) )  ! canopy sen. heat [w/m2]   [+ to atm]
  ALLOCATE ( SHBXY     (XSTART:XEND,YSTART:YEND) )  ! bare sensible heat [w/m2]  [+ to atm]
  ALLOCATE ( EVGXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground evap. heat [w/m2]  [+ to atm]
  ALLOCATE ( EVBXY     (XSTART:XEND,YSTART:YEND) )  ! bare soil evaporation [w/m2]  [+ to atm]
  ALLOCATE ( GHVXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground heat flux [w/m2]  [+ to soil]
  ALLOCATE ( GHBXY     (XSTART:XEND,YSTART:YEND) )  ! bare ground heat flux [w/m2] [+ to soil]
  ALLOCATE ( IRGXY     (XSTART:XEND,YSTART:YEND) )  ! veg ground net LW rad. [w/m2] [+ to atm]
  ALLOCATE ( IRCXY     (XSTART:XEND,YSTART:YEND) )  ! canopy net LW rad. [w/m2] [+ to atm]
  ALLOCATE ( IRBXY     (XSTART:XEND,YSTART:YEND) )  ! bare net longwave rad. [w/m2] [+ to atm]
  ALLOCATE ( TRXY      (XSTART:XEND,YSTART:YEND) )  ! transpiration [w/m2]  [+ to atm]
  ALLOCATE ( EVCXY     (XSTART:XEND,YSTART:YEND) )  ! canopy evaporation heat [w/m2]  [+ to atm]
  ALLOCATE ( CHLEAFXY  (XSTART:XEND,YSTART:YEND) )  ! leaf exchange coefficient 
  ALLOCATE ( CHUCXY    (XSTART:XEND,YSTART:YEND) )  ! under canopy exchange coefficient 
  ALLOCATE ( CHV2XY    (XSTART:XEND,YSTART:YEND) )  ! veg 2m exchange coefficient 
  ALLOCATE ( CHB2XY    (XSTART:XEND,YSTART:YEND) )  ! bare 2m exchange coefficient 
  ALLOCATE ( Z0        (XSTART:XEND,YSTART:YEND) )  ! roughness length output to WRF 
  ALLOCATE ( ZNT       (XSTART:XEND,YSTART:YEND) )  ! roughness length output to WRF 

  ALLOCATE ( XLONIN    (XSTART:XEND,YSTART:YEND) )  ! longitude
  ALLOCATE ( TERRAIN   (XSTART:XEND,YSTART:YEND) )  ! terrain height
  ALLOCATE ( GVFMIN    (XSTART:XEND,YSTART:YEND) )  ! annual minimum in vegetation fraction
  ALLOCATE ( GVFMAX    (XSTART:XEND,YSTART:YEND) )  ! annual maximum in vegetation fraction

!------------------------------------------------------------------------
! Needed for Water tracer capability
!------------------------------------------------------------------------
if (OPT_WT .eq. 1) then 
  ALLOCATE ( QFX_TR       (XSTART:XEND,YSTART:YEND) )  ! latent heat flux [kg s-1 m-2]
  ALLOCATE ( LH_TR        (XSTART:XEND,YSTART:YEND) )  ! latent heat flux [W m-2]
  ALLOCATE ( SFCRUNOFF_TR (XSTART:XEND,YSTART:YEND) )  ! accumulated surface runoff [m]
  ALLOCATE ( UDRUNOFF_TR  (XSTART:XEND,YSTART:YEND) )  ! accumulated sub-surface runoff [m]
  ALLOCATE ( SMOIS_TR     (XSTART:XEND,1:NTRACER,YSTART:YEND) )     ! volumetric soil moisture [m3/m3]
  ALLOCATE ( SH2OXY_TR    (XSTART:XEND,1:NTRACER,YSTART:YEND) )     ! volumetric liquid soil moisture [m3/m3]
  ALLOCATE ( SNOWXY_TR    (XSTART:XEND,YSTART:YEND) )  ! snow water equivalent [mm]
  ALLOCATE ( SNOWHXY_TR   (XSTART:XEND,YSTART:YEND) )  ! physical snow depth [m]
  ALLOCATE ( CANWAT_TR    (XSTART:XEND,YSTART:YEND) )  ! total canopy water + ice [mm]
  ALLOCATE ( ACSNOM_TR    (XSTART:XEND,YSTART:YEND) )  ! accumulated snow melt leaving pack
  ALLOCATE ( ACSNOW_TR    (XSTART:XEND,YSTART:YEND) )  ! accumulated snow on grid
  ALLOCATE ( CANICEXY_TR  (XSTART:XEND,YSTART:YEND) )  ! canopy-intercepted ice (mm)
  ALLOCATE ( CANLIQXY_TR  (XSTART:XEND,YSTART:YEND) )  ! canopy-intercepted liquid water (mm)
  ALLOCATE ( QSNOWXY_TR   (XSTART:XEND,YSTART:YEND) )  ! snowfall on the ground [mm/s]
  ALLOCATE ( WSLAKEXY_TR  (XSTART:XEND,YSTART:YEND) )  ! lake water storage [mm]
  ALLOCATE ( WAXY_TR      (XSTART:XEND,YSTART:YEND) )  ! water in the "aquifer" [mm]
  ALLOCATE ( WTXY_TR      (XSTART:XEND,YSTART:YEND) )  ! groundwater storage [mm]
  ALLOCATE ( SMCWTDXY_TR  (XSTART:XEND,YSTART:YEND) )  ! soil moisture below the bottom of the column (m3m-3)
  ALLOCATE ( DEEPRECHXY_TR(XSTART:XEND,YSTART:YEND) )  ! recharge to the water table when deep (m)
  ALLOCATE ( RECHXY_TR    (XSTART:XEND,YSTART:YEND) )  ! recharge to the water table (diagnostic) (m)
  ALLOCATE ( ZSNSOXY_TR   (XSTART:XEND,-NSNOW+1:NTRACER,YSTART:YEND) )  ! snow layer depth [m]
  ALLOCATE ( SNICEXY_TR   (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow layer ice [mm]
  ALLOCATE ( SNLIQXY_TR   (XSTART:XEND,-NSNOW+1:0,    YSTART:YEND) )  ! snow layer liquid water [mm]
  ALLOCATE ( RUNSFXY_TR   (XSTART:XEND,YSTART:YEND) )  ! surface runoff [mm/s]
  ALLOCATE ( RUNSBXY_TR   (XSTART:XEND,YSTART:YEND) )  ! subsurface runoff [mm/s]
  ALLOCATE ( ECANXY_TR    (XSTART:XEND,YSTART:YEND) )  ! evaporation of intercepted water (mm/s)
  ALLOCATE ( EDIRXY_TR    (XSTART:XEND,YSTART:YEND) )  ! soil surface evaporation rate (mm/s]
  ALLOCATE ( ETRANXY_TR   (XSTART:XEND,YSTART:YEND) )  ! transpiration rate (mm/s)
  ALLOCATE ( FLUXXY_TR     (XSTART:XEND,1:NTRACER,YSTART:YEND) )     ! tracer flux at soil bottom [m]
  ALLOCATE ( SMOIS_SUB     (XSTART:XEND,1:NTRACER,YSTART:YEND) )     ! volumetric soil moisture in sublayers [m3/m3]
  ALLOCATE ( SH2OXY_SUB    (XSTART:XEND,1:NTRACER,YSTART:YEND) )     ! volumetric liquid soil moisture in sublayers [m3/m3]
end if

!------------------------------------------------------------------------
! Needed for MMF_RUNOFF (IOPT_RUN = 5); not part of MP driver in WRF
!------------------------------------------------------------------------

  ALLOCATE ( MSFTX       (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( MSFTY       (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( EQZWT       (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( RIVERBEDXY  (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( RIVERCONDXY (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( PEXPXY      (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( FDEPTHXY    (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( AREAXY      (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( QRFSXY      (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( QSPRINGSXY  (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( QRFXY       (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( QSPRINGXY   (XSTART:XEND,YSTART:YEND) )  ! 
  ALLOCATE ( QSLATXY     (XSTART:XEND,YSTART:YEND) )  ! 

!------------------------------------------------------------------------

  ALLOCATE ( CHSTARXY  (XSTART:XEND,YSTART:YEND) )  ! for consistency with MP_init; delete later
  ALLOCATE ( SEAICE    (XSTART:XEND,YSTART:YEND) )  ! seaice fraction

   allocate( greenfrac   (XSTART:XEND,YSTART:YEND))
   greenfrac = 0

   ALLOCATE ( ACCPRCP  (XSTART:XEND,YSTART:YEND) )  ! accumulated precip [mm]
   ALLOCATE ( ACCECAN  (XSTART:XEND,YSTART:YEND) )  ! accumulated canopy evap [mm]
   ALLOCATE ( ACCETRAN  (XSTART:XEND,YSTART:YEND) )  ! accumulated transpiration [mm]
   ALLOCATE ( ACCEDIR  (XSTART:XEND,YSTART:YEND) )  ! accumulated direct soil evap [mm]

   ALLOCATE ( SOILSAT_TOP  (XSTART:XEND,YSTART:YEND) )  ! top 2 layer soil saturation [fraction]
   ALLOCATE ( SOILSAT  (XSTART:XEND,YSTART:YEND) )  ! column integrated soil saturation [fraction]
   ALLOCATE ( SOILICE  (XSTART:XEND,YSTART:YEND) )  ! fraction of soil moisture that is ice [fraction]
   ALLOCATE ( SNOWT_AVG  (XSTART:XEND,YSTART:YEND) )  ! snowpack average temperature (by layer mass) [K]

 ! Initialize accumulator variables to 0
   ACCPRCP = 0.0
   ACCECAN = 0.0
   ACCETRAN = 0.0
   ACCEDIR = 0.0

  COSZEN     = undefined_real
  XLAT_URB2D = undefined_real
  DZ8W       = undefined_real
  DZS        = undefined_real
  IVGTYP     = undefined_int
  ISLTYP     = undefined_int
  VEGFRA     = undefined_real
  GVFMAX     = undefined_real
  TMN        = undefined_real
  XLAND      = undefined_real
  XICE       = undefined_real
  T_PHY      = undefined_real
  QV_CURR    = undefined_real
  U_PHY      = undefined_real
  V_PHY      = undefined_real
  SWDOWN     = undefined_real
  GLW        = undefined_real
  P8W        = undefined_real
  RAINBL     = undefined_real
  RAINBL_tmp = undefined_real
  SR         = undefined_real
  RAINCV     = undefined_real
  RAINNCV    = undefined_real
  RAINSHV    = undefined_real
  SNOWNCV    = undefined_real
  GRAUPELNCV = undefined_real
  HAILNCV    = undefined_real
  TSK        = undefined_real
  QFX        = undefined_real
  SMSTAV     = undefined_real
  SMSTOT     = undefined_real
  SMOIS      = undefined_real
  SH2O       = undefined_real
  TSLB       = undefined_real
  SNOW       = undefined_real
  SNOWH      = undefined_real
  CANWAT     = undefined_real
  ACSNOM     = 0.0
  ACSNOW     = 0.0
  QSFC       = undefined_real
  SFCRUNOFF  = 0.0
  UDRUNOFF   = 0.0
  SMOISEQ    = undefined_real
  ALBEDO     = undefined_real
  ISNOWXY    = undefined_int
  TVXY       = undefined_real
  TGXY       = undefined_real
  CANICEXY   = undefined_real
  CANLIQXY   = undefined_real
  EAHXY      = undefined_real
  TAHXY      = undefined_real
  CMXY       = undefined_real
  CHXY       = undefined_real
  FWETXY     = undefined_real
  SNEQVOXY   = undefined_real
  ALBOLDXY   = undefined_real
  QSNOWXY    = undefined_real
  WSLAKEXY   = undefined_real
  ZWTXY      = undefined_real
  WAXY       = undefined_real
  WTXY       = undefined_real
  TSNOXY     = undefined_real
  SNICEXY    = undefined_real
  SNLIQXY    = undefined_real
  LFMASSXY   = undefined_real
  RTMASSXY   = undefined_real
  STMASSXY   = undefined_real
  WOODXY     = undefined_real
  STBLCPXY   = undefined_real
  FASTCPXY   = undefined_real
  LAI        = undefined_real
  XSAIXY     = undefined_real
  TAUSSXY    = undefined_real
  XLONIN     = undefined_real
  SEAICE     = undefined_real
  SMCWTDXY   = undefined_real
  DEEPRECHXY = 0.0
  RECHXY     = 0.0
  ZSNSOXY    = undefined_real
  GRDFLX     = undefined_real
  HFX        = undefined_real
  LH         = undefined_real
  EMISS      = undefined_real
  SNOWC      = undefined_real
  T2MVXY     = undefined_real
  T2MBXY     = undefined_real
  Q2MVXY     = undefined_real
  Q2MBXY     = undefined_real
  TRADXY     = undefined_real
  NEEXY      = undefined_real
  GPPXY      = undefined_real
  NPPXY      = undefined_real
  FVEGXY     = undefined_real
  RUNSFXY    = undefined_real
  RUNSBXY    = undefined_real
  ECANXY     = undefined_real
  EDIRXY     = undefined_real
  ETRANXY    = undefined_real
  FSAXY      = undefined_real
  FIRAXY     = undefined_real
  APARXY     = undefined_real
  PSNXY      = undefined_real
  SAVXY      = undefined_real
  FIRAXY     = undefined_real
  SAGXY      = undefined_real
  RSSUNXY    = undefined_real
  RSSHAXY    = undefined_real
  ALBSNDXY   = undefined_real
  ALBSNIXY   = undefined_real
  BGAPXY     = undefined_real
  WGAPXY     = undefined_real
  TGVXY      = undefined_real
  TGBXY      = undefined_real
  CHVXY      = undefined_real
  CHBXY      = undefined_real
  SHGXY      = undefined_real
  SHCXY      = undefined_real
  SHBXY      = undefined_real
  EVGXY      = undefined_real
  EVBXY      = undefined_real
  GHVXY      = undefined_real
  GHBXY      = undefined_real
  IRGXY      = undefined_real
  IRCXY      = undefined_real
  IRBXY      = undefined_real
  TRXY       = undefined_real
  EVCXY      = undefined_real
  CHLEAFXY   = undefined_real
  CHUCXY     = undefined_real
  CHV2XY     = undefined_real
  CHB2XY     = undefined_real
  TERRAIN    = undefined_real
  GVFMIN     = undefined_real
  GVFMAX     = undefined_real
  MSFTX      = undefined_real
  MSFTY      = undefined_real
  EQZWT      = undefined_real
  RIVERBEDXY = undefined_real
  RIVERCONDXY= undefined_real
  PEXPXY     = undefined_real
  FDEPTHXY   = undefined_real
  AREAXY     = undefined_real
  QRFSXY     = undefined_real
  QSPRINGSXY = undefined_real
  QRFXY      = undefined_real
  QSPRINGXY  = undefined_real
  QSLATXY    = undefined_real
  CHSTARXY   = undefined_real
  Z0         = undefined_real
  ZNT        = undefined_real
 
  SOILSAT_TOP = undefined_real 
  SOILSAT     = undefined_real
  SOILICE     = undefined_real
  SNOWT_AVG   = undefined_real

  XLAND          = 1.0   ! water = 2.0, land = 1.0
  XICE           = 0.0   ! fraction of grid that is seaice
  XICE_THRESHOLD = 0.5   ! fraction of grid determining seaice (from WRF)

!------initialized water tracer variables---------
if (OPT_WT .eq. 1) then 
  QFX_TR        = undefined_real
  LH_TR         = undefined_real
  SFCRUNOFF_TR  = 0.0
  UDRUNOFF_TR   = 0.0
  SMOIS_TR      = undefined_real
  SH2OXY_TR       = undefined_real
  SNOWXY_TR       = undefined_real
  SNOWHXY_TR      = undefined_real
  CANWAT_TR     = undefined_real
  ACSNOM_TR     = 0.0
  ACSNOW_TR     = 0.0
  CANICEXY_TR   = undefined_real
  CANLIQXY_TR   = undefined_real
  QSNOWXY_TR    = undefined_real
  WSLAKEXY_TR   = undefined_real
  WAXY_TR       = undefined_real
  WTXY_TR       = undefined_real
  SMCWTDXY_TR   = undefined_real
  DEEPRECHXY_TR = 0.0
  RECHXY_TR     = 0.0
  ZSNSOXY_TR    = undefined_real
  SNICEXY_TR    = undefined_real
  SNLIQXY_TR    = undefined_real
  RUNSFXY_TR    = undefined_real
  RUNSBXY_TR    = undefined_real
  ECANXY_TR     = undefined_real
  EDIRXY_TR     = undefined_real
  ETRANXY_TR    = undefined_real
  FLUXXY_TR     = undefined_real
  SMOIS_SUB      = undefined_real
  SH2OXY_SUB     = undefined_real
end if
print *, "FLAG2"
print *, "NTRACER", NTRACER
!-------------------------------------------------

!----------------------------------------------------------------------
! Read Landuse Type and Soil Texture and Other Information
!----------------------------------------------------------------------
 
  CALL READLAND_HRLDAS(HRLDAS_SETUP_FILE, XSTART, XEND, YSTART, YEND,     &
       ISWATER, IVGTYP, ISLTYP, TERRAIN, TMN, XLAT_URB2D, XLONIN, XLAND, SEAICE,MSFTX,MSFTY)
  
  WHERE(SEAICE > 0.0) XICE = 1.0

!------------------------------------------------------------------------
! For spatially-varying soil parameters, read in necessary extra fields
!------------------------------------------------------------------------

  
!------------------------------------------------------------------------
! For IOPT_RUN = 5 (MMF groundwater), read in necessary extra fields
! This option is not tested for parallel use in the offline driver
!------------------------------------------------------------------------

  if (runoff_option == 5) then
    CALL READ_MMF_RUNOFF(MMF_RUNOFF_FILE, XSTART, XEND, YSTART, YEND,&
                         ZWTXY,EQZWT,RIVERBEDXY,RIVERCONDXY,PEXPXY,FDEPTHXY)
  end if

!----------------------------------------------------------------------
! Initialize Model State
!----------------------------------------------------------------------

  SLOPETYP = 2
  DZS       =  SOIL_THICK_INPUT(1:NSOIL)

  ITIMESTEP = 1

  if (restart_filename_requested /= " ") then
     restart_flag = .TRUE.

  if(rst_bi_in .eq. 0) then

      tmpdate = olddate

     call find_restart_file(rank, trim(restart_filename_requested), startdate, khour, olddate, restart_flnm)

     call read_restart(trim(restart_flnm), xstart, xend, xstart, ixfull, jxfull, nsoil, olddate)

      olddate = tmpdate 

       ITIMESTEP = 2

     call mpp_land_bcast_char(19,OLDDATE(1:19))

     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SOIL_T"  , TSLB     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNOW_T"  , TSNOXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMC"     , SMOIS    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SH2O"    , SH2O     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ZSNSO"   , ZSNSOXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNICE"   , SNICEXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNLIQ"   , SNLIQXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSNOW"   , QSNOWXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FWET"    , FWETXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNEQVO"  , SNEQVOXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "EAH"     , EAHXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TAH"     , TAHXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ALBOLD"  , ALBOLDXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CM"      , CMXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CH"      , CHXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ISNOW"   , ISNOWXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANLIQ"  , CANLIQXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANICE"  , CANICEXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNEQV"   , SNOW     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNOWH"   , SNOWH    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TV"      , TVXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TG"      , TGXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ZWT"     , ZWTXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WA"      , WAXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WT"      , WTXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WSLAKE"  , WSLAKEXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "LFMASS"  , LFMASSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RTMASS"  , RTMASSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "STMASS"  , STMASSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WOOD"    , WOODXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "STBLCP"  , STBLCPXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FASTCP"  , FASTCPXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "LAI"     , LAI      )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SAI"     , XSAIXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "VEGFRA"  , VEGFRA   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GVFMIN"  , GVFMIN   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "GVFMAX"  , GVFMAX   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACMELT"  , ACSNOM   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACSNOW"  , ACSNOW   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "TAUSS"   , TAUSSXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSFC"    , QSFC     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SFCRUNOFF",SFCRUNOFF   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "UDRUNOFF" ,UDRUNOFF    )
!------------------------------------------------------------------------------------------
! Needed for water tracer capabilit
!------------------------------------------------------------------------------------------
  if (OPT_WT .eq. 1) then
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMC_TR"     , SMOIS_TR    )  !1
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SH2O_TR"    , SH2OXY_TR   )  !2
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNICE_TR"   , SNICEXY_TR  )  !3
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNLIQ_TR"   , SNLIQXY_TR  )  !4
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSNOW_TR"   , QSNOWXY_TR  )  !5
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANLIQ_TR"  , CANLIQXY_TR )  !6
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "CANICE_TR"  , CANICEXY_TR )  !7
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNEQV_TR"   , SNOWXY_TR   )  !8
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SNOWH_TR"   , SNOWHXY_TR  )  !9
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WA_TR"      , WAXY_TR     )  !10
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "WT_TR"      , WTXY_TR     )  !11
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACMELT_TR"  , ACSNOM_TR   )  !12
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACSNOW_TR"  , ACSNOW_TR   )  !13
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SFCRUNOFF_TR",SFCRUNOFF_TR   ) !14
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "UDRUNOFF_TR" ,UDRUNOFF_TR    ) !15
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMC_SUB"     , SMOIS_SUB    )  !16
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SH2O_SUB"    , SH2OXY_SUB   )  !17
  end if
!------------------------------------------------------------------------------------------

     if(checkRstV("ACCPRCP") .eq. 0) call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACCPRCP"  ,ACCPRCP    )
     if(checkRstV("ACCECAN") .eq. 0) call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACCECAN"  ,ACCECAN    )
     if(checkRstV("ACCEDIR") .eq. 0) call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACCEDIR"  ,ACCEDIR    )
     if(checkRstV("ACCETRAN") .eq. 0) call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "ACCETRAN" ,ACCETRAN    )
! below for opt_run = 5
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMOISEQ"   , SMOISEQ    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "AREAXY"    , AREAXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "SMCWTDXY"  , SMCWTDXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QRFXY"     , QRFXY      )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "DEEPRECHXY", DEEPRECHXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSPRINGXY" , QSPRINGXY  )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSLATXY"   , QSLATXY    )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QRFSXY"    , QRFSXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "QSPRINGSXY", QSPRINGSXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RECHXY"    , RECHXY     )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "FDEPTHXY"   ,FDEPTHXY   )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RIVERCONDXY",RIVERCONDXY)
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "RIVERBEDXY" ,RIVERBEDXY )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "EQZWT"      ,EQZWT      )
     call get_from_restart(xstart, xend, xstart, ixfull, jxfull, "PEXPXY"     ,PEXPXY     )
  else  !  rst_bi_in
     call mpp_land_bcast_char(19,OLDDATE(1:19))
     ITIMESTEP = 2
     call lsm_rst_bi_in()
  endif

     STEPWTD = nint(WTDDT*60./DTBL)
     STEPWTD = max(STEPWTD,1)

! Must still call NOAHMP_INIT even in restart to set up parameter arrays (also done in WRF)

     CALL NOAHMP_INIT(    LLANDUSE,     SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, &   ! call from WRF phys_init
                    TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH, &
                     TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                  FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, WSLAKEXY,    ZWTXY,     WAXY, &
                    WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                  T2MVXY,   T2MBXY, CHSTARXY,                                         &
                   NSOIL,  .true.,                                                   &
                  .true.,runoff_option,                                                   &
                  ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                  ims,ime, jms,jme, kms,kme,                &  ! memory
                  its,ite, jts,jte, kts,kte                 &  ! tile
                     ,smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy ,dx, dy, msftx, msfty,&
                     wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                     fdepthxy ,terrain ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy              &
                     )
  else  ! for none restart

     restart_flag = .FALSE.

     SMOIS     =  undefined_real
     TSLB      =  undefined_real
     SH2O      =  undefined_real
     CANLIQXY  =  undefined_real
     TSK       =  undefined_real
     RAINBL_tmp    =  undefined_real
     SNOW      =  undefined_real
     SNOWH     =  undefined_real

! LRK - Remove HRLDAS_ini_typ for WRF-Hydro. Originally, there was a fork 
!       for difference forcing types, or if the user desired a different
!       reading from wrfinput. However, for WRF-Hydro uses, we only use
!       option 1 at this point, so removing all other potential options,
!       which was 0 before. This will force the user to place SHDMIN/SHDMAX
!       into the wrfinput file now as opposed to allowing it to be read 
!       in from the geogrid file.
!#ifdef 1
!  if((forc_typ .gt. 2) .and. (forc_typ .ne. 6) ) HRLDAS_ini_typ = 0

!  if(HRLDAS_ini_typ .eq. 1) then
!     ! read initial parameters and conditions from the HRLDAS forcing data
!     if(forc_typ .eq. 2) then
!          inflnm = trim(indir)//"/"//&
!          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
!          startdate(15:16)//".LDASIN_DOMAIN"//hgrid
!     else
!          inflnm = trim(indir)//"/"//&
!          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
!          ".LDASIN_DOMAIN"//hgrid
!     endif

!#else

!     inflnm = trim(indir)//"/"//&
!          startdate(1:4)//startdate(6:7)//startdate(9:10)//startdate(12:13)//&
!          ".LDASIN_DOMAIN"//hgrid
!#endif

     CALL READINIT_HRLDAS(HRLDAS_SETUP_FILE, xstart, xend, ystart, yend,  &
          NSOIL, DZS, OLDDATE, LDASIN_VERSION, SMOIS,       &
          TSLB, CANWAT, TSK, SNOW, SNOWH, FNDSNOWH)
	  
!yw     VEGFRA    =  undefined_real
!yw     LAI       =  undefined_real
!yw     GVFMIN    =  undefined_real
!yw     GVFMAX    =  undefined_real

     CALL READVEG_HRLDAS(HRLDAS_SETUP_FILE, xstart, xend, ystart, yend,  &
          OLDDATE, IVGTYP, VEGFRA, LAI, GVFMIN, GVFMAX)

!#ifdef 1
!  else   !  HRLDAS_ini_typ        

!#ifdef 1
!    call HYDRO_HRLDAS_ini_mpp   &
!#else
!    call HYDRO_HRLDAS_ini   &
!#endif
!       (trim(hrldas_setup_file), xend-xstart+1,yend-ystart+1, &
!       nsoil,SMOIS(:,1:NSOIL,:),TSLB(:,1:NSOIL,:),SH2O(:,1:NSOIL,:), CANWAT, TSK,SNOW,SNOWH,lai,VEGFRA,IVGTYP,FNDSNOWH)
!       if(maxval(VEGFRA) .le. 1)  VEGFRA = VEGFRA * 100

!       greenfrac = 0.0
!#ifdef 1
!        call get_greenfrac_mpp &
!#else
!      call get_greenfrac &
!#endif
!         (trim(GEO_STATIC_FLNM),greenfrac, ix, jx, olddate, GVFMAX)
!        !yw GVFMAX = maxval(greenfrac)
!	if(maxval(GVFMAX) .le. 1)  GVFMAX = GVFMAX * 100
!  endif   ! initialization type
!#endif

!     SNOW = SNOW * 1000.    ! Convert snow water equivalent to mm. MB: remove v3.7

     FNDSOILW = .FALSE.
     CALL NOAHMP_INIT(    LLANDUSE,     SNOW,    SNOWH,   CANWAT,   ISLTYP,   IVGTYP, &   ! call from WRF phys_init
                    TSLB,    SMOIS,     SH2O,      DZS, FNDSOILW, FNDSNOWH, &
                     TSK,  ISNOWXY,     TVXY,     TGXY, CANICEXY,      TMN,     XICE, &
                CANLIQXY,    EAHXY,    TAHXY,     CMXY,     CHXY,                     &
                  FWETXY, SNEQVOXY, ALBOLDXY,  QSNOWXY, WSLAKEXY,    ZWTXY,     WAXY, &
                    WTXY,   TSNOXY,  ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY, RTMASSXY, &
                STMASSXY,   WOODXY, STBLCPXY, FASTCPXY,   XSAIXY, LAI,                    &
                  T2MVXY,   T2MBXY, CHSTARXY,                                         &
                   NSOIL,  .false.,                                                   &
                  .true.,runoff_option,                                                   &
                  ids,ide+1, jds,jde+1, kds,kde,                &  ! domain
                  ims,ime, jms,jme, kms,kme,                &  ! memory
                  its,ite, jts,jte, kts,kte                 &  ! tile
                     ,smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy ,dx, dy, msftx, msfty,&
                     wtddt    ,stepwtd  ,dtbl  ,qrfsxy ,qspringsxy  ,qslatxy,                  &
                     fdepthxy ,terrain ,riverbedxy ,eqzwt ,rivercondxy ,pexpxy              &
                     )

      TAUSSXY = 0.0   ! Need to be added to _INIT later


  endif  ! end of restart if block
  
  NTIME=(KHOUR)*3600./nint(dtbl)

  print*, "NTIME = ", NTIME , "KHOUR=",KHOUR,"dtbl = ", dtbl

  ! assinging the KHOUR to be used in the NWM output routing for global metadata
  nlst_rt(did)%khour = khour

  call system_clock(count=clock_count_1)   ! Start a timer


   allocate( infxsrt   (ix,jx) )
   allocate( sfcheadrt (ix,jx) )
   allocate( soldrain  (ix,jx) )
   allocate( etpnd     (ix,jx) )
   allocate( prcp0     (ix,jx) )
   
   !-----water tracer vars-----
   if (OPT_WT .eq. 1) then
     allocate( infxsrt_tr   (ix,jx) )
     allocate( sfcheadrt_tr (ix,jx) )
     allocate( soldrain_tr  (ix,jx) )
     infxsrt_tr  = 0.0
     sfcheadrt_tr= 0.0
     soldrain_tr = 0.0
   end if
   !-----water tracer end------   

   prcp0     = 0
   sfcheadrt = 0.0
   infxsrt   = 0.0
   etpnd     = 0.0
   soldrain  = 0.0

   allocate(zsoil (NSOIL))
   zsoil = 0

   zsoil(1) = -1*soil_thick_input(1)
   do kk = 2, NSOIL
      zsoil(kk) = zsoil(kk-1)-soil_thick_input(kk)
   end do
   print*, "zsoil/soil_thick_input = ", soil_thick_input(1:NSOIL) 

  !-------water tracer calc------
  if (OPT_WT .eq. 1) then
   call hrldas_drv_HYDRO_ini(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:), &
         infxsrt,sfcheadrt,soldrain,ix,jx,NSOIL,NSOIL,real(noah_timestep),          &
         olddate,zsoil(1:NSOIL)      &
         , NTRACER, OPT_WT=OPT_WT                           &
         , SMOIS_TR=SMOIS_TR, SH2OXY_TR=SH2OXY_TR           &
         , SMOIS_SUB=SMOIS_SUB, SH2OXY_SUB=SH2OXY_SUB       &
         , infxsrt_tr=infxsrt_tr                            &
         , sfcheadrt_tr=sfcheadrt_tr                        &
         , soldrain_tr =soldrain_tr)
  else
   call hrldas_drv_HYDRO_ini(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:), &
         infxsrt,sfcheadrt,soldrain,ix,jx,NSOIL,NSOIL,real(noah_timestep),          &
         olddate,zsoil(1:NSOIL),NTRACER)
  end if

call get_iocflag(1, io_config_outputs)

if (io_config_outputs .gt. 0) then
     ACCPRCP = 0.0
     ACCECAN = 0.0
     ACCEDIR = 0.0
     ACCETRAN = 0.0
     SFCRUNOFF = 0.0
     UDRUNOFF = 0.0
     ACSNOM = 0.0
     ACSNOW = 0.0
endif

   !!--- Setup variable list. Change as needed; rest of parsing should adapt automatically.
   if (io_config_outputs .eq. 0) then
        VARLIST = 'IVGTYP,ISLTYP,FVEG,LAI,SAI,SWFORC,COSZ,LWFORC,RAINRATE,EMISS,FSA,FIRA,GRDFLX,HFX,LH,ECAN,EDIR,ALBEDO,' // &
                  'ETRAN,UGDRNOFF,SFCRNOFF,CANLIQ,CANICE,ZWT,WA,WT,ACCPRCP,ACCECAN,ACCEDIR,ACCETRAN,SAV,TR,EVC,IRC,SHC,' // &
                  'IRG,SHG,EVG,GHV,SAG,IRB,SHB,EVB,GHB,TRAD,TG,TV,TAH,TGV,TGB,T2MV,T2MB,Q2MV,Q2MB,EAH,FWET,ZSNSO_SN,SNICE,' // &
                  'SNLIQ,SOIL_T,SOIL_W,SNOW_T,SOIL_M,SNOWH,SNEQV,QSNOW,ISNOW,FSNO,ACSNOW,ACSNOM,CM,CH,CHV,CHB,CHLEAF,CHUC,' // &
                  'CHV2,CHB2,LFMASS,RTMASS,STMASS,WOOD,STBLCP,FASTCP,NEE,GPP,NPP,PSN,APAR,ACCET,CANWAT,SOILICE,SOILSAT_TOP,'// &
                  'SOILSAT,SNOWT_AVG'
   endif
   if (io_config_outputs .eq. 1) then 
        VARLIST = 'SNOWH,SNEQV,FSNO,SOILSAT_TOP,SNOWT_AVG,ACCET'
   endif
   if (io_config_outputs .eq. 2) then 
        VARLIST = 'SNOWH,SNEQV,FSNO,SOILSAT_TOP,SNOWT_AVG,ACCET'
   endif
   if (io_config_outputs .eq. 3) then
        VARLIST = 'UGDRNOFF,ACSNOM,SNOWH,SNEQV,ACCECAN,ACCETRAN,ACCEDIR,SNLIQ,ISNOW,SOIL_T,FSNO,SOIL_M,GRDFLX,HFX,LH,FIRA,FSA,TRAD,SOILSAT_TOP,SNOWT_AVG,SOILICE,ACCET,CANWAT' 
   endif
   if (io_config_outputs .eq. 4) then 
        VARLIST = 'UGDRNOFF,SFCRNOFF,ACSNOM,SNEQV,SOILSAT_TOP,SOILSAT,ACCET,CANWAT,PET'
   endif
   if (io_config_outputs .eq. 5) then
        VARLIST = 'UGDRNOFF,SFCRNOFF,ACCET,SNEQV,SNOWH,FSNO,SOIL_M,SOIL_W,TRAD,FIRA,FSA,LH,HFX'
   endif
   if (io_config_outputs .eq. 6) then
        VARLIST = 'UGDRNOFF,SFCRNOFF,ACSNOM,SNOWH,SNEQV,ACCECAN,ACCETRAN,ACCEDIR,SNLIQ,ISNOW,SOIL_T,FSNO,SOIL_M,GRDFLX,HFX,LH,FIRA,FSA,TRAD,SOILSAT_TOP,SOILSAT,SNOWT_AVG,SOILICE,ACCET,CANWAT'
   endif

   !!--- Parse into character array. Constructor not valid with uneven 
   !!--- strings in f90 so using brute force. 
   do while (brkflag .eq. 0)
        if (index(VARLIST, ',') .eq. 0) then 
                IOCVARS(varind) = adjustl(VARLIST)
                brkflag = 1
        else 
                IOCVARS(varind) = adjustl(VARLIST(1:(index(VARLIST, ',')-1)))
                VARLIST = VARLIST((index(VARLIST, ',')+1):)
                varind = varind + 1
        endif
   end do

   if(finemesh .ne. 0 ) then
       if(restart_flag) then
          NTIME_out =  10
       else
          NTIME_out =  1 
       endif
       return
   endif

   NTIME_out = NTIME 

   call get_t0OutputFlag(1, t0OutputFlag)
   if(my_id .eq. io_id) &
       print*, "t0OutputFlag: ", t0OutputFlag  

! ldas_output subroutine will be called when 
! the t0OutputFlag is 1 in the hydro.namelist
! the ldas_output requires one variables ITIME 
! which is the LSM timestep, we declare it here
! since it does not exist at this point
  ITIME = 0
  if(t0OutputFlag .eq. 1) call ldas_output(ITIME, wvtflag)
!#else
!  if (restart_filename_requested == " ") then
!     if(t0OutputFlag .eq. 1) call ldas_output()
!  endif
!#endif
print *, "FLAG3"
print *, "NTRACER", NTRACER
end subroutine land_driver_ini

!===============================================================================
  subroutine land_driver_exe(itime)
  use module_hydro_io, only: read_channel_only
     implicit  none
     INTERFACE
        subroutine hrldas_drv_HYDRO(STC_io,SMC_io,SH2OX_io,infxsrt,sfcheadrt,soldrain,ii,jj,kk, &
                                   NTRACER, WVTFLAG, SMOIS_SUB, SH2OXY_SUB, &      !wt vars
                                   SMOIS_TR, SH2OXY_TR,infxsrt_tr, sfcheadrt_tr,soldrain_tr )   !wt vars
          implicit none
          integer,                                intent(INOUT) :: ii,jj,kk
          real,dimension(ii,kk,jj),               intent(INOUT) :: STC_io,SMC_io,SH2OX_io
          real,dimension(ii,jj),                  intent(INOUT) :: infxsrt,sfcheadrt, soldrain
         !-----water tracer vars------
          integer,                                intent(INOUT) :: NTRACER
          integer,                      optional, intent(INOUT) :: WVTFLAG
          real,dimension(ii,NTRACER,jj),optional, intent(INOUT) :: SMOIS_SUB, SH2OXY_SUB, SMOIS_TR,SH2OXY_TR
          real,dimension(ii,jj),        optional, intent(INOUT) :: infxsrt_tr,sfcheadrt_tr, soldrain_tr
        end subroutine hrldas_drv_HYDRO
     END INTERFACE


     integer :: itime          ! timestep loop

!---------------------------------------------------------------------------------
! Read the forcing data.
!---------------------------------------------------------------------------------

     call mpp_land_bcast_char(19,OLDDATE(1:19))

!      if(forc_typ .eq. 8) then
!          call read_forc_ldasout(olddate,hgrid, indir, dtbl,ix,jx,infxsrt,soldrain)
!          call hrldas_drv_HYDRO(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:),infxsrt,sfcheadrt,soldrain,ix,jx,NSOIL)
!          return
!      endif

      if(forc_typ .eq. 8) then
          call read_forc_ldasout(olddate,hgrid, indir, dtbl,ix,jx,infxsrt,soldrain)
          call hrldas_drv_HYDRO(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:),infxsrt,sfcheadrt,soldrain,ix,jx,NSOIL,NTRACER)
          call geth_newdate(newdate, olddate, nint(dtbl))
          olddate = newdate
          return
      endif
      if(forc_typ .eq. 9 .or. forc_typ .eq. 10) then
         !! JLM:: fix hgrid: This becomes 1 eventhough 3 is specified in hydro.namelist
         !! JLM: This is initalized by read_hrldas_hdrinfo
         !! JLM: Appears that we should differentiate the LSM and HYDRO igrids, define a local 
         !! JLM: igrid for this purpose.
         !! JLM: ?*?* hrldas_drv_HYDRO.F should be made a module *?*?
         !! JLM:  Simple modification which forces type, rank and kind checking...
          call getNameList('igrid', igrid_hydro)  !! get hydro namelist info :: case sensitive
          write(hgrid_hydro,'(I1)') igrid_hydro
          call read_channel_only(olddate, hgrid_hydro, indir, forcing_timestep)
          call hrldas_drv_HYDRO(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:),infxsrt,sfcheadrt,soldrain,ix,jx,NSOIL,NTRACER)
          call geth_newdate(newdate, olddate, nint(dtbl))
          olddate = newdate
          return
      endif

! For HRLDAS, we're assuming (for now) that each time period is in a 
! separate file.  So we can open a new one right now.

     inflnm = trim(indir)//"/"//&
          olddate(1:4)//olddate(6:7)//olddate(9:10)//olddate(12:13)//&
          ".LDASIN_DOMAIN"//hgrid

     ! Build a filename template
     inflnm_template = trim(indir)//"/<date>.LDASIN_DOMAIN"//hgrid



     if(finemesh .ne. 0) goto 991
    
     if(forc_typ .eq. 0) then
        CALL READFORC_HRLDAS(INFLNM_TEMPLATE, FORCING_TIMESTEP, OLDDATE,  &
             XSTART, XEND, YSTART, YEND,                                  &
             T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),          &
	       P8W(:,1,:), GLW       ,SWDOWN      ,RAINBL_tmp, VEGFRA, update_veg, LAI, update_lai)
     else
        if(olddate == forcDate) then
           CALL HYDRO_frocing_drv(trim(indir), forc_typ,snow_assim,olddate,                      &
               xstart, xend, ystart, yend,    &
               T_PHY(:,1,:),QV_CURR(:,1,:),U_PHY(:,1,:),V_PHY(:,1,:),P8W(:,1,:),    &
               GLW,SWDOWN,RAINBL_tmp,LAI,VEGFRA,SNOWH,ITIME,FORCING_TIMESTEP,prcp0)

               if(maxval(VEGFRA) .le. 1)  VEGFRA = VEGFRA * 100

           call geth_newdate(newdate, forcDate, FORCING_TIMESTEP)
           forcDate = newdate
        endif
     endif

!------water tracer variables initiation-------
if (OPT_WT .eq. 1) then 
  if (olddate < wvtstart) then
    wvtflag = 0
  else
    if (olddate(1:15) == wvtstart(1:15)) then
      wvtflag = 1
    else
      if (olddate < wvtend) then
        wvtflag = 2
      else
        wvtflag = 3
      end if
    end if
  end if
!print *, "FLAG4, SFCRUNOFF(32,26),SFCRUNOFF_TR(32,26):",SFCRUNOFF(32,26),SFCRUNOFF_TR(32,26)
print *, "wvtflag", wvtflag
print *, "NTRACER", NTRACER
  ! call tell_wt_flag(olddate, wvtstart, wvtend, wvtflag)
  if (wvtflag == 1) then
        QFX_TR        = 0.0
        LH_TR         = 0.0
        SFCRUNOFF_TR  = 0.0
        UDRUNOFF_TR   = 0.0
        SMOIS_TR      = 0.0
        SH2OXY_TR       = 0.0
        SNOWXY_TR       = 0.0
        SNOWHXY_TR      = 0.0
        CANWAT_TR     = 0.0
        ACSNOM_TR     = 0.0
        ACSNOW_TR     = 0.0
        CANICEXY_TR   = 0.0
        CANLIQXY_TR   = 0.0
        QSNOWXY_TR    = 0.0
        WSLAKEXY_TR   = 0.0
        WAXY_TR       = 0.0
        WTXY_TR       = 0.0
        ZSNSOXY_TR    = 0.0
	K = 0
        do I = 1, NSOIL
           do J = 1,NSUB(I)
              K = K + 1
              ZSNSOXY_TR(:,K,:) = ZSNSOXY(:,I,:) + DZS(I)*(NSUB(I)-J)/NSUB(I)
           end do
        end do
	smcwtdxy_TR   = 0.0
        rechxy_TR     = 0.0
        deeprechxy_TR = 0.0
        SNICEXY_TR    = 0.0
        SNLIQXY_TR    = 0.0
        RUNSFXY_TR    = 0.0
        RUNSBXY_TR    = 0.0
        ECANXY_TR     = 0.0
        EDIRXY_TR     = 0.0
        ETRANXY_TR    = 0.0
        SMOIS_SUB      = 0.0
        SH2OXY_SUB       = 0.0
        FLUXXY_TR      = 0.0
        ! specicial care for SMOIS_SUB...
        K = 0
        DO I = 1, NSOIL
              DO J =  NSUB(I), 1, -1
                 K = K + 1
                 SMOIS_SUB(:,K,:) = SMOIS(:,I,:)
                 SH2OXY_SUB(:,K,:) = SH2O(:,I,:)
              END DO
         !     if (I .eq. 1) then
         !        SMOIS_SUB(:,2:NSUB(1),:) = SH2O_SUB(:,2:NSUB(1),:)
         !        SMOIS_SUB(:,1,:) = SMOIS(:,I,:)*NSUB(1)-SH2O(:,I,:)*(NSUB(1)-1)
         !     end if
        END DO
	! call initialized_wt(xxx)
	wvtflag = 2
  end if
  
end if
!print *, "FLAG5, SFCRUNOFF(44,23),SFCRUNOFF_TR(44,23):",SFCRUNOFF(44,23),SFCRUNOFF_TR(44,23)
print *, "NTRACER", NTRACER
!----------------------------------------------

991  continue

     where(XLAND > 1.5)   T_PHY(:,1,:) = 0.0  ! Prevent some overflow problems with ifort compiler [MB:20150812]
     where(XLAND > 1.5)   U_PHY(:,1,:) = 0.0
     where(XLAND > 1.5)   V_PHY(:,1,:) = 0.0
     where(XLAND > 1.5) QV_CURR(:,1,:) = 0.0
     where(XLAND > 1.5)     P8W(:,1,:) = 0.0
     where(XLAND > 1.5)     GLW        = 0.0
     where(XLAND > 1.5)  SWDOWN        = 0.0
     where(XLAND > 1.5) RAINBL_tmp     = 0.0

     QV_CURR(:,1,:) = QV_CURR(:,1,:)/(1.0 - QV_CURR(:,1,:))  ! Assuming input forcing are specific hum.;
                                                             ! WRF wants mixing ratio at driver level
     P8W(:,2,:)     = P8W(:,1,:)      ! WRF uses lowest two layers
     T_PHY(:,2,:)   = T_PHY(:,1,:)    ! Only pressure is needed in two layer but fill the rest
     U_PHY(:,2,:)   = U_PHY(:,1,:)    ! 
     V_PHY(:,2,:)   = V_PHY(:,1,:)    ! 
     QV_CURR(:,2,:) = QV_CURR(:,1,:)  ! 
     RAINBL = RAINBL_tmp * DTBL       ! RAINBL in WRF is [mm]
     SR         = 0.0                 ! Will only use component if opt_snf=4
     RAINCV     = 0.0
     RAINNCV    = RAINBL
     RAINSHV    = 0.0
     SNOWNCV    = 0.0
     GRAUPELNCV = 0.0
     HAILNCV    = 0.0
     DZ8W = 2*ZLVL                    ! 2* to be consistent with WRF model level
!------------------------------------------------------------------------
! Noah-MP updates we can do before spatial loop.
!------------------------------------------------------------------------

   ! create a few fields that are IN in WRF - coszen, julian_in,yr

    DO J = YSTART,YEND
    DO I = XSTART,XEND
      CALL CALC_DECLIN(OLDDATE(1:19),XLAT_URB2D(I,J), XLONIN(I,J),COSZEN(I,J),JULIAN_IN)
    END DO
    END DO

    READ(OLDDATE(1:4),*)  YR
    YEARLEN = 365                      ! find length of year for phenology (also S Hemisphere)
    if (mod(YR,4) == 0) then
       YEARLEN = 366
       if (mod(YR,100) == 0) then
          YEARLEN = 365
          if (mod(YR,400) == 0) then
             YEARLEN = 366
          endif
       endif
    endif

    IF (ITIME == 1 .AND. .NOT. RESTART_FLAG ) THEN
      EAHXY = (P8W(:,1,:)*QV_CURR(:,1,:))/(0.622+QV_CURR(:,1,:)) ! Initial guess only.
      TAHXY = T_PHY(:,1,:)                                       ! Initial guess only.
      CHXY = 0.1
      CMXY = 0.1
    ENDIF

!------------------------------------------------------------------------
! Skip model call at t=1 since initial conditions are at start time; First model time is +1
!------------------------------------------------------------------------

   IF (ITIME > 0) THEN

!------------------------------------------------------------------------
! Call to Noah-MP driver same as surface_driver
!------------------------------------------------------------------------
     sflx_count_sum = 0 ! Timing

   ! Timing information for SFLX:

    call system_clock(count=count_before_sflx, count_rate=clock_rate)
    
    if (OPT_WT .eq. 1) then
      if (wvtflag == 2) then
        wvt_ratio = 1.0
      else
        wvt_ratio = 0.0
      end if
    end if

print *, "FLAG6"

    if (wvtflag .le. 0) then
print *, "FLAG7"
print *, "NTRACER", NTRACER
         CALL noahmplsm(ITIMESTEP,       YR, JULIAN_IN,   COSZEN, XLAT_URB2D, &
	           DZ8W,     DTBL,      DZS,     NUM_SOIL_LAYERS,         DX, &
		 IVGTYP,   ISLTYP,   VEGFRA,   GVFMAX,       TMN,             &
		  XLAND,     XICE,     XICE_THRESHOLD,                        &
                  IDVEG, IOPT_CRS, IOPT_BTR, IOPT_RUN,  IOPT_SFC,   IOPT_FRZ, &
	       IOPT_INF, IOPT_RAD, IOPT_ALB, IOPT_SNF, IOPT_TBOT,   IOPT_STC, &
	       IOPT_GLA, IOPT_RSF, IZ0TLND,                                   &
		  T_PHY,  QV_CURR,    U_PHY,    V_PHY,    SWDOWN,        GLW, &
		    P8W,   RAINBL,       SR,                                  &
		    TSK,      HFX,      QFX,       LH,    GRDFLX,     SMSTAV, &
		 SMSTOT,SFCRUNOFF, UDRUNOFF,   ALBEDO,     SNOWC,      SMOIS, &
		   SH2O,     TSLB,     SNOW,    SNOWH,    CANWAT,     ACSNOM, &
		 ACSNOW,    EMISS,     QSFC,                                  &
 		     Z0,      ZNT,                                            & ! IN/OUT LSM eqv
		ISNOWXY,     TVXY,     TGXY, CANICEXY,  CANLIQXY,      EAHXY, &
		  TAHXY,     CMXY,     CHXY,   FWETXY,  SNEQVOXY,   ALBOLDXY, &
		QSNOWXY, WSLAKEXY,    ZWTXY,     WAXY,      WTXY,     TSNOXY, &
		ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY,  RTMASSXY,   STMASSXY, &
		 WOODXY, STBLCPXY, FASTCPXY,      LAI,    XSAIXY,    TAUSSXY, &
	        SMOISEQ, SMCWTDXY,DEEPRECHXY,  RECHXY,                        & ! IN/OUT Noah MP only
	         T2MVXY,   T2MBXY,   Q2MVXY,   Q2MBXY,                        &
                 TRADXY,    NEEXY,    GPPXY,    NPPXY,    FVEGXY,    RUNSFXY, &
	        RUNSBXY,   ECANXY,   EDIRXY,  ETRANXY,     FSAXY,     FIRAXY, &
                 APARXY,    PSNXY,    SAVXY,    SAGXY,   RSSUNXY,    RSSHAXY, &
               ALBSNDXY, ALBSNIXY,                                            & ! OUT Noah MP only
                 BGAPXY,   WGAPXY,    TGVXY,    TGBXY,     CHVXY,      CHBXY, &
		  SHGXY,    SHCXY,    SHBXY,    EVGXY,     EVBXY,      GHVXY, &
		  GHBXY,    IRGXY,    IRCXY,    IRBXY,      TRXY,      EVCXY, &
	       CHLEAFXY,   CHUCXY,   CHV2XY,   CHB2XY,                        &                          
                 sfcheadrt,INFXSRT,soldrain,                          &    !O
                ids,ide, jds,jde, kds,kde,                      &
                ims,ime, jms,jme, kms,kme,                      &
                its,ite, jts,jte, kts,kte,        &
! variables below are optional
                MP_RAINC =  RAINCV, MP_RAINNC =    RAINNCV, MP_SHCV = RAINSHV,&
		MP_SNOW  = SNOWNCV, MP_GRAUP  = GRAUPELNCV, MP_HAIL = HAILNCV &
                , ACCPRCP=ACCPRCP,  ACCECAN=ACCECAN, ACCETRAN=ACCETRAN,  ACCEDIR=ACCEDIR  &
                , SOILSAT_TOP=SOILSAT_TOP, SOILSAT=SOILSAT, SOILICE=SOILICE, SNOWT_AVG=SNOWT_AVG      &
                , NTRACER=NTRACER)

print *, "NTRACER", NTRACER
!print *, "FLAG8,SFCRUNOFF(28,22):",SFCRUNOFF(28,22)
    else
!print *, "FLAG9,SFCRUNOFF(44,23),SFCRUNOFF_TR(44,23):",SFCRUNOFF(44,23),SFCRUNOFF_TR(44,23)
print *, "NTRACER", NTRACER
         CALL noahmplsm(ITIMESTEP,       YR, JULIAN_IN,   COSZEN, XLAT_URB2D, &
	           DZ8W,     DTBL,      DZS,     NUM_SOIL_LAYERS,         DX, &
		 IVGTYP,   ISLTYP,   VEGFRA,   GVFMAX,       TMN,             &
		  XLAND,     XICE,     XICE_THRESHOLD,                        &
                  IDVEG, IOPT_CRS, IOPT_BTR, IOPT_RUN,  IOPT_SFC,   IOPT_FRZ, &
	       IOPT_INF, IOPT_RAD, IOPT_ALB, IOPT_SNF, IOPT_TBOT,   IOPT_STC, &
	       IOPT_GLA, IOPT_RSF, IZ0TLND,                                   &
		  T_PHY,  QV_CURR,    U_PHY,    V_PHY,    SWDOWN,        GLW, &
		    P8W,   RAINBL,       SR,                                  &
		    TSK,      HFX,      QFX,       LH,    GRDFLX,     SMSTAV, &
		 SMSTOT,SFCRUNOFF, UDRUNOFF,   ALBEDO,     SNOWC,      SMOIS, &
		   SH2O,     TSLB,     SNOW,    SNOWH,    CANWAT,     ACSNOM, &
		 ACSNOW,    EMISS,     QSFC,                                  &
 		     Z0,      ZNT,                                            & ! IN/OUT LSM eqv
		ISNOWXY,     TVXY,     TGXY, CANICEXY,  CANLIQXY,      EAHXY, &
		  TAHXY,     CMXY,     CHXY,   FWETXY,  SNEQVOXY,   ALBOLDXY, &
		QSNOWXY, WSLAKEXY,    ZWTXY,     WAXY,      WTXY,     TSNOXY, &
		ZSNSOXY,  SNICEXY,  SNLIQXY, LFMASSXY,  RTMASSXY,   STMASSXY, &
		 WOODXY, STBLCPXY, FASTCPXY,      LAI,    XSAIXY,    TAUSSXY, &
	        SMOISEQ, SMCWTDXY,DEEPRECHXY,  RECHXY,                        & ! IN/OUT Noah MP only
	         T2MVXY,   T2MBXY,   Q2MVXY,   Q2MBXY,                        &
                 TRADXY,    NEEXY,    GPPXY,    NPPXY,    FVEGXY,    RUNSFXY, &
	        RUNSBXY,   ECANXY,   EDIRXY,  ETRANXY,     FSAXY,     FIRAXY, &
                 APARXY,    PSNXY,    SAVXY,    SAGXY,   RSSUNXY,    RSSHAXY, &
               ALBSNDXY, ALBSNIXY,                                            & ! OUT Noah MP only
                 BGAPXY,   WGAPXY,    TGVXY,    TGBXY,     CHVXY,      CHBXY, &
		  SHGXY,    SHCXY,    SHBXY,    EVGXY,     EVBXY,      GHVXY, &
		  GHBXY,    IRGXY,    IRCXY,    IRBXY,      TRXY,      EVCXY, &
	       CHLEAFXY,   CHUCXY,   CHV2XY,   CHB2XY,                        &                          
                 sfcheadrt,INFXSRT,soldrain,                          &    !O
                ids,ide, jds,jde, kds,kde,                      &
                ims,ime, jms,jme, kms,kme,                      &
                its,ite, jts,jte, kts,kte,        &
! variables below are optional
                MP_RAINC =  RAINCV, MP_RAINNC =    RAINNCV, MP_SHCV = RAINSHV,&
		MP_SNOW  = SNOWNCV, MP_GRAUP  = GRAUPELNCV, MP_HAIL = HAILNCV &
                , ACCPRCP=ACCPRCP,  ACCECAN=ACCECAN, ACCETRAN=ACCETRAN,  ACCEDIR=ACCEDIR  &
                , SOILSAT_TOP=SOILSAT_TOP, SOILSAT=SOILSAT, SOILICE=SOILICE, SNOWT_AVG=SNOWT_AVG      &
                , NTRACER=NTRACER, wvtflag=wvtflag, iopt_pm=iopt_pm, wvt_ratio=wvt_ratio, NSUB=NSUB           & ! IN WVT
		, SMOIS_SUB=SMOIS_SUB, SH2OXY_SUB=SH2OXY_SUB                                              & ! IN WVT
                , QFX_TR=QFX_TR, LH_TR=LH_TR, SFCRUNOFF_TR=SFCRUNOFF_TR, UDRUNOFF_TR=UDRUNOFF_TR      & !IN/OUT WVT
		, SMOIS_TR=SMOIS_TR, SH2OXY_TR=SH2OXY_TR, SNOWXY_TR=SNOWXY_TR, SNOWHXY_TR=SNOWHXY_TR, CANWAT_TR=CANWAT_TR &
		, ACSNOM_TR=ACSNOM_TR, ACSNOW_TR=ACSNOW_TR, CANICEXY_TR=CANICEXY_TR, CANLIQXY_TR=CANLIQXY_TR  &
		, QSNOWXY_TR=QSNOWXY_TR, WSLAKEXY_TR=WSLAKEXY_TR, WAXY_TR=WAXY_TR, WTXY_TR=WTXY_TR   &
		, ZSNSOXY_TR=ZSNSOXY_TR, SNICEXY_TR=SNICEXY_TR, SNLIQXY_TR=SNLIQXY_TR                & ! IN/OUT WVT variables
                , SMCWTDXY_TR=SMCWTDXY_TR, DEEPRECHXY_TR=DEEPRECHXY_TR, RECHXY_TR=RECHXY_TR          & ! INOUT
                , RUNSFXY_TR=RUNSFXY_TR, RUNSBXY_TR=RUNSBXY_TR, ECANXY_TR=ECANXY_TR                  &
		, EDIRXY_TR=EDIRXY_TR, ETRANXY_TR=ETRANXY_TR, FLUXXY_TR=FLUXXY_TR                    & ! IN WVT variables
                , INFXSRT_TR=INFXSRT_TR, sfcheadrt_tr=sfcheadrt_tr,soldrain_tr=soldrain_tr           &    !O water tracer
                )
!print *, "FLAG9-1, SMOIS(44,23,:),SMOIS_SUB(44,23,:),SH2OX(44,23),SH2OXY_SUB(44,23):",SMOIS(44,:,23),SMOIS_SUB(44,:,23),SH2O(44,:,23),SH2OXY_SUB(44,:,23)
!print *, "FLAG9-1,SFCRUNOFF(32,26),SFCRUNOFF_TR(32,26),RUNSFXY(32,26),RUNSFXY_TR(32,26),INFXSRT(32,26),INFXSRT_TR(32,26),sfcheadrt(32,26),sfcheadrt_tr(32,26):",SFCRUNOFF(32,26),SFCRUNOFF_TR(32,26),RUNSFXY(32,26),RUNSFXY_TR(32,26),INFXSRT(32,26),INFXSRT_TR(32,26),sfcheadrt(32,26),sfcheadrt_tr(32,26)

    end if

          call system_clock(count=count_after_sflx, count_rate=clock_rate)
          sflx_count_sum = sflx_count_sum + ( count_after_sflx - count_before_sflx )

  IF(RUNOFF_OPTION.EQ.5.AND.MOD(ITIME,STEPWTD).EQ.0)THEN
           CALL wrf_message('calling WTABLE' )

!gmm update wtable from lateral flow and shed water to rivers
           CALL WTABLE_MMF_NOAHMP(                                        &
	       NUM_SOIL_LAYERS,  XLAND, XICE,       XICE_THRESHOLD, ISICE,    &
               ISLTYP,      SMOISEQ,    DZS,        WTDDT,                &
               FDEPTHXY,    AREAXY,     TERRAIN,    ISURBAN,    IVGTYP,   &
               RIVERCONDXY, RIVERBEDXY, EQZWT,      PEXPXY,               &
               SMOIS,       SH2O,       SMCWTDXY,   ZWTXY,                &
	       QRFXY,       DEEPRECHXY, QSPRINGXY,                        &
               QSLATXY,     QRFSXY,     QSPRINGSXY, RECHXY,               &
               IDS,IDE, JDS,JDE, KDS,KDE,                                 &
               IMS,IME, JMS,JME, KMS,KME,                                 &
               ITS,ITE, JTS,JTE, KTS,KTE )

 ENDIF

!------------------------------------------------------------------------
! END of surface_driver consistent code
!------------------------------------------------------------------------

 ENDIF   ! SKIP FIRST TIMESTEP

!print *, "FLAG9-2,SFCRUNOFF(28,22):",SFCRUNOFF(28,22)
     call geth_newdate(newdate, olddate, nint(dtbl))
     olddate = newdate
     if (OPT_WT .le. 0) then
       call hrldas_drv_HYDRO(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:),infxsrt,sfcheadrt,soldrain,ix,jx,NSOIL,NTRACER) 
     else
       call hrldas_drv_HYDRO(TSLB(:,1:NSOIL,:),SMOIS(:,1:NSOIL,:),SH2O(:,1:NSOIL,:),infxsrt,sfcheadrt,soldrain,ix,jx,NSOIL, &
            NTRACER, WVTFLAG=WVTFLAG, SMOIS_SUB=SMOIS_SUB, SH2OXY_SUB=SH2OXY_SUB, &    ! water tracer vars      
            SMOIS_TR=SMOIS_TR, SH2OXY_TR=SH2OXY_TR,              &    ! water tracer vars
            infxsrt_tr=infxsrt_tr, sfcheadrt_tr=sfcheadrt_tr,soldrain_tr=soldrain_tr)  ! water tracer vars
!print *, "FLAG9-3, SMOIS(44,:,23),SMOIS_SUB(44,:,23),SH2OX(44,:,23),SH2OXY_SUB(44,:,23):",SMOIS(44,:,23),SMOIS_SUB(44,:,23),SH2O(44,:,23),SH2OXY_SUB(44,:,23)
     end if

! Output for history
     OUTPUT_FOR_HISTORY: if (output_timestep > 0) then
        if (mod(ITIME*noah_timestep, output_timestep) == 0) then

           ! convert RAINRATE back to mm/s for output
           RAINBL = RAINBL_tmp
!print *, "FLAG9-3,SFCRUNOFF(28,22):",SFCRUNOFF(28,22)
           call ldas_output(ITIME, wvtflag)


        endif
     endif OUTPUT_FOR_HISTORY

!print *, "FLAG9-4"
     if (IVGTYP(xstart,ystart)==ISWATER) then
       write(*,'(" ***DATE=", A19)', advance="NO") olddate
     else
       write(*,'(" ***DATE=", A19, 6F10.5)', advance="NO") olddate, TSLB(xstart,1,ystart), LAI(xstart,ystart)
     endif

!------------------------------------------------------------------------
! Write Restart - timestamp equal to output will have same states
!------------------------------------------------------------------------

      if ( (restart_frequency_hours .gt. 0) .and. &
           (mod(ITIME, int(restart_frequency_hours*3600./nint(dtbl))) == 0)) then
!print *, "FLAG9-5"
       if(rst_bi_out .eq. 0) then
           call lsm_restart()
       else
           call lsm_rst_bi_out()
       endif
      else
!print *, "FLAG9-6"
       if (restart_frequency_hours <= 0) then
          if ( (olddate( 9:10) == "01") .and. (olddate(12:13) == "00") .and. &
               (olddate(15:16) == "00") .and. (olddate(18:19) == "00") ) then
               if(rst_bi_out .eq. 0) then
                   call lsm_restart()  ! jlm - i moved all the restart code to a subroutine. 
               else
                   call lsm_rst_bi_out()
               endif
          endif
       endif
      endif

!------------------------------------------------------------------------
! Advance the time 
!------------------------------------------------------------------------

! update the timer
     call system_clock(count=clock_count_2, count_rate=clock_rate)
     timing_sum = timing_sum + float(clock_count_2-clock_count_1)/float(clock_rate)
     write(*,'("    Timing: ",f6.2," Cumulative:  ", f10.2, "  SFLX: ", f6.2 )') &
          float(clock_count_2-clock_count_1)/float(clock_rate), &
          timing_sum, real(sflx_count_sum) / real(clock_rate)
     clock_count_1 = clock_count_2


      ITIMESTEP = ITIMESTEP + 1
!print *, "FLAG10,SFCRUNOFF(28,22):",SFCRUNOFF(28,22)

end subroutine land_driver_exe

!!===============================================================================
subroutine  ldas_output(itime, wvtflag)
integer, intent(in)  :: itime ! time step of the LSM  
integer, intent(in)  :: wvtflag ! -999-no tracer capability, 0-2-with tracer capability

!#ifdef 1
!if ( (io_config_outputs .eq. 0) ) then
!#endif
!#ifndef 1
!           call prepare_output_file (trim(outdir), version, &
!                igrid, output_timestep, llanduse, split_output_count, hgrid,                &
!                ixfull, jxfull, ixpar, jxpar, xstartpar, ystartpar,                         &
!                iswater, mapproj, lat1, lon1, dx, dy, truelat1, truelat2, cen_lon,          &
!                nsoil, nsnow, dzs, startdate, olddate, IVGTYP, ISLTYP)
!
!           DEFINE_MODE_LOOP : do imode = 1, 2
!
!              call set_output_define_mode(imode)
!
!              ! For 3D arrays, we need to know whether the Z dimension is snow layers, or soil layers.
!
!        ! Properties - Assigned or predicted
!              call add_to_output(IVGTYP     , "IVGTYP"  , "Dominant vegetation category"         , "category"              )
!              call add_to_output(ISLTYP     , "ISLTYP"  , "Dominant soil category"               , "category"              )
!              call add_to_output(FVEGXY     , "FVEG"    , "Green Vegetation Fraction"              , "-"                   )
!              call add_to_output(LAI        , "LAI"     , "Leaf area index"                      , "-"                     )
!              call add_to_output(XSAIXY     , "SAI"     , "Stem area index"                      , "-"                     )
!        ! Forcing
!              call add_to_output(SWDOWN     , "SWFORC"  , "Shortwave forcing"                    , "W m{-2}"               )
!              call add_to_output(COSZEN     , "COSZ"    , "Cosine of zenith angle"                    , "W m{-2}"               )
!              call add_to_output(GLW        , "LWFORC"  , "Longwave forcing"                    , "W m{-2}"               )
!              call add_to_output(RAINBL     , "RAINRATE", "Precipitation rate"                   , "kg m{-2} s{-1}"        )
!        ! Grid energy budget terms
!              call add_to_output(EMISS      , "EMISS"   , "Grid emissivity"                    , ""               )
!              call add_to_output(FSAXY      , "FSA"     , "Total absorbed SW radiation"          , "W m{-2}"               )         
!              call add_to_output(FIRAXY     , "FIRA"    , "Total net LW radiation to atmosphere" , "W m{-2}"               )
!              call add_to_output(GRDFLX     , "GRDFLX"  , "Heat flux into the soil"              , "W m{-2}"               )
!              call add_to_output(HFX        , "HFX"     , "Total sensible heat to atmosphere"    , "W m{-2}"               )
!              call add_to_output(LH         , "LH"      , "Total latent heat to atmosphere"    , "W m{-2}"               )
!              call add_to_output(ECANXY     , "ECAN"    , "Canopy water evaporation rate"        , "kg m{-2} s{-1}"        )
!              call add_to_output(EDIRXY     , "EDIR"    , "Direct from soil evaporation rate"    , "kg m{-2} s{-1}"        )
!              call add_to_output(ALBEDO     , "ALBEDO"  , "Surface albedo"                         , "-"                   )
!              call add_to_output(ETRANXY    , "ETRAN"   , "Transpiration rate"                   , "kg m{-2} s{-1}"        )
!        ! Grid water budget terms - in addition to above
!              call add_to_output(UDRUNOFF   , "UGDRNOFF", "Accumulated underground runoff"       , "mm"                    )
!              call add_to_output(SFCRUNOFF  , "SFCRNOFF", "Accumulatetd surface runoff"          , "mm"                    )
!              call add_to_output(CANLIQXY   , "CANLIQ"  , "Canopy liquid water content"          , "mm"                    )
!              call add_to_output(CANICEXY   , "CANICE"  , "Canopy ice water content"             , "mm"                    )
!              call add_to_output(ZWTXY      , "ZWT"     , "Depth to water table"                 , "m"                     )
!              call add_to_output(WAXY       , "WA"      , "Water in aquifer"                     , "kg m{-2}"              )
!              call add_to_output(WTXY       , "WT"      , "Water in aquifer and saturated soil"  , "kg m{-2}"              )
!        ! Additional needed to close the canopy energy budget
!              call add_to_output(SAVXY      , "SAV"     , "Solar radiative heat flux absorbed by vegetation", "W m{-2}"    )
!              call add_to_output(TRXY       , "TR"      , "Transpiration heat"                     , "W m{-2}"             )
!              call add_to_output(EVCXY      , "EVC"     , "Canopy evap heat"                       , "W m{-2}"             )
!              call add_to_output(IRCXY      , "IRC"     , "Canopy net LW rad"                      , "W m{-2}"             )
!              call add_to_output(SHCXY      , "SHC"     , "Canopy sensible heat"                   , "W m{-2}"             )
!        ! Additional needed to close the under canopy ground energy budget
!              call add_to_output(IRGXY      , "IRG"     , "Ground net LW rad"                      , "W m{-2}"             )
!              call add_to_output(SHGXY      , "SHG"     , "Ground sensible heat"                   , "W m{-2}"             )
!              call add_to_output(EVGXY      , "EVG"     , "Ground evap heat"                       , "W m{-2}"             )
!              call add_to_output(GHVXY      , "GHV"     , "Ground heat flux + to soil vegetated"   , "W m{-2}"             )
!        ! Needed to close the bare ground energy budget
!              call add_to_output(SAGXY      , "SAG"     , "Solar radiative heat flux absorbed by ground", "W m{-2}"        )
!              call add_to_output(IRBXY      , "IRB"     , "Net LW rad to atm bare"                 , "W m{-2}"             )
!              call add_to_output(SHBXY      , "SHB"     , "Sensible heat to atm bare"              , "W m{-2}"             )
!              call add_to_output(EVBXY      , "EVB"     , "Evaporation heat to atm bare"           , "W m{-2}"             )
!              call add_to_output(GHBXY      , "GHB"     , "Ground heat flux + to soil bare"        , "W m{-2}"             )
!        ! Above-soil temperatures
!              call add_to_output(TRADXY     , "TRAD"    , "Surface radiative temperature"        , "K"                     )
!              call add_to_output(TGXY       , "TG"      , "Ground temperature"                   , "K"                     )
!              call add_to_output(TVXY       , "TV"      , "Vegetation temperature"               , "K"                     )
!              call add_to_output(TAHXY      , "TAH"     , "Canopy air temperature"               , "K"                     )
!              call add_to_output(TGVXY      , "TGV"     , "Ground surface Temp vegetated"          , "K"                   )
!              call add_to_output(TGBXY      , "TGB"     , "Ground surface Temp bare"               , "K"                   )
!              call add_to_output(T2MVXY     , "T2MV"    , "2m Air Temp vegetated"                  , "K"                   )
!              call add_to_output(T2MBXY     , "T2MB"    , "2m Air Temp bare"                       , "K"                   )
!        ! Above-soil moisture
!              call add_to_output(Q2MVXY     , "Q2MV"    , "2m mixing ratio vegetated"              , "kg/kg"               )
!              call add_to_output(Q2MBXY     , "Q2MB"    , "2m mixing ratio bare"                   , "kg/kg"               )
!              call add_to_output(EAHXY      , "EAH"     , "Canopy air vapor pressure"            , "Pa"                    )
!              call add_to_output(FWETXY     , "FWET"    , "Wetted or snowed fraction of canopy"  , "fraction"              )
!        ! Snow and soil - 3D terms
!              call add_to_output(ZSNSOXY(:,-nsnow+1:0,:),  "ZSNSO_SN" , "Snow layer depths from snow surface", "m", "SNOW")
!              call add_to_output(SNICEXY    , "SNICE"   , "Snow layer ice"                       , "mm"             , "SNOW")
!              call add_to_output(SNLIQXY    , "SNLIQ"   , "Snow layer liquid water"              , "mm"             , "SNOW")
!              call add_to_output(TSLB       , "SOIL_T"  , "soil temperature"                     , "K"              , "SOIL")
!              call add_to_output(SH2O       , "SOIL_W"  , "liquid volumetric soil moisture"      , "m3 m-3"         , "SOIL")
!              call add_to_output(TSNOXY     , "SNOW_T"  , "snow temperature"                     , "K"              , "SNOW")
!              call add_to_output(SMOIS      , "SOIL_M"  , "volumetric soil moisture"             , "m{3} m{-3}"     , "SOIL")
!        ! Snow - 2D terms
!              call add_to_output(SNOWH      , "SNOWH"   , "Snow depth"                           , "m"                     )
!              call add_to_output(SNOW       , "SNEQV"   , "Snow water equivalent"                , "kg m{-2}"              )
!              call add_to_output(QSNOWXY    , "QSNOW"   , "Snowfall rate"                        , "mm s{-1}"              )
!              call add_to_output(ISNOWXY    , "ISNOW"   , "Number of snow layers"                , "count"                 )
!              call add_to_output(SNOWC      , "FSNO"    , "Snow-cover fraction on the ground"      , ""                    )
!              call add_to_output(ACSNOW     , "ACSNOW"  , "accumulated snow fall"                  , "mm"                  )
!              call add_to_output(ACSNOM     , "ACSNOM"  , "accumulated melting water out of snow bottom" , "mm"            )
!        ! Exchange coefficients
!              call add_to_output(CMXY       , "CM"      , "Momentum drag coefficient"            , ""                      )
!              call add_to_output(CHXY       , "CH"      , "Sensible heat exchange coefficient"   , ""                      )
!              call add_to_output(CHVXY      , "CHV"     , "Exchange coefficient vegetated"         , "m s{-1}"             )
!              call add_to_output(CHBXY      , "CHB"     , "Exchange coefficient bare"              , "m s{-1}"             )
!              call add_to_output(CHLEAFXY   , "CHLEAF"  , "Exchange coefficient leaf"              , "m s{-1}"             )
!              call add_to_output(CHUCXY     , "CHUC"    , "Exchange coefficient bare"              , "m s{-1}"             )
!              call add_to_output(CHV2XY     , "CHV2"    , "Exchange coefficient 2-meter vegetated" , "m s{-1}"             )
!              call add_to_output(CHB2XY     , "CHB2"    , "Exchange coefficient 2-meter bare"      , "m s{-1}"             )
!        ! Carbon allocation model
!              call add_to_output(LFMASSXY   , "LFMASS"  , "Leaf mass"                            , "g m{-2}"               )
!              call add_to_output(RTMASSXY   , "RTMASS"  , "Mass of fine roots"                   , "g m{-2}"               )
!              call add_to_output(STMASSXY   , "STMASS"  , "Stem mass"                            , "g m{-2}"               )
!              call add_to_output(WOODXY     , "WOOD"    , "Mass of wood and woody roots"         , "g m{-2}"               )
!              call add_to_output(STBLCPXY   , "STBLCP"  , "Stable carbon in deep soil"           , "g m{-2}"               )
!              call add_to_output(FASTCPXY   , "FASTCP"  , "Short-lived carbon in shallow soil"   , "g m{-2}"               )
!              call add_to_output(NEEXY      , "NEE"     , "Net ecosystem exchange"                 , "g m{-2} s{-1} CO2"   )
!              call add_to_output(GPPXY      , "GPP"     , "Net instantaneous assimilation"         , "g m{-2} s{-1} C"     )
!              call add_to_output(NPPXY      , "NPP"     , "Net primary productivity"               , "g m{-2} s{-1} C"     )
!              call add_to_output(PSNXY      , "PSN"     , "Total photosynthesis"                   , "umol CO@ m{-2} s{-1}")
!              call add_to_output(APARXY     , "APAR"    , "Photosynthesis active energy by canopy" , "W m{-2}"             )
!
!        ! Carbon allocation model
!	    IF(RUNOFF_OPTION == 5) THEN
!              call add_to_output(SMCWTDXY   , "SMCWTD"   , "Leaf mass"                            , "g m{-2}"               )
!              call add_to_output(RECHXY     , "RECH"     , "Mass of fine roots"                   , "g m{-2}"               )
!              call add_to_output(QRFSXY     , "QRFS"     , "Stem mass"                            , "g m{-2}"               )
!              call add_to_output(QSPRINGSXY , "QSPRINGS" , "Mass of wood and woody roots"         , "g m{-2}"               )
!              call add_to_output(QSLATXY    , "QSLAT"    , "Stable carbon in deep soil"           , "g m{-2}"               )
!	    ENDIF
!
!           enddo DEFINE_MODE_LOOP
!
!           call finalize_output_file(split_output_count)
!#ifdef 1
!else
!#endif
!#endif 

   ! Logan add begin
   ! Go through each variable. For the first time the NWM output routine is
   ! called, the file is created and all output variables (desired per flags),
   ! are created in define mode. 

   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,ITIME,startdate,olddate,ixpar,jxpar,1,float(IVGTYP),IVGTYP,1)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,float(ISLTYP),IVGTYP,2)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,FVEGXY,IVGTYP,3)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,LAI,IVGTYP,4)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,XSAIXY,IVGTYP,5)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SWDOWN,IVGTYP,6)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,COSZEN,IVGTYP,7)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,GLW,IVGTYP,8)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,RAINBL,IVGTYP,9)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,EMISS,IVGTYP,10)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,FSAXY,IVGTYP,11)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,FIRAXY,IVGTYP,12)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,GRDFLX,IVGTYP,13)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,HFX,IVGTYP,14)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,LH,IVGTYP,15)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ECANXY,IVGTYP,16)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,EDIRXY,IVGTYP,17)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ALBEDO,IVGTYP,18)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ETRANXY,IVGTYP,19)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,UDRUNOFF,IVGTYP,20)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SFCRUNOFF,IVGTYP,21)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CANLIQXY,IVGTYP,22)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CANICEXY,IVGTYP,23)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ZWTXY,IVGTYP,24)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,WAXY,IVGTYP,25)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,WTXY,IVGTYP,26)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACCPRCP,IVGTYP,27)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACCECAN,IVGTYP,28)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACCEDIR,IVGTYP,29)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACCETRAN,IVGTYP,30)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SAVXY,IVGTYP,31)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,TRXY,IVGTYP,32)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,EVCXY,IVGTYP,33)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,IRCXY,IVGTYP,34)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SHCXY,IVGTYP,35)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,IRGXY,IVGTYP,36)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SHGXY,IVGTYP,37)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,EVGXY,IVGTYP,38)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,GHVXY,IVGTYP,39)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SAGXY,IVGTYP,40)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,IRBXY,IVGTYP,41)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SHBXY,IVGTYP,42)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,EVBXY,IVGTYP,43)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,GHBXY,IVGTYP,44)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,TRADXY,IVGTYP,45)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,TGXY,IVGTYP,46)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,TVXY,IVGTYP,47)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,TAHXY,IVGTYP,48)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,TGVXY,IVGTYP,49)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,TGBXY,IVGTYP,50)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,T2MVXY,IVGTYP,51)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,T2MBXY,IVGTYP,52)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,Q2MVXY,IVGTYP,53)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,Q2MBXY,IVGTYP,54)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,EAHXY,IVGTYP,55)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,FWETXY,IVGTYP,56)
   !call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,float(ZSNOXY(:,-nsnow+1:0,:)),IVGTYP,57)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,3,SNICEXY,IVGTYP,58)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,3,SNLIQXY,IVGTYP,59)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,4,TSLB,IVGTYP,60)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,4,SH2O,IVGTYP,61)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,3,TSNOXY,IVGTYP,62)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,4,SMOIS,IVGTYP,63)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SNOWH,IVGTYP,64)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SNOW,IVGTYP,65)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,QSNOWXY,IVGTYP,66)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,float(ISNOWXY),IVGTYP,67)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SNOWC,IVGTYP,68)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACSNOW,IVGTYP,69)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACSNOM,IVGTYP,70)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CMXY,IVGTYP,71)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CHXY,IVGTYP,72)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CHVXY,IVGTYP,73)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CHBXY,IVGTYP,74)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CHLEAFXY,IVGTYP,75)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CHUCXY,IVGTYP,76)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CHV2XY,IVGTYP,77)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CHB2XY,IVGTYP,78)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,LFMASSXY,IVGTYP,79)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,RTMASSXY,IVGTYP,80)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,STMASSXY,IVGTYP,81)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,WOODXY,IVGTYP,82)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,STBLCPXY,IVGTYP,83)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,FASTCPXY,IVGTYP,84)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,NEEXY,IVGTYP,85)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,GPPXY,IVGTYP,86)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,NPPXY,IVGTYP,87)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,PSNXY,IVGTYP,88)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,APARXY,IVGTYP,89)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,(ACCECAN+ACCEDIR+ACCETRAN),IVGTYP,90)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,(CANLIQXY+CANICEXY),IVGTYP,91)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,(SOILICE),IVGTYP,92)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,(SOILSAT_TOP),IVGTYP,93)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,(SOILSAT),IVGTYP,94)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,(SNOWT_AVG),IVGTYP,95)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,2,ALBSNDXY,IVGTYP,96)
   call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,2,ALBSNIXY,IVGTYP,97)
   
   if (wvtflag .ge. 0) then
print *, "FLAG009-4"
print *, "wvtflag", wvtflag
print *, "LH_TR, LH:", LH_TR(0,0),LH(0,0)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,LH_TR,IVGTYP,98)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ECANXY_TR,IVGTYP,99)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,EDIRXY_TR,IVGTYP,100)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ETRANXY_TR,IVGTYP,101)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,UDRUNOFF_TR,IVGTYP,102)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SFCRUNOFF_TR,IVGTYP,103)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CANLIQXY_TR,IVGTYP,104)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,CANICEXY_TR,IVGTYP,105)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,WAXY_TR,IVGTYP,106)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,WTXY_TR,IVGTYP,107)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,3,SNICEXY_TR,IVGTYP,108)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,3,SNLIQXY_TR,IVGTYP,109)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,4,SH2OXY_TR,IVGTYP,110)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,4,SMOIS_TR,IVGTYP,111)
!     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,4,FLUXXY_TR,IVGTYP,112)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SNOWHXY_TR,IVGTYP,113)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,SNOWXY_TR,IVGTYP,114)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,QSNOWXY_TR,IVGTYP,115)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACSNOW_TR,IVGTYP,116)
     call output_NoahMP_NWM(trim(outdir),igrid,output_timestep,itime,startdate,olddate,ixpar,jxpar,1,ACSNOM_TR,IVGTYP,117)
   end if
   

!#ifdef 1
!endif
!#endif

end subroutine  ldas_output

!!===============================================================================
subroutine lsm_restart()
  implicit none 
  character(len=256):: tmpStr
  integer :: ncid
  
  print*, 'Write restart at '//olddate(1:13)

  call prepare_restart_file (trim(outdir), version, igrid, llanduse, olddate, startdate,                         & 
       ixfull, jxfull, ixpar, jxpar, xstartpar, ystartpar,                   &
       nsoil, nsnow, dx, dy, truelat1, truelat2, mapproj, lat1, lon1,        &
       cen_lon, iswater, ivgtyp, NTRACER)

   write(tmpStr, '(A,"/RESTART.",A10,"_DOMAIN",I1)') trim(outdir), olddate(1:4)//olddate(6:7)//olddate(9:10)//olddate(12:13), igrid


  if(my_id .eq. io_id) then
     ierr = nf90_open(trim(tmpStr),  NF90_WRITE, ncid)
print *, 'flag rst 1'
     call define_rst_variables(ncid)
print *, 'flag rst 3'
  endif
  
  call add_to_restart(ncid,TSLB      , "SOIL_T", layers="SOIL")
  call add_to_restart(ncid,TSNOXY    , "SNOW_T", layers="SNOW")
  call add_to_restart(ncid,SMOIS     , "SMC"   , layers="SOIL")
  call add_to_restart(ncid,SH2O      , "SH2O"  , layers="SOIL")
  call add_to_restart(ncid,ZSNSOXY   , "ZSNSO" , layers="SOSN")
  call add_to_restart(ncid,SNICEXY   , "SNICE" , layers="SNOW")
  call add_to_restart(ncid,SNLIQXY   , "SNLIQ" , layers="SNOW")
  call add_to_restart(ncid,QSNOWXY   , "QSNOW" )
  call add_to_restart(ncid,FWETXY    , "FWET"  )
  call add_to_restart(ncid,SNEQVOXY  , "SNEQVO")
  call add_to_restart(ncid,EAHXY     , "EAH"   )
  call add_to_restart(ncid,TAHXY     , "TAH"   )
  call add_to_restart(ncid,ALBOLDXY  , "ALBOLD")
  call add_to_restart(ncid,CMXY      , "CM"    )
  call add_to_restart(ncid,CHXY      , "CH"    )
  call add_to_restart(ncid,ISNOWXY   , "ISNOW" )
  call add_to_restart(ncid,CANLIQXY  , "CANLIQ")
  call add_to_restart(ncid,CANICEXY  , "CANICE")
  call add_to_restart(ncid,SNOW      , "SNEQV" )
  call add_to_restart(ncid,SNOWH     , "SNOWH" )
  call add_to_restart(ncid,TVXY      , "TV"    )
  call add_to_restart(ncid,TGXY      , "TG"    )
  call add_to_restart(ncid,ZWTXY     , "ZWT"   )
  call add_to_restart(ncid,WAXY      , "WA"    )
  call add_to_restart(ncid,WTXY      , "WT"    )
  call add_to_restart(ncid,WSLAKEXY  , "WSLAKE")
  call add_to_restart(ncid,LFMASSXY  , "LFMASS")
  call add_to_restart(ncid,RTMASSXY  , "RTMASS")
  call add_to_restart(ncid,STMASSXY  , "STMASS")
  call add_to_restart(ncid,WOODXY    , "WOOD"  )
  call add_to_restart(ncid,STBLCPXY  , "STBLCP")
  call add_to_restart(ncid,FASTCPXY  , "FASTCP")
  call add_to_restart(ncid,LAI       , "LAI"   )
  call add_to_restart(ncid,XSAIXY    , "SAI"   )
  call add_to_restart(ncid,VEGFRA    , "VEGFRA")
  call add_to_restart(ncid,GVFMIN    , "GVFMIN")
  call add_to_restart(ncid,GVFMAX    , "GVFMAX")
  call add_to_restart(ncid,ACSNOM    , "ACMELT")
  call add_to_restart(ncid,ACSNOW    , "ACSNOW")
  call add_to_restart(ncid,TAUSSXY   , "TAUSS" )
  call add_to_restart(ncid,QSFC      , "QSFC"  )
  call add_to_restart(ncid,SFCRUNOFF , "SFCRUNOFF")
  call add_to_restart(ncid,UDRUNOFF  , "UDRUNOFF" )
!-----------------------------------------------------------------
! Needed for water tracer capability
!-----------------------------------------------------------------
  call add_to_restart(ncid,SMOIS_TR     , "SMC_TR"   , layers="TRAC")   !1
  call add_to_restart(ncid,SH2OXY_TR    , "SH2O_TR"  , layers="TRAC")   !2
  call add_to_restart(ncid,SNICEXY_TR   , "SNICE_TR" , layers="SNOW")   !3
  call add_to_restart(ncid,SNLIQXY_TR   , "SNLIQ_TR" , layers="SNOW")   !4
  call add_to_restart(ncid,QSNOWXY_TR   , "QSNOW_TR" )                  !5
  call add_to_restart(ncid,CANLIQXY_TR  , "CANLIQ_TR")                  !6 
  call add_to_restart(ncid,CANICEXY_TR  , "CANICE_TR")                  !7
  call add_to_restart(ncid,SNOWXY_TR    , "SNEQV_TR" )                  !8
  call add_to_restart(ncid,SNOWHXY_TR   , "SNOWH_TR" )                  !9
  call add_to_restart(ncid,WAXY_TR      , "WA_TR"    )                  !10
  call add_to_restart(ncid,WTXY_TR      , "WT_TR"    )                  !11
  call add_to_restart(ncid,ACSNOM_TR    , "ACMELT_TR")                  !12
  call add_to_restart(ncid,ACSNOW_TR    , "ACSNOW_TR")                  !13
  call add_to_restart(ncid,SFCRUNOFF_TR , "SFCRUNOFF_TR")               !14
  call add_to_restart(ncid,UDRUNOFF_TR  , "UDRUNOFF_TR" )               !15
  call add_to_restart(ncid,SMOIS_SUB    , "SMC_SUB"  , layers="TRAC")   !16
  call add_to_restart(ncid,SH2OXY_SUB   , "SH2O_SUB" , layers="TRAC")   !17
!-----------------------------------------------------------------

  call add_to_restart(ncid,ACCPRCP   , "ACCPRCP" )
  call add_to_restart(ncid,ACCECAN   , "ACCECAN" )
  call add_to_restart(ncid,ACCEDIR   , "ACCEDIR" )
  call add_to_restart(ncid,ACCETRAN  , "ACCETRAN" )
! below for opt_run = 5
  call add_to_restart(ncid,SMOISEQ   , "SMOISEQ"  , layers="SOIL"  )
  call add_to_restart(ncid,AREAXY    , "AREAXY"     )
  call add_to_restart(ncid,SMCWTDXY  , "SMCWTDXY"   )
  call add_to_restart(ncid,DEEPRECHXY, "DEEPRECHXY" )
  call add_to_restart(ncid,QSLATXY   , "QSLATXY"    )
  call add_to_restart(ncid,QRFSXY    , "QRFSXY"     )
  call add_to_restart(ncid,QSPRINGSXY, "QSPRINGSXY" )
  call add_to_restart(ncid,RECHXY    , "RECHXY"     )
  call add_to_restart(ncid,QRFXY     , "QRFXY"      )
  call add_to_restart(ncid,QSPRINGXY , "QSPRINGXY"  )
  call add_to_restart(ncid,FDEPTHXY , "FDEPTHXY"  )
  call add_to_restart(ncid,RIVERCONDXY , "RIVERCONDXY"  )
  call add_to_restart(ncid,RIVERBEDXY , "RIVERBEDXY"  )
  call add_to_restart(ncid,EQZWT , "EQZWT"  )
  call add_to_restart(ncid,PEXPXY , "PEXPXY"  )

  if(my_id .eq. io_id) then
     ierr = nf90_close(ncid)
  endif

  call finalize_restart_file()

end subroutine lsm_restart

subroutine lsm_rst_bi_out()
  implicit none 
  integer :: iunit, ierr
  character(len=256) :: output_flnm, str_tmp
  integer  :: i0,ie, i, istep, mkdirStatus


  call mpp_land_sync()


 i0 = 0
 istep = 64
 ie = istep
 do i = 0, numprocs,istep
   if(my_id .ge. i0 .and. my_id .lt. ie) then

  write(output_flnm, '(A,"/RESTART.",A10,"_DOMAIN",I1)') trim(outdir), olddate(1:4)//olddate(6:7)//olddate(9:10)//olddate(12:13),igrid
  iunit =56 

             if(my_id .lt. 10) then
                  write(str_tmp,'(I1)') my_id
             else if(my_id .lt. 100) then
                  write(str_tmp,'(I2)') my_id
             else if(my_id .lt. 1000) then
                  write(str_tmp,'(I3)') my_id
             else if(my_id .lt. 10000) then
                  write(str_tmp,'(I4)') my_id
             else if(my_id .lt. 100000) then
                  write(str_tmp,'(I5)') my_id
             else
                continue
             endif
  open(iunit,file="restart/"//trim(output_flnm)//"."//trim(str_tmp),form="unformatted",ERR=102, access="sequential")

  write(iunit,ERR=101) TSLB  
  write(iunit,ERR=101) TSNOXY
  write(iunit,ERR=101) SMOIS     
  write(iunit,ERR=101) SH2O      
  write(iunit,ERR=101) ZSNSOXY   
  write(iunit,ERR=101) SNICEXY   
  write(iunit,ERR=101) SNLIQXY   
  write(iunit,ERR=101) QSNOWXY   
  write(iunit,ERR=101) FWETXY    
  write(iunit,ERR=101) SNEQVOXY  
  write(iunit,ERR=101) EAHXY     
  write(iunit,ERR=101) TAHXY     
  write(iunit,ERR=101) ALBOLDXY  
  write(iunit,ERR=101) CMXY      
  write(iunit,ERR=101) CHXY      
  write(iunit,ERR=101) ISNOWXY   
  write(iunit,ERR=101) CANLIQXY  
  write(iunit,ERR=101) CANICEXY  
  write(iunit,ERR=101) SNOW      
  write(iunit,ERR=101) SNOWH     
  write(iunit,ERR=101) TVXY      
  write(iunit,ERR=101) TGXY      
  write(iunit,ERR=101) ZWTXY     
  write(iunit,ERR=101) WAXY      
  write(iunit,ERR=101) WTXY      
  write(iunit,ERR=101) WSLAKEXY  
  write(iunit,ERR=101) LFMASSXY  
  write(iunit,ERR=101) RTMASSXY  
  write(iunit,ERR=101) STMASSXY  
  write(iunit,ERR=101) WOODXY    
  write(iunit,ERR=101) STBLCPXY  
  write(iunit,ERR=101) FASTCPXY  
  write(iunit,ERR=101) LAI       
  write(iunit,ERR=101) XSAIXY    
  write(iunit,ERR=101) VEGFRA    
  write(iunit,ERR=101) GVFMIN    
  write(iunit,ERR=101) GVFMAX    
  write(iunit,ERR=101) ACSNOM    
  write(iunit,ERR=101) ACSNOW    
  write(iunit,ERR=101) TAUSSXY   
  write(iunit,ERR=101) QSFC      
  write(iunit,ERR=101) SFCRUNOFF    
  write(iunit,ERR=101) UDRUNOFF     
! #ifndef REALTIME
! #ifdef 1
!   write(iunit,ERR=101) ACCPRCP   
!   write(iunit,ERR=101) ACCECAN   
!   write(iunit,ERR=101) ACCEDIR   
!   write(iunit,ERR=101) ACCETRAN  
! #endif
! #endif
! ! below for opt_run = 5
!   if(IOPT_RUN .eq. 5) then
!      write(iunit,ERR=101) SMOISEQ   
!      write(iunit,ERR=101) AREAXY    
!      write(iunit,ERR=101) SMCWTDXY  
!      write(iunit,ERR=101) DEEPRECHXY
!      write(iunit,ERR=101) QSLATXY   
!      write(iunit,ERR=101) QRFSXY    
!      write(iunit,ERR=101) QSPRINGSXY
!      write(iunit,ERR=101) RECHXY    
!      write(iunit,ERR=101) QRFXY     
!      write(iunit,ERR=101) QSPRINGXY 
!      write(iunit,ERR=101) FDEPTHXY 
!      write(iunit,ERR=101) RIVERCONDXY 
!      write(iunit,ERR=101) RIVERBEDXY 
!      write(iunit,ERR=101) EQZWT 
!      write(iunit,ERR=101) PEXPXY 
!   endif

  close(iunit)

    endif
    call mpp_land_sync()
    i0 = i0 + istep
    ie = ie + istep
  end do ! end do of i loop

  return
101  continue
  call fatal_error_stop("FATAL ERROR: failed to write lsm restartfile")
102  continue
  call fatal_error_stop("FATAL ERROR: failed to open lsm restartfile")
end subroutine lsm_rst_bi_out

subroutine lsm_rst_bi_in()
  implicit none 
  integer :: iunit, ierr
  character(len=256):: str_tmp
  integer  :: i0,ie, i, istep

  iunit = 56
 i0 = 0
 istep = 64
 ie = istep
 do i = 0, numprocs,istep
   if(my_id .ge. i0 .and. my_id .lt. ie) then

             if(my_id .lt. 10) then
                  write(str_tmp,'(I1)') my_id
             else if(my_id .lt. 100) then
                  write(str_tmp,'(I2)') my_id
             else if(my_id .lt. 1000) then
                  write(str_tmp,'(I3)') my_id
             else if(my_id .lt. 10000) then
                  write(str_tmp,'(I4)') my_id
             else if(my_id .lt. 100000) then
                  write(str_tmp,'(I5)') my_id
             else
                continue
             endif
  open(iunit,file=trim(restart_filename_requested)//"."//trim(str_tmp),form="unformatted",ERR=101, access="sequential")


  read(iunit,ERR=101) TSLB  
  read(iunit,ERR=101) TSNOXY
  read(iunit,ERR=101) SMOIS     
  read(iunit,ERR=101) SH2O      
  read(iunit,ERR=101) ZSNSOXY   
  read(iunit,ERR=101) SNICEXY   
  read(iunit,ERR=101) SNLIQXY   
  read(iunit,ERR=101) QSNOWXY   
  read(iunit,ERR=101) FWETXY    
  read(iunit,ERR=101) SNEQVOXY  
  read(iunit,ERR=101) EAHXY     
  read(iunit,ERR=101) TAHXY     
  read(iunit,ERR=101) ALBOLDXY  
  read(iunit,ERR=101) CMXY      
  read(iunit,ERR=101) CHXY      
  read(iunit,ERR=101) ISNOWXY   
  read(iunit,ERR=101) CANLIQXY  
  read(iunit,ERR=101) CANICEXY  
  read(iunit,ERR=101) SNOW      
  read(iunit,ERR=101) SNOWH     
  read(iunit,ERR=101) TVXY      
  read(iunit,ERR=101) TGXY      
  read(iunit,ERR=101) ZWTXY     
  read(iunit,ERR=101) WAXY      
  read(iunit,ERR=101) WTXY      
  read(iunit,ERR=101) WSLAKEXY  
  read(iunit,ERR=101) LFMASSXY  
  read(iunit,ERR=101) RTMASSXY  
  read(iunit,ERR=101) STMASSXY  
  read(iunit,ERR=101) WOODXY    
  read(iunit,ERR=101) STBLCPXY  
  read(iunit,ERR=101) FASTCPXY  
  read(iunit,ERR=101) LAI       
  read(iunit,ERR=101) XSAIXY    
  read(iunit,ERR=101) VEGFRA    
  read(iunit,ERR=101) GVFMIN    
  read(iunit,ERR=101) GVFMAX    
  read(iunit,ERR=101) ACSNOM    
  read(iunit,ERR=101) ACSNOW    
  read(iunit,ERR=101) TAUSSXY   
  read(iunit,ERR=101) QSFC      
  read(iunit,ERR=101) SFCRUNOFF    
  read(iunit,ERR=101) UDRUNOFF     
!-----------------------------------
! Needed for water tracer capability
!----------------------------------
! not needed for now 
!----------------------------------

! #ifndef REALTIME
! #ifdef 1
!   read(iunit,ERR=101) ACCPRCP   
!   read(iunit,ERR=101) ACCECAN   
!   read(iunit,ERR=101) ACCEDIR   
!   read(iunit,ERR=101) ACCETRAN  
! #endif
! #endif
! ! below for opt_run = 5
!   if(IOPT_RUN .eq. 5) then
!       read(iunit,ERR=101) SMOISEQ   
!       read(iunit,ERR=101) AREAXY    
!       read(iunit,ERR=101) SMCWTDXY  
!       read(iunit,ERR=101) DEEPRECHXY
!       read(iunit,ERR=101) QSLATXY   
!       read(iunit,ERR=101) QRFSXY    
!       read(iunit,ERR=101) QSPRINGSXY
!       read(iunit,ERR=101) RECHXY    
!       read(iunit,ERR=101) QRFXY     
!       read(iunit,ERR=101) QSPRINGXY 
!       read(iunit,ERR=101) FDEPTHXY 
!       read(iunit,ERR=101) RIVERCONDXY 
!       read(iunit,ERR=101) RIVERBEDXY 
!       read(iunit,ERR=101) EQZWT 
!       read(iunit,ERR=101) PEXPXY 
!   endif
  close(iunit)

    endif
    call mpp_land_sync()
    i0 = i0 + istep
    ie = ie + istep
  end do ! end do of i loop

  return

101  continue
  call fatal_error_stop("FATAL ERROR: failed to read in lsm restartfile "   &
          //trim(restart_filename_requested)//"."//trim(str_tmp))
end subroutine lsm_rst_bi_in


   subroutine define_rst_variables(ncid)
      implicit none
      integer ncid
!      character(len=*) :: tmpStr

      call error_handler(ierr, "In module_hrldas_netcdf_io.F add_to_restart_2d_float() - "// &
                             "Problem nf90_open")
      ierr = nf90_redef(ncid)
      call error_handler(ierr, "In module_hrldas_netcdf_io.F add_to_restart_2d_float() - "// & 
                             "Problem nf90_redef")
print *, 'flag rst 01'
! add the variables

  call define_rst_var(ncid,TSLB      , "SOIL_T", layers="SOIL")
  call define_rst_var(ncid,TSNOXY    , "SNOW_T", layers="SNOW")
  call define_rst_var(ncid,SMOIS     , "SMC"   , layers="SOIL")
  call define_rst_var(ncid,SH2O      , "SH2O"  , layers="SOIL")
  call define_rst_var(ncid,ZSNSOXY   , "ZSNSO" , layers="SOSN")
  call define_rst_var(ncid,SNICEXY   , "SNICE" , layers="SNOW")
  call define_rst_var(ncid,SNLIQXY   , "SNLIQ" , layers="SNOW")
  call define_rst_var(ncid,QSNOWXY   , "QSNOW" )
  call define_rst_var(ncid,FWETXY    , "FWET"  )
  call define_rst_var(ncid,SNEQVOXY  , "SNEQVO")
  call define_rst_var(ncid,EAHXY     , "EAH"   )
  call define_rst_var(ncid,TAHXY     , "TAH"   )
  call define_rst_var(ncid,ALBOLDXY  , "ALBOLD")
  call define_rst_var(ncid,CMXY      , "CM"    )
  call define_rst_var(ncid,CHXY      , "CH"    )
  call define_rst_var(ncid,ISNOWXY   , "ISNOW" )
  call define_rst_var(ncid,CANLIQXY  , "CANLIQ")
  call define_rst_var(ncid,CANICEXY  , "CANICE")
  call define_rst_var(ncid,SNOW      , "SNEQV" )
  call define_rst_var(ncid,SNOWH     , "SNOWH" )
  call define_rst_var(ncid,TVXY      , "TV"    )
  call define_rst_var(ncid,TGXY      , "TG"    )
  call define_rst_var(ncid,ZWTXY     , "ZWT"   )
  call define_rst_var(ncid,WAXY      , "WA"    )
  call define_rst_var(ncid,WTXY      , "WT"    )
  call define_rst_var(ncid,WSLAKEXY  , "WSLAKE")
  call define_rst_var(ncid,LFMASSXY  , "LFMASS")
  call define_rst_var(ncid,RTMASSXY  , "RTMASS")
  call define_rst_var(ncid,STMASSXY  , "STMASS")
  call define_rst_var(ncid,WOODXY    , "WOOD"  )
  call define_rst_var(ncid,STBLCPXY  , "STBLCP")
  call define_rst_var(ncid,FASTCPXY  , "FASTCP")
  call define_rst_var(ncid,LAI       , "LAI"   )
  call define_rst_var(ncid,XSAIXY    , "SAI"   )
  call define_rst_var(ncid,VEGFRA    , "VEGFRA")
  call define_rst_var(ncid,GVFMIN    , "GVFMIN")
  call define_rst_var(ncid,GVFMAX    , "GVFMAX")
  call define_rst_var(ncid,ACSNOM    , "ACMELT")
  call define_rst_var(ncid,ACSNOW    , "ACSNOW")
  call define_rst_var(ncid,TAUSSXY   , "TAUSS" )
  call define_rst_var(ncid,QSFC      , "QSFC"  )
  call define_rst_var(ncid,SFCRUNOFF , "SFCRUNOFF")
  call define_rst_var(ncid,UDRUNOFF  , "UDRUNOFF" )
!-----------------------------------------------
! Neede for water tracer capaiblity
!----------------------------------------------
if (OPT_WT .eq. 1) then
  call define_rst_var(ncid,SMOIS_TR     , "SMC_TR"   , layers="TRAC")    !1
  call define_rst_var(ncid,SH2OXY_TR    , "SH2O_TR"  , layers="TRAC")    !2
  call define_rst_var(ncid,SNICEXY_TR   , "SNICE_TR" , layers="SNOW")    !3
  call define_rst_var(ncid,SNLIQXY_TR   , "SNLIQ_TR" , layers="SNOW")    !4
  call define_rst_var(ncid,QSNOWXY_TR   , "QSNOW_TR" )                   !5
  call define_rst_var(ncid,CANLIQXY_TR  , "CANLIQ_TR")                   !6
  call define_rst_var(ncid,CANICEXY_TR  , "CANICE_TR")                   !7
  call define_rst_var(ncid,SNOWXY_TR    , "SNEQV_TR" )                   !8
  call define_rst_var(ncid,SNOWHXY_TR   , "SNOWH_TR" )                   !9
  call define_rst_var(ncid,WAXY_TR      , "WA_TR"    )                   !10
  call define_rst_var(ncid,WTXY_TR      , "WT_TR"    )                   !11
  call define_rst_var(ncid,ACSNOM_TR    , "ACMELT_TR")                   !12
  call define_rst_var(ncid,ACSNOW_TR    , "ACSNOW_TR")                   !13
  call define_rst_var(ncid,SFCRUNOFF_TR , "SFCRUNOFF_TR")                !14
  call define_rst_var(ncid,UDRUNOFF_TR  , "UDRUNOFF_TR" )                !15
  call define_rst_var(ncid,SMOIS_SUB    , "SMC_SUB"   , layers="TRAC")   !16
  call define_rst_var(ncid,SH2OXY_SUB   , "SH2O_SUB"  , layers="TRAC")   !17
end if
!-----------------------------------------------
  call define_rst_var(ncid,ACCPRCP   , "ACCPRCP" )
  call define_rst_var(ncid,ACCECAN   , "ACCECAN" )
  call define_rst_var(ncid,ACCEDIR   , "ACCEDIR" )
  call define_rst_var(ncid,ACCETRAN  , "ACCETRAN" )
! below for opt_run = 5
  call define_rst_var(ncid,SMOISEQ   , "SMOISEQ"  , layers="SOIL"  )
  call define_rst_var(ncid,AREAXY    , "AREAXY"     )
  call define_rst_var(ncid,SMCWTDXY  , "SMCWTDXY"   )
  call define_rst_var(ncid,DEEPRECHXY, "DEEPRECHXY" )
  call define_rst_var(ncid,QSLATXY   , "QSLATXY"    )
  call define_rst_var(ncid,QRFSXY    , "QRFSXY"     )
  call define_rst_var(ncid,QSPRINGSXY, "QSPRINGSXY" )
  call define_rst_var(ncid,RECHXY    , "RECHXY"     )
  call define_rst_var(ncid,QRFXY     , "QRFXY"      )
  call define_rst_var(ncid,QSPRINGXY , "QSPRINGXY"  )
  call define_rst_var(ncid,FDEPTHXY , "FDEPTHXY"  )
  call define_rst_var(ncid,RIVERCONDXY , "RIVERCONDXY"  )
  call define_rst_var(ncid,RIVERBEDXY , "RIVERBEDXY"  )
  call define_rst_var(ncid,EQZWT , "EQZWT"  )
  call define_rst_var(ncid,PEXPXY , "PEXPXY"  )
      ierr = nf90_enddef(ncid)
      call error_handler(ierr, "In module_hrldas_netcdf_io.F add_to_restart_2d_float() - "// &
                             "Problem nf90_enddef")
      call error_handler(ierr, "In module_hrldas_netcdf_io.F add_to_restart_3d() - "// &
                             "Problem nf90_close")

   end subroutine define_rst_variables

end module module_NoahMP_hrldas_driver

!subroutine wrf_message(msg)
!  implicit none
!  character(len=*), intent(in) :: msg
!  print*, msg
!end subroutine wrf_message

logical function wrf_dm_on_monitor() result(l)
  l = .TRUE.
  return
end function wrf_dm_on_monitor


!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

SUBROUTINE CALC_DECLIN ( NOWDATE, LATITUDE, LONGITUDE, COSZ, JULIAN )

  USE MODULE_DATE_UTILITIES
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

   REAL, PARAMETER :: DEGRAD = 3.14159265/180.
   REAL, PARAMETER :: DPD    = 360./365.
! !ARGUMENTS:
   CHARACTER(LEN=19), INTENT(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
   REAL,              INTENT(IN)  :: LATITUDE
   REAL,              INTENT(IN)  :: LONGITUDE
   REAL,              INTENT(OUT) :: COSZ
   REAL,              INTENT(OUT) :: JULIAN
   REAL                           :: HRANG
   REAL                           :: DECLIN
   REAL                           :: OBECL
   REAL                           :: SINOB
   REAL                           :: SXLONG
   REAL                           :: ARG
   REAL                           :: TLOCTIM
   INTEGER                        :: IDAY
   INTEGER                        :: IHOUR
   INTEGER                        :: IMINUTE
   INTEGER                        :: ISECOND

   CALL GETH_IDTS(NOWDATE(1:10), NOWDATE(1:4)//"-01-01", IDAY)
   READ(NOWDATE(12:13), *) IHOUR
   READ(NOWDATE(15:16), *) IMINUTE
   READ(NOWDATE(18:19), *) ISECOND
   JULIAN = REAL(IDAY) + REAL(IHOUR)/24.

!
! FOR SHORT WAVE RADIATION

   DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! LOCAL TIME IN HOURS
   TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)
   HRANG=15.*(TLOCTIM-12.)*DEGRAD
   COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)

 END SUBROUTINE CALC_DECLIN

!
!------------------------------------------------------------------------------------------
