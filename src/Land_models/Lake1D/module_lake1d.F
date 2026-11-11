!===============================================================================
! module_lake1d: driver for the 1-D lake model (module_sf_lake) in WRF-Hydro.
!
! Holds the lake state and runs the lake model on ISLAKE cells, coupled to the
! WRF-Hydro level-pool reservoirs (water_temperature branch): river inflow heat,
! outlet-depth (withdrawal) temperature for the stream temperature model, lake
! precipitation/evaporation for the reservoir water balance, and the lake depth
! following the reservoir level.  Used by both drivers:
!   - the offline HRLDAS/Noah-MP driver (module_NoahMP_hrldas_driver), and
!   - the NUOPC cap (coupled to an external LSM, e.g. LIS in the NASA Land Coupler).
! The code was moved here unchanged from module_NoahMP_hrldas_driver.
!===============================================================================
module module_lake1d
  use module_sf_lake, only: lake, lakeini, lakeini_rst, nlevsoil, nlevlake, nlevsnow, &
                            lake_mixfact2d, lake_kbg2d
  use config_base,    only: noah_lsm, nlst
#ifdef WRF_HYDRO
  use module_rt_data, only: rt_domain
#endif
#ifdef MPP_LAND
  use module_mpp_land, only: sum_real1d, sum_int1d, write_IO_real, my_id, IO_id, &
                             global_nx, global_ny, mpp_land_sync
  use module_cpl_land, only: fatal_error_stop
#endif
  use netcdf
  implicit none

  public
!------------------------------------------------------------------------
! Needed for lake model (sf_lake_physics = 1)
!------------------------------------------------------------------------

  INTEGER                             :: lakeflag
  INTEGER                             :: lake_depth_flag
  LOGICAL, ALLOCATABLE, DIMENSION(:,:):: lake2d

  REAL, ALLOCATABLE, DIMENSION(:,:)   :: lake_depth
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: lakedepth2d
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: savedtke12d
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: snowdp2d
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: h2osno2d
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: snl2d
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: t_grnd2d
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: lakemask
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: t2        ! temporary add these since they dont exist
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: th2       ! temporary add these since they dont exist
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: q2        ! temporary add these since they dont exist

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: t_lake3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: lake_icefrac3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z_lake3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dz_lake3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: t_soisno3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: h2osoi_ice3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: h2osoi_liq3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: h2osoi_vol3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dz3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: zi3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: watsat3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: csol3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tkmg3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tkdry3d
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tksatu3d

!------------------------------------------------------------------------
! Needed for wrf_hydro lake model
!------------------------------------------------------------------------
#ifdef WRF_HYDRO
  INTEGER                                :: nlakes

  REAL   , ALLOCATABLE, DIMENSION(:)     :: lake_precip
  REAL   , ALLOCATABLE, DIMENSION(:)     :: lake_evap
  REAL   , ALLOCATABLE, DIMENSION(:)     :: lake_netoutflow
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: lake_grids
  REAL   , ALLOCATABLE, DIMENSION(:)     :: lake_dh_avg      ! change in lake water level this LSM step, from the level-pool reservoir (m)
  REAL   , ALLOCATABLE, DIMENSION(:)     :: lake_resht_prev  ! level-pool water elevation at the previous LSM step (m)
  LOGICAL                                :: lake_resht_init = .false.
  REAL   , PARAMETER                     :: lake_min_depth = 1.0  ! floor on 1-D lake model depth (m)
  REAL   , ALLOCATABLE, DIMENSION(:)     :: lake_t_sum       ! per-lake sum of outlet-depth temperature over this rank's cells (K)
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: lake_t_ncell     ! per-lake count of this rank's lake cells
  REAL   , ALLOCATABLE, DIMENSION(:)     :: lake_t_sum_dp    ! as lake_t_sum, cells whose column reaches the outlet depth
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: lake_t_ncell_dp  ! as lake_t_ncell, cells whose column reaches the outlet depth
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: lake_ix          ! routing lake index of each lake cell (0 = none)
  ! 1-D lake model sub-stepping: the lake model (with its multi-layer snow)
  ! is run with a time step no longer than lake_dt_max; fluxes are averaged
  ! over the sub-steps.  Rain-on-snow at a 3600 s step produced runaway
  ! lake snow-layer temperatures (t_soisno > 500 K).
  REAL   , PARAMETER                     :: lake_dt_max = 600.0  ! maximum lake sub-step (s)
  REAL   , ALLOCATABLE, DIMENSION(:,:)   :: rain_lk, hfx_lk, lh_lk, grdflx_lk, qfx_lk
  LOGICAL                                :: lake_heat_warned = .false.

  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: lake_id
#endif

#ifdef WRF_HYDRO
  ! --- Coupled use (NUOPC cap, e.g. LIS as the land model): land-surface fields
  ! from the setup (wrfinput) file, and the lake-model inputs/outputs that the
  ! offline Noah-MP driver otherwise provides (lake1d_cpl_init/step/post).
  integer :: cpl_iswater, cpl_islake
  integer, allocatable, dimension(:,:)   :: cpl_ivgtyp, cpl_isltyp
  real,    allocatable, dimension(:,:)   :: cpl_terrain, cpl_xland, cpl_xice, cpl_tsk, cpl_snow, cpl_xlat, &
                                            cpl_albedo, cpl_emiss, cpl_hfx, cpl_lh, cpl_grdflx, cpl_qfx,    &
                                            cpl_glw, cpl_swdown, cpl_rainbl
  real,    allocatable, dimension(:,:,:) :: cpl_t_phy, cpl_p8w, cpl_dz8w, cpl_qv, cpl_u, cpl_v
  ! Forcing over lakes: the land model provides no forcing on lake cells, so each
  ! lake cell gets the mean over the land cells on that lake's shoreline.
  integer, allocatable, dimension(:,:)   :: cpl_shore_ix    ! routing lake index of a shoreline land cell (0 = none)
  integer, allocatable, dimension(:,:)   :: cpl_lake_ix     ! routing lake index of a lake (LAKE_ID) cell (0 = none)
  integer :: cpl_nlkgrp = 0        ! number of lake groups (distinct LAKE_IDs) for the lake-cell forcing
  integer :: cpl_out_seconds = 0   ! model seconds since the start of the run (LAKE1D output timing)
#endif

  ! Lake state variables that must all be present in a restart file to restart the lake model
  character(len=16), dimension(22), parameter :: lake_rst_vars = [character(len=16) :: &
       'LAKE_DEPTH', 'LAKE_SAVETKE1', 'LAKE_SNOWP', 'LAKE_H2OSNO', 'LAKE_SNL', 'LAKE_T_GRND', &
       'LAKE_T_SOILSNOW', 'LAKE_H2OSOI_ICE', 'LAKE_H2OSOI_LIQ', 'LAKE_H2OSOI_VOL', 'LAKE_Z', 'LAKE_DZ', &
       'LAKE_T_LAKE', 'LAKE_ICEFRAC', 'LAKE_Z_LAKE', 'LAKE_DZ_LAKE', &
       'LAKE_WATSAT', 'LAKE_CSOL', 'LAKE_TKMG', 'LAKE_TKDRY', 'LAKE_TKSATU', 'LAKE_ZI']
  ! Lake surface fields in LAKE_RST files (boundary conditions for a coupled
  ! atmosphere over lake cells); optional on input
  character(len=16), dimension(10), parameter :: lake_rst_sfc_vars = [character(len=16) :: &
       'LAKE_TSK', 'LAKE_HFX', 'LAKE_LH', 'LAKE_QFX', 'LAKE_GRDFLX', 'LAKE_ALBEDO', 'LAKE_EMISS', &
       'LAKE_T2', 'LAKE_TH2', 'LAKE_Q2']

contains

!-------------------------------------------------------------------------------
! Allocate the lake state on the tile xstart:xend, ystart:yend (global indices).
! The lake-model arrays are allocated only with the lake model on
! (sf_lake_physics = 1); lake_id is always allocated.
!-------------------------------------------------------------------------------
  subroutine lake1d_allocate(XSTART, XEND, YSTART, YEND)
    implicit none
    integer, intent(in) :: XSTART, XEND, YSTART, YEND

  IF (noah_lsm%sf_lake_physics == 1) THEN

  ALLOCATE ( lake2d(XSTART:XEND,YSTART:YEND))

  ALLOCATE (  lake_depth(XSTART:XEND,YSTART:YEND))
  ALLOCATE ( lakedepth2d(XSTART:XEND,YSTART:YEND))
  ALLOCATE ( savedtke12d(XSTART:XEND,YSTART:YEND))
  ALLOCATE (    snowdp2d(XSTART:XEND,YSTART:YEND))
  ALLOCATE (    h2osno2d(XSTART:XEND,YSTART:YEND))
  ALLOCATE (       snl2d(XSTART:XEND,YSTART:YEND))
  ALLOCATE (    t_grnd2d(XSTART:XEND,YSTART:YEND))
  ALLOCATE (    lakemask(XSTART:XEND,YSTART:YEND))
  ALLOCATE (          t2(XSTART:XEND,YSTART:YEND)) ! temporary add these since they dont exist
  ALLOCATE (         th2(XSTART:XEND,YSTART:YEND)) ! temporary add these since they dont exist
  ALLOCATE (          q2(XSTART:XEND,YSTART:YEND)) ! temporary add these since they dont exist

  ALLOCATE (       t_lake3d(XSTART:XEND,          1:nlevlake,YSTART:YEND)) 
  ALLOCATE ( lake_icefrac3d(XSTART:XEND,          1:nlevlake,YSTART:YEND)) 
  ALLOCATE (       z_lake3d(XSTART:XEND,          1:nlevlake,YSTART:YEND)) 
  ALLOCATE (      dz_lake3d(XSTART:XEND,          1:nlevlake,YSTART:YEND)) 
  ALLOCATE (     t_soisno3d(XSTART:XEND,-nlevsnow+1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (   h2osoi_ice3d(XSTART:XEND,-nlevsnow+1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (   h2osoi_liq3d(XSTART:XEND,-nlevsnow+1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (   h2osoi_vol3d(XSTART:XEND,-nlevsnow+1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (            z3d(XSTART:XEND,-nlevsnow+1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (           dz3d(XSTART:XEND,-nlevsnow+1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (           zi3d(XSTART:XEND,-nlevsnow+0:nlevsoil,YSTART:YEND)) 
  ALLOCATE (       watsat3d(XSTART:XEND,          1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (         csol3d(XSTART:XEND,          1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (         tkmg3d(XSTART:XEND,          1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (        tkdry3d(XSTART:XEND,          1:nlevsoil,YSTART:YEND)) 
  ALLOCATE (       tksatu3d(XSTART:XEND,          1:nlevsoil,YSTART:YEND)) 

  ! Defined values on all cells: lakeini / lakeini_rst set the state on lake
  ! cells only, and restart files would otherwise hold uninitialized memory
  ! on the other cells.
  lake2d = .false.
  lake_depth = 0.0;  lakedepth2d = 0.0;  savedtke12d = 0.0;  snowdp2d = 0.0
  h2osno2d = 0.0;    snl2d = 0.0;        t_grnd2d = 0.0;     lakemask = 0.0
  t2 = 0.0;          th2 = 0.0;          q2 = 0.0
  t_lake3d = 0.0;    lake_icefrac3d = 0.0;  z_lake3d = 0.0;  dz_lake3d = 0.0
  t_soisno3d = 0.0;  h2osoi_ice3d = 0.0;    h2osoi_liq3d = 0.0;  h2osoi_vol3d = 0.0
  z3d = 0.0;         dz3d = 0.0;            zi3d = 0.0
  watsat3d = 0.0;    csol3d = 0.0;  tkmg3d = 0.0;  tkdry3d = 0.0;  tksatu3d = 0.0

  ENDIF  ! sf_lake_physics == 1

  ALLOCATE ( lake_id(XSTART:XEND,YSTART:YEND))

  end subroutine lake1d_allocate

!-------------------------------------------------------------------------------
! Lake mask, lake depth and LAKE_ID from the setup (wrfinput) file.
! Call with the lake model on (sf_lake_physics = 1).
!-------------------------------------------------------------------------------
  subroutine lake1d_setup(XSTART, XEND, YSTART, YEND, ivgtyp, iswater, islake)
    implicit none
    integer, intent(in) :: XSTART, XEND, YSTART, YEND
    integer, intent(in) :: ivgtyp(XSTART:XEND,YSTART:YEND)
    integer, intent(in) :: iswater, islake
    integer            :: lake_id_bad(1)               ! LAKE_ID cells whose land use is not ISLAKE

     lakemask = 0
     where(ivgtyp ==  islake) lakemask = 1

     lakeflag = 1
     if(iswater < 0) lakeflag = 0

     lake_depth_flag = 0
     if(noah_lsm%use_lakedepth == 1) then
       call read_lake_depth(noah_lsm%HRLDAS_SETUP_FILE, XSTART, XEND, YSTART, YEND,LAKE_DEPTH)
       lake_depth_flag = 1
     end if

     lake_id = 0
     call read_lake_id(noah_lsm%hrldas_setup_file, xstart, xend, ystart, yend, lake_id)

     ! The 1-D lake model runs only on ISLAKE cells, but lake outflow
     ! temperature, precipitation and evaporation are averaged over all cells
     ! of a lake's LAKE_ID footprint.  A LAKE_ID cell with another land use
     ! would add unset lake temperatures and land P/E to the lake means, so
     ! LAKE_ID and IVGTYP must be consistent (a pre-processing requirement).
     lake_id_bad(1) = count(lake_id > 0 .and. ivgtyp /= islake)
#ifdef MPP_LAND
     call sum_int1d(lake_id_bad, 1)
#endif
     if (lake_id_bad(1) > 0) then
        write(6,'(a,i0,a,i0,a)') 'FATAL ERROR: ', lake_id_bad(1), &
             ' cells in wrfinput have LAKE_ID > 0 but IVGTYP /= ISLAKE (', islake, &
             '). LAKE_ID lake footprints must have the lake land-use class; fix the domain pre-processing.'
        call flush(6)
        call fatal_error_stop('FATAL ERROR: LAKE_ID cells without lake land use (IVGTYP /= ISLAKE) in wrfinput')
     end if

  end subroutine lake1d_setup

!-------------------------------------------------------------------------------
! Initialize the lake state: lakeini (cold start) or lakeini_rst (state read
! from a restart file into the module arrays).  ivgtyp, xland, xice may be
! changed on lake cells (as lakeini/lakeini_rst do).
!-------------------------------------------------------------------------------
  subroutine lake1d_init(from_restart, ivgtyp, isltyp, terrain, snow, restart_flag,        &
                         iswater, xice, xice_threshold, xland, tsk,                       &
                         its, ite, jts, jte, ims, ime, jms, jme)
    implicit none
    logical, intent(in)    :: from_restart, restart_flag
    integer, intent(in)    :: its, ite, jts, jte, ims, ime, jms, jme, iswater
    integer, intent(inout) :: ivgtyp(ims:ime,jms:jme)
    integer, intent(inout) :: isltyp(ims:ime,jms:jme)
    real,    intent(in)    :: terrain(ims:ime,jms:jme), tsk(ims:ime,jms:jme), xice_threshold
    real,    intent(inout) :: snow(ims:ime,jms:jme), xice(ims:ime,jms:jme), xland(ims:ime,jms:jme)

    if (from_restart) then
       call  lakeini_rst(ivgtyp,          isltyp,        terrain,           snow,           & !i
                         noah_lsm%lake_min_elev,   restart_flag,        noah_lsm%lakedepth_default, lake_depth,     &
                         lakedepth2d,     savedtke12d,   snowdp2d,          h2osno2d,       & !o
                         snl2d,           t_grnd2d,      t_lake3d,          lake_icefrac3d, &
                         z_lake3d,        dz_lake3d,     t_soisno3d,        h2osoi_ice3d,   &
                         h2osoi_liq3d,    h2osoi_vol3d,  z3d,               dz3d,           &
                         zi3d,            watsat3d,      csol3d,            tkmg3d,         &
                         iswater,         xice,          xice_threshold,     xland,  tsk,   &
                         lakemask,        lakeflag,                                         &
                         lake_depth_flag, noah_lsm%use_lakedepth,                           &
                         tkdry3d,         tksatu3d,      lake2d,                            &
                         its,ite,jts,jte,ims,ime,jms,jme)
    else
       call  lakeini(ivgtyp,          isltyp,        terrain,           snow,           & !i
                     noah_lsm%lake_min_elev,   restart_flag,       noah_lsm%lakedepth_default, lake_depth,     &
                     lakedepth2d,     savedtke12d,   snowdp2d,          h2osno2d,       & !o
                     snl2d,           t_grnd2d,      t_lake3d,          lake_icefrac3d, &
                     z_lake3d,        dz_lake3d,     t_soisno3d,        h2osoi_ice3d,   &
                     h2osoi_liq3d,    h2osoi_vol3d,  z3d,               dz3d,           &
                     zi3d,            watsat3d,      csol3d,            tkmg3d,         &
                     iswater,         xice,          xice_threshold,     xland,  tsk,   &
                     lakemask,        lakeflag,                                         &
                     lake_depth_flag, noah_lsm%use_lakedepth,                           &
                     tkdry3d,         tksatu3d,      lake2d,                            &
                     its,ite,jts,jte,ims,ime,jms,jme)
    end if

  end subroutine lake1d_init

#ifdef WRF_HYDRO
!-------------------------------------------------------------------------------
! Per-lake arrays exchanged with the level-pool reservoirs (nlakes from routing).
!-------------------------------------------------------------------------------
  subroutine lake1d_allocate_lakes(n)
    implicit none
    integer, intent(in) :: n
    nlakes = n
   allocate(lake_precip(nlakes))
   allocate(lake_evap(nlakes))
   allocate(lake_netoutflow(nlakes))
   allocate(lake_grids(nlakes))
   allocate(lake_dh_avg(nlakes))
   allocate(lake_resht_prev(nlakes))
   allocate(lake_t_sum(nlakes))
   allocate(lake_t_ncell(nlakes))
   allocate(lake_t_sum_dp(nlakes))
   allocate(lake_t_ncell_dp(nlakes))
   lake_precip = 0.0
   lake_evap   = 0.0
   lake_grids  = 0
   lake_dh_avg = 0.0
   lake_resht_prev = 0.0
   lake_netoutflow = 0.0
   lake_t_sum   = 0.0
   lake_t_ncell = 0
   lake_t_sum_dp   = 0.0
   lake_t_ncell_dp = 0
  end subroutine lake1d_allocate_lakes
#endif

!-------------------------------------------------------------------------------
! One LSM time step of the lake model on the tile, with its coupling to the
! level-pool reservoirs: sub-stepped 1-D lake physics, river inflow heat,
! whole-lake precipitation/evaporation and outlet-depth temperature (summed
! across ranks), lake depth following the reservoir level, LakeOutflowTemp.
! Forcing as in the HRLDAS driver: 3-D fields at level 1 (p8w also level 2),
! qv_curr mixing ratio [kg/kg], rainbl precipitation [mm per step].
! hfx, lh, grdflx, tsk, qfx are set on lake cells.
!-------------------------------------------------------------------------------
  subroutine lake1d_step(XSTART, XEND, YSTART, YEND, dtbl,                               &
                         t_phy, p8w, dz8w, qv_curr, u_phy, v_phy, glw, emiss, rainbl,     &
                         swdown, albedo, xlat_urb2d, ivgtyp, terrain, xland, iswater,     &
                         xice, xice_threshold,                                            &
                         ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,      &
                         its, ite, jts, jte, kts, kte,                                    &
                         hfx, lh, grdflx, tsk, qfx)
    implicit none
    integer, intent(in) :: XSTART, XEND, YSTART, YEND, iswater
    integer, intent(in) :: ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte
    real,    intent(in) :: dtbl, xice_threshold
    real,    dimension(ims:ime,kms:kme,jms:jme), intent(in) :: t_phy, p8w, dz8w, qv_curr, u_phy, v_phy
    real,    dimension(ims:ime,jms:jme), intent(in)    :: glw, emiss, rainbl, swdown, xlat_urb2d, terrain, xland
    real,    dimension(ims:ime,jms:jme), intent(inout) :: albedo, xice
    integer, dimension(ims:ime,jms:jme), intent(inout) :: ivgtyp
    real,    dimension(ims:ime,jms:jme), intent(inout) :: hfx, lh, grdflx, tsk, qfx

    integer :: isub_lk, nsub_lk   ! 1-D lake model sub-steps
    real    :: dt_lk              ! 1-D lake model sub-step (s)
#ifdef WRF_HYDRO
    integer :: ilake
    integer :: lake_idx
    integer, parameter :: did=1
    ! lake outlet-depth temperature sampling
    integer :: il, jl, kl, ncell_lk, ncell_dp
    real    :: z_out, t_sum, t_col, wgt_lk, t_sum_dp
#endif


         ! Sub-step the 1-D lake model (see lake_dt_max).  Each sub-step gets an
         ! equal share of the step's precipitation (rainbl is an amount per
         ! step); fluxes are averaged over the sub-steps and the surface state
         ! (TSK, T2, Q2, albedo) is the last sub-step's.  With one sub-step this
         ! reproduces the single call exactly.
         nsub_lk = max(1, ceiling(dtbl / lake_dt_max - 1.0e-4))
         dt_lk   = dtbl / real(nsub_lk)
         if (.not. allocated(rain_lk)) then
            allocate(rain_lk(XSTART:XEND,YSTART:YEND), hfx_lk(XSTART:XEND,YSTART:YEND),    &
                     lh_lk(XSTART:XEND,YSTART:YEND),   grdflx_lk(XSTART:XEND,YSTART:YEND), &
                     qfx_lk(XSTART:XEND,YSTART:YEND))
         end if
         rain_lk   = rainbl / real(nsub_lk)
         hfx_lk    = 0.0
         lh_lk     = 0.0
         grdflx_lk = 0.0
         qfx_lk    = 0.0

#ifdef WRF_HYDRO
         ! Per-cell mixing parameters for the 1-D lake model from the per-lake
         ! LAKEPARM LakeMixFact and LakeKbg (filled once; cells of lakes not in
         ! LAKEPARM keep 0, i.e. the lake model's built-in default).
         if (.not. allocated(lake_mixfact2d) .and. nlakes .gt. 0) then
            if (allocated(rt_domain(did)%LakeMixFact)) then
               allocate(lake_mixfact2d(XSTART:XEND,YSTART:YEND))
               lake_mixfact2d = 0.0
               do ilake = 1, nlakes
                  where (lake_id == rt_domain(did)%lakeidm(ilake)) &
                     lake_mixfact2d = rt_domain(did)%LakeMixFact(ilake)
               end do
            end if
         end if
         if (.not. allocated(lake_kbg2d) .and. nlakes .gt. 0) then
            if (allocated(rt_domain(did)%LakeKbg)) then
               allocate(lake_kbg2d(XSTART:XEND,YSTART:YEND))
               lake_kbg2d = 0.0
               do ilake = 1, nlakes
                  where (lake_id == rt_domain(did)%lakeidm(ilake)) &
                     lake_kbg2d = rt_domain(did)%LakeKbg(ilake)
               end do
            end if
         end if
#endif

         do isub_lk = 1, nsub_lk
         CALL Lake(  t_phy        ,p8w            ,dz8w         ,qv_curr         ,&  !i
                     u_phy        ,v_phy          ,glw          ,emiss           ,&
                     rain_lk      ,dt_lk          ,swdown       ,albedo          ,&
                     xlat_urb2d   ,z_lake3d       ,dz_lake3d    ,lakedepth2d     ,&
                     watsat3d     ,csol3d         ,tkmg3d       ,tkdry3d         ,&
                     tksatu3d     ,ivgtyp         ,terrain      ,xland           ,&
                     iswater      ,xice           ,xice_threshold, noah_lsm%lake_min_elev ,&
                     ids          ,ide            ,jds          ,jde             ,&
                     kds          ,kde            ,ims          ,ime             ,&
                     jms          ,jme            ,kms          ,kme             ,&
                     its          ,ite            ,jts          ,jte             ,&
                     kts          ,kte                                           ,&
                     h2osno2d     ,snowdp2d       ,snl2d        ,z3d             ,&  !h
                     dz3d         ,zi3d           ,h2osoi_vol3d ,h2osoi_liq3d    ,&
                     h2osoi_ice3d ,t_grnd2d       ,t_soisno3d   ,t_lake3d        ,&
                     savedtke12d  ,lake_icefrac3d                                ,&
                     lakemask                                                    ,&
                     hfx          ,lh             ,grdflx       ,tsk             ,&  !o
                     qfx          ,t2             ,th2          ,q2 )
            where (lakemask == 1)
               hfx_lk    = hfx_lk    + hfx
               lh_lk     = lh_lk     + lh
               grdflx_lk = grdflx_lk + grdflx
               qfx_lk    = qfx_lk    + qfx
            end where
         end do
         where (lakemask == 1)
            hfx    = hfx_lk    / real(nsub_lk)
            lh     = lh_lk     / real(nsub_lk)
            grdflx = grdflx_lk / real(nsub_lk)
            qfx    = qfx_lk    / real(nsub_lk)
         end where

#ifdef WRF_HYDRO
         ! Heat carried into the lakes by river inflow (water_temperature = 2).
         if (nlst(did)%water_temperature .eq. 2 .and. nlakes .gt. 0) then
            if (allocated(rt_domain(did)%LakeInflowTemp) .and. allocated(rt_domain(did)%QLAKEI) &
                .and. allocated(rt_domain(did)%LakeArea)) then
               if (.not. allocated(lake_ix)) then
                  allocate(lake_ix(XSTART:XEND,YSTART:YEND))
                  lake_ix = 0
                  do ilake = 1, nlakes
                     where (lake_id == rt_domain(did)%lakeidm(ilake)) lake_ix = ilake
                  end do
               end if
               call lake_inflow_heat(XSTART, XEND, YSTART, YEND, nlakes, lake_ix, lakemask,   &
                    dtbl, rt_domain(did)%QLAKEI(1:nlakes), rt_domain(did)%LakeInflowTemp(1:nlakes), &
                    rt_domain(did)%LakeArea(1:nlakes), rt_domain(did)%LakeOutletDepth(1:nlakes),  &
                    z_lake3d, dz_lake3d, lake_icefrac3d, t_lake3d)
            else if (.not. lake_heat_warned) then
               write(6,*) 'WARNING: lake inflow heat disabled: lake inflow, inflow temperature ', &
                          'or lake area arrays are not allocated.'
               lake_heat_warned = .true.
            end if
         end if
#endif

         ! fake values for testing:
         !where(ivgtyp == islake) rainbl  = 2.0
         !where(ivgtyp == islake) qfx = 1.0
#ifdef WRF_HYDRO
#endif
         do ilake = 1, nlakes
           
           ! find ID for this lake.  These are sums over this rank's tile only;
           ! they are reduced across ranks after the loop.
           lake_idx = rt_domain(did)%lakeidm(ilake)
           lake_grids(ilake) = count(mask = lake_id == lake_idx)
           lake_precip(ilake) = sum(rainbl, mask = lake_id == lake_idx)
           lake_evap(ilake)   = sum(qfx, mask = lake_id == lake_idx)

           ! Pass water temperature from 1D lake physics to routing [K], sampled
           ! at the lake's outlet (withdrawal) depth below the surface.  Each
           ! column's profile t_lake3d/z_lake3d is interpolated at the outlet
           ! depth (clamped to the column), then averaged over the lake's cells
           ! whose column (bottom of the lowest lake layer) reaches the outlet
           ! depth, so shallow coves and arms do not contribute their warm
           ! bottom water to a deep withdrawal.  If no cell reaches the outlet
           ! depth, all cells are averaged (each at its clamped depth).
           ! LakeOutletDepth = 0 reproduces the legacy surface-layer withdrawal.
           z_out    = rt_domain(did)%LakeOutletDepth(ilake)
           t_sum    = 0.0
           ncell_lk = 0
           t_sum_dp = 0.0
           ncell_dp = 0
           do jl = YSTART, YEND
              do il = XSTART, XEND
                 if (lake_id(il,jl) == lake_idx) then
                    if (z_out .le. z_lake3d(il,1,jl)) then
                       t_col = t_lake3d(il,1,jl)
                    else if (z_out .ge. z_lake3d(il,nlevlake,jl)) then
                       t_col = t_lake3d(il,nlevlake,jl)
                    else
                       t_col = t_lake3d(il,1,jl)
                       do kl = 1, nlevlake-1
                          if (z_out .le. z_lake3d(il,kl+1,jl)) then
                             wgt_lk = (z_out - z_lake3d(il,kl,jl)) / &
                                      max(1.0e-6, z_lake3d(il,kl+1,jl) - z_lake3d(il,kl,jl))
                             t_col  = (1.0 - wgt_lk) * t_lake3d(il,kl,jl) &
                                    +        wgt_lk  * t_lake3d(il,kl+1,jl)
                             exit
                          end if
                       end do
                    end if
                    t_sum    = t_sum + t_col
                    ncell_lk = ncell_lk + 1
                    if (z_out .le. z_lake3d(il,nlevlake,jl) + 0.5 * dz_lake3d(il,nlevlake,jl)) then
                       t_sum_dp = t_sum_dp + t_col
                       ncell_dp = ncell_dp + 1
                    end if
                 end if
              end do
           end do
           lake_t_sum(ilake)      = t_sum
           lake_t_ncell(ilake)    = ncell_lk
           lake_t_sum_dp(ilake)   = t_sum_dp
           lake_t_ncell_dp(ilake) = ncell_dp

         end do

         ! Reduce the outlet-temperature sums across ranks: a lake's cells live
         ! on the land tiles of a few ranks, but every rank needs the resulting
         ! LakeOutflowTemp (the rank owning the outlet reach in the channel
         ! decomposition generally has no lake cells and would otherwise keep
         ! the 10 K init value, tripping the stream model's fallback).
#ifdef MPP_LAND
         call sum_real1d(lake_t_sum, nlakes)
         call sum_int1d(lake_t_ncell, nlakes)
         call sum_real1d(lake_t_sum_dp, nlakes)
         call sum_int1d(lake_t_ncell_dp, nlakes)
#endif

         ! Reduce lake cell counts and P/E sums across ranks so every rank holds
         ! whole-lake totals.  The coupler (hrldas_cpl_HYDRO) divides these by
         ! lake_grids to get lake-mean P and E for the level-pool reservoirs;
         ! the reservoir arrays must be identical on all ranks because
         ! updateLake_seq adds each rank's deviation from the IO rank's values.
#ifdef MPP_LAND
         call sum_int1d(lake_grids, nlakes)
         call sum_real1d(lake_precip, nlakes)
         call sum_real1d(lake_evap, nlakes)
#endif

         ! Update the 1-D lake model depth to follow the level-pool reservoir.
         ! The level pool is the single lake water balance (inflow - outflow
         ! + P - E, with P and E from this lake model passed via the coupler);
         ! the lake model's depth is shifted by the reservoir's water-elevation
         ! change since the previous LSM step.  RESHT is synchronized across
         ! ranks after routing, so every rank applies the same change.  The
         ! reference elevation (previous step's RESHT) is kept in
         ! rt_domain%LakeReshtLSM, which is written to and read from HYDRO_RST,
         ! so a restart applies the same depth change as a continuous run.  On a
         ! cold start, or a restart file without it, it is initialized on the
         ! first call with no depth change.
         if (allocated(rt_domain(did)%RESHT)) then
            if (size(rt_domain(did)%RESHT) .ge. nlakes) then
               if (.not. lake_resht_init) then
                  lake_resht_prev(1:nlakes) = rt_domain(did)%RESHT(1:nlakes)
                  if (allocated(rt_domain(did)%LakeReshtLSM)) then
                     if (all(rt_domain(did)%LakeReshtLSM(1:nlakes) .gt. -9000.0)) &
                        lake_resht_prev(1:nlakes) = rt_domain(did)%LakeReshtLSM(1:nlakes)
                  end if
                  lake_resht_init = .true.
               end if
               do ilake = 1, nlakes
                  lake_dh_avg(ilake) = rt_domain(did)%RESHT(ilake) - lake_resht_prev(ilake)
                  lake_resht_prev(ilake) = rt_domain(did)%RESHT(ilake)
                  if (lake_grids(ilake) .gt. 0) then
                     lake_idx = rt_domain(did)%lakeidm(ilake)
                     where(lake_id == lake_idx) &
                        lakedepth2d = max(lake_min_depth, lakedepth2d + lake_dh_avg(ilake))
                  end if
               end do
               if (allocated(rt_domain(did)%LakeReshtLSM)) &
                  rt_domain(did)%LakeReshtLSM(1:nlakes) = lake_resht_prev(1:nlakes)
            end if
         end if

         do ilake = 1, nlakes
           if (lake_t_ncell_dp(ilake) .gt. 0) then
              rt_domain(did)%LakeOutflowTemp(ilake) = lake_t_sum_dp(ilake) / real(lake_t_ncell_dp(ilake))
           else if (lake_t_ncell(ilake) .gt. 0) then
              rt_domain(did)%LakeOutflowTemp(ilake) = lake_t_sum(ilake) / real(lake_t_ncell(ilake))
           end if
#ifdef HYDRO_D
           print *, "LakeOutflowTemp ilake z_out T", ilake, &
                rt_domain(did)%LakeOutletDepth(ilake), rt_domain(did)%LakeOutflowTemp(ilake)
#endif
         end do


  end subroutine lake1d_step

!===============================================================================
! Heat carried into each 1-D lake by its river inflow.
!
! Over one LSM step the lake receives V = Q_in * dt of water at T_in (the
! flow-weighted temperature of the reaches entering the lake, from the stream
! temperature model).  Spread uniformly over the lake surface area A (LAKEPARM
! LkArea), this is a layer of thickness h = V / A in every column.  The inflow
! enters at its level of neutral buoyancy: the first ice-free layer whose
! density is at least that of the inflow (the bottom layer if none is).  The
! level-pool reservoir carries the lake water balance, so an equal volume
! leaves through the outlet: water between the insertion layer and the outlet
! layer (layer containing LakeOutletDepth) is displaced by h toward the outlet
! (first-order upwind, in passes of at most one layer thickness).  The heat
! change of each column is therefore h * (T_in - T_outlet layer), and outflow
! leaves at the outlet temperature.  If h exceeds the ice-free column depth the
! column is flushed to T_in.
!===============================================================================
  subroutine lake_inflow_heat(xs, xe, ys, ye, nlakes, lake_ix, lakemask, dt, q_in, t_in, &
                              area_km2, z_outlet, z_lake3d, dz_lake3d, icefrac3d, t_lake3d)
    use module_sf_lake, only: nlevlake
    implicit none
    integer, intent(in) :: xs, xe, ys, ye, nlakes
    integer, intent(in) :: lake_ix(xs:xe,ys:ye)
    real,    intent(in) :: lakemask(xs:xe,ys:ye)
    real,    intent(in) :: dt                      ! LSM step (s)
    real,    intent(in) :: q_in(nlakes)            ! lake inflow (m3/s)
    real,    intent(in) :: t_in(nlakes)            ! lake inflow temperature (K)
    real,    intent(in) :: area_km2(nlakes)        ! lake surface area (km2)
    real,    intent(in) :: z_outlet(nlakes)        ! outlet depth below the surface (m)
    real,    intent(in) :: z_lake3d(xs:xe,1:nlevlake,ys:ye), dz_lake3d(xs:xe,1:nlevlake,ys:ye)
    real,    intent(in) :: icefrac3d(xs:xe,1:nlevlake,ys:ye)
    real,    intent(inout) :: t_lake3d(xs:xe,1:nlevlake,ys:ye)

    integer :: i, j, k, lk, kf, k_in, k_out, kmin, kmax, npass, ip
    real    :: h, hp, tin, rho_in, dzmin

    do j = ys, ye
       do i = xs, xe
          lk = lake_ix(i,j)
          if (lk .le. 0) cycle
          if (lakemask(i,j) .ne. 1) cycle
          if (q_in(lk) .le. 0.0 .or. area_km2(lk) .le. 0.0) cycle
          if (t_in(lk) .lt. 200.0 .or. t_in(lk) .gt. 400.0) cycle      ! unset / invalid
          h   = q_in(lk) * dt / (area_km2(lk) * 1.0e6)
          tin = max(273.15, t_in(lk))

          ! first ice-free layer
          kf = 0
          do k = 1, nlevlake
             if (icefrac3d(i,k,j) .le. 0.0) then
                kf = k
                exit
             end if
          end do
          if (kf .eq. 0) cycle

          ! flushed column
          if (h .ge. sum(dz_lake3d(i,kf:nlevlake,j))) then
             t_lake3d(i,kf:nlevlake,j) = tin
             cycle
          end if

          ! insertion layer (neutral buoyancy)
          rho_in = lake_water_density(tin)
          k_in = nlevlake
          do k = kf, nlevlake
             if (lake_water_density(t_lake3d(i,k,j)) .ge. rho_in) then
                k_in = k
                exit
             end if
          end do

          ! outlet layer
          k_out = nlevlake
          do k = 1, nlevlake
             if (z_outlet(lk) .le. z_lake3d(i,k,j) + 0.5 * dz_lake3d(i,k,j)) then
                k_out = k
                exit
             end if
          end do
          k_out = max(k_out, kf)

          kmin  = min(k_in, k_out)
          kmax  = max(k_in, k_out)
          dzmin = max(1.0e-3, minval(dz_lake3d(i,kmin:kmax,j)))
          npass = max(1, ceiling(h / dzmin))
          hp    = h / real(npass)
          do ip = 1, npass
             if (k_in .lt. k_out) then          ! displace downward toward the outlet
                do k = k_out, k_in+1, -1
                   t_lake3d(i,k,j) = t_lake3d(i,k,j) + hp / dz_lake3d(i,k,j) * (t_lake3d(i,k-1,j) - t_lake3d(i,k,j))
                end do
             else if (k_in .gt. k_out) then     ! displace upward toward the outlet
                do k = k_out, k_in-1
                   t_lake3d(i,k,j) = t_lake3d(i,k,j) + hp / dz_lake3d(i,k,j) * (t_lake3d(i,k+1,j) - t_lake3d(i,k,j))
                end do
             end if
             t_lake3d(i,k_in,j) = t_lake3d(i,k_in,j) + hp / dz_lake3d(i,k_in,j) * (tin - t_lake3d(i,k_in,j))
          end do
       end do
    end do
  end subroutine lake_inflow_heat

  ! Fresh-water density (kg/m3) as in the 1-D lake model (maximum at 277 K).
  pure real function lake_water_density(t)
    implicit none
    real, intent(in) :: t
    lake_water_density = 1000.0 * (1.0 - 1.9549e-5 * abs(t - 277.0)**1.68)
  end function lake_water_density

#ifdef WRF_HYDRO
!===============================================================================
! Coupled use (NUOPC cap).  With an external land model the offline Noah-MP
! driver is not used, so the lake model gets its land-surface fields from the
! setup (wrfinput) file and its forcing from the coupler.
!===============================================================================

!-------------------------------------------------------------------------------
! Read the land-surface fields of the setup file as the HRLDAS driver does
! (READLAND_HRLDAS with the lake model on), and cold-start the lake model
! (lakeini with the setup-file skin temperature, as an HRLDAS cold start).
! XSTART:XEND, YSTART:YEND: global indices of this rank's LSM tile.
!-------------------------------------------------------------------------------
  subroutine lake1d_cpl_init(XSTART, XEND, YSTART, YEND, did)
    implicit none
    integer, intent(in) :: XSTART, XEND, YSTART, YEND, did
    character(len=256)  :: flnm
    real, allocatable   :: tmp(:,:), gvg(:,:)
    integer, allocatable :: gid(:,:)
    integer :: ncid, dimid, nxg, nyg, ierr, i, j, ii, jj, k
    logical :: from_rst

    flnm = noah_lsm%hrldas_setup_file
    allocate(tmp(XSTART:XEND,YSTART:YEND))
    allocate(cpl_ivgtyp(XSTART:XEND,YSTART:YEND), cpl_isltyp(XSTART:XEND,YSTART:YEND))
    allocate(cpl_terrain(XSTART:XEND,YSTART:YEND), cpl_xland(XSTART:XEND,YSTART:YEND),  &
             cpl_xice(XSTART:XEND,YSTART:YEND),    cpl_tsk(XSTART:XEND,YSTART:YEND),    &
             cpl_snow(XSTART:XEND,YSTART:YEND),    cpl_xlat(XSTART:XEND,YSTART:YEND),   &
             cpl_albedo(XSTART:XEND,YSTART:YEND),  cpl_emiss(XSTART:XEND,YSTART:YEND),  &
             cpl_hfx(XSTART:XEND,YSTART:YEND),     cpl_lh(XSTART:XEND,YSTART:YEND),     &
             cpl_grdflx(XSTART:XEND,YSTART:YEND),  cpl_qfx(XSTART:XEND,YSTART:YEND),    &
             cpl_glw(XSTART:XEND,YSTART:YEND),     cpl_swdown(XSTART:XEND,YSTART:YEND), &
             cpl_rainbl(XSTART:XEND,YSTART:YEND))
    allocate(cpl_t_phy(XSTART:XEND,1:2,YSTART:YEND), cpl_p8w(XSTART:XEND,1:2,YSTART:YEND),  &
             cpl_dz8w(XSTART:XEND,1:2,YSTART:YEND),  cpl_qv(XSTART:XEND,1:2,YSTART:YEND),   &
             cpl_u(XSTART:XEND,1:2,YSTART:YEND),     cpl_v(XSTART:XEND,1:2,YSTART:YEND))

    ierr = nf90_open(trim(flnm), NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_cpl_init: cannot open setup file '//trim(flnm))
    ierr = nf90_get_att(ncid, NF90_GLOBAL, 'ISWATER', cpl_iswater)
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_cpl_init: no ISWATER attribute in '//trim(flnm))
    ierr = nf90_get_att(ncid, NF90_GLOBAL, 'ISLAKE', cpl_islake)
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_cpl_init: no ISLAKE attribute in '//trim(flnm))
    ierr = nf90_inq_dimid(ncid, 'west_east', dimid);   ierr = nf90_inquire_dimension(ncid, dimid, len=nxg)
    ierr = nf90_inq_dimid(ncid, 'south_north', dimid); ierr = nf90_inquire_dimension(ncid, dimid, len=nyg)
    ierr = nf90_close(ncid)

    call lake1d_read_field(flnm, 'IVGTYP', XSTART, XEND, YSTART, YEND, tmp); cpl_ivgtyp = nint(tmp)
    call lake1d_read_field(flnm, 'ISLTYP', XSTART, XEND, YSTART, YEND, tmp); cpl_isltyp = nint(tmp)
    call lake1d_read_field(flnm, 'HGT',    XSTART, XEND, YSTART, YEND, cpl_terrain)
    call lake1d_read_field(flnm, 'XLAND',  XSTART, XEND, YSTART, YEND, cpl_xland)
    call lake1d_read_field(flnm, 'XLAT',   XSTART, XEND, YSTART, YEND, cpl_xlat)
    call lake1d_read_field(flnm, 'TSK',    XSTART, XEND, YSTART, YEND, cpl_tsk)
    call lake1d_read_field(flnm, 'SNOW',   XSTART, XEND, YSTART, YEND, cpl_snow)
    ! as READLAND_HRLDAS with the lake model on
    where (cpl_ivgtyp == cpl_iswater) cpl_isltyp = 14
    where (cpl_isltyp == 14 .and. cpl_ivgtyp /= cpl_iswater) cpl_isltyp = 8
    cpl_xice = 0.0

    call lake1d_allocate(XSTART, XEND, YSTART, YEND)
    call lake1d_setup(XSTART, XEND, YSTART, YEND, cpl_ivgtyp, cpl_iswater, cpl_islake)
    ! Lake state from LAKE_RESTART_FILE (namelist.hrldas): a LAKE_RST file or
    ! an HRLDAS RESTART file with the lake state; cold start if not set.  The
    ! land model holds the land state, so RESTART_FILENAME_REQUESTED is unused.
    if (len_trim(noah_lsm%restart_filename_requested) > 0) &
       write(6,*) 'WARNING: lake1d_cpl_init: RESTART_FILENAME_REQUESTED is not used in coupled ', &
                  'runs; the 1-D lake model state is read from LAKE_RESTART_FILE'
    from_rst = len_trim(noah_lsm%lake_restart_file) > 0
    if (from_rst) then
       if (.not. lake1d_read_restart(trim(noah_lsm%lake_restart_file), XSTART, XEND, YSTART, YEND)) &
          call lake1d_fatal('lake1d_cpl_init: no complete 1-D lake model state in LAKE_RESTART_FILE '// &
                            trim(noah_lsm%lake_restart_file))
    end if
    call lake1d_init(from_rst, cpl_ivgtyp, cpl_isltyp, cpl_terrain, cpl_snow, from_rst,   &
                     cpl_iswater, cpl_xice, 0.5, cpl_xland, cpl_tsk,                        &
                     XSTART, XEND, YSTART, YEND, XSTART, XEND, YSTART, YEND)
    cpl_albedo = 0.0
    cpl_emiss  = 0.0
    where (cpl_ivgtyp == cpl_islake) cpl_emiss  = 0.98   ! as the HRLDAS driver
    where (cpl_ivgtyp == cpl_islake) cpl_albedo = 0.08
    cpl_hfx = 0.0; cpl_lh = 0.0; cpl_grdflx = 0.0; cpl_qfx = 0.0
    if (from_rst) call lake1d_read_restart_sfc(trim(noah_lsm%lake_restart_file), &
                         XSTART, XEND, YSTART, YEND, cpl_tsk, cpl_hfx, cpl_lh, cpl_qfx,  &
                         cpl_grdflx, cpl_albedo, cpl_emiss)
    call lake1d_allocate_lakes(rt_domain(did)%nlakes)

    ! Lake groups for the lake-cell forcing: one per LAKE_ID on the domain,
    ! whether or not the lake is a routing reservoir (LAKEPARM), numbered in
    ! increasing LAKE_ID order (the same on every rank).  cpl_lake_ix: group of
    ! a lake-model (ISLAKE) cell, -1 if it has no LAKE_ID; cpl_shore_ix: group
    ! of the lake next to a land cell (8 neighbors), 0 if none.
    allocate(gid(nxg,nyg), gvg(nxg,nyg))
    call read_lake_id(flnm, 1, nxg, 1, nyg, gid)
    call lake1d_read_field(flnm, 'IVGTYP', 1, nxg, 1, nyg, gvg)
    cpl_nlkgrp = 0
    k = 0
    do
       ii = minval(gid, mask = gid > k)            ! next larger LAKE_ID
       if (ii == huge(ii) .or. ii <= k) exit
       cpl_nlkgrp = cpl_nlkgrp + 1
       where (gid == ii) gid = -cpl_nlkgrp - 1     ! group g stored as -(g+1) until done
       k = 0
    end do
    where (gid < -1) gid = -gid - 1                ! now gid = group (1..), others <= 0
    where (gid < 0)  gid = 0
    allocate(cpl_shore_ix(XSTART:XEND,YSTART:YEND), cpl_lake_ix(XSTART:XEND,YSTART:YEND))
    cpl_shore_ix = 0
    cpl_lake_ix  = 0
    do j = YSTART, YEND
       do i = XSTART, XEND
          if (cpl_ivgtyp(i,j) == cpl_islake) then
             cpl_lake_ix(i,j) = gid(i,j)
             if (gid(i,j) == 0) cpl_lake_ix(i,j) = -1
          else if (gid(i,j) == 0 .and. nint(gvg(i,j)) /= cpl_iswater .and. &
                   nint(gvg(i,j)) /= cpl_islake) then
             nbr: do jj = max(1,j-1), min(nyg,j+1)
                do ii = max(1,i-1), min(nxg,i+1)
                   if (gid(ii,jj) > 0) then
                      cpl_shore_ix(i,j) = gid(ii,jj)
                      exit nbr
                   end if
                end do
             end do nbr
          end if
       end do
    end do
    deallocate(gid, gvg, tmp)
    if (from_rst) then
       write(6,*) 'lake1d_cpl_init: 1-D lake model state from ', &
                  trim(noah_lsm%lake_restart_file)
    else
       write(6,*) 'lake1d_cpl_init: 1-D lake model cold start'
    end if
    write(6,*) 'lake1d_cpl_init: lake cells on this tile: ', &
               count(cpl_lake_ix > 0), ', shoreline cells: ', count(cpl_shore_ix > 0)

  end subroutine lake1d_cpl_init

!-------------------------------------------------------------------------------
! Lake cells of each lake (LAKE_ID) get the mean of a forcing field over that
! lake's shoreline land cells with valid values (summed across ranks).  Lake
! cells without a LAKE_ID, or of a lake without shoreline land cells, get the
! mean over all land cells of the domain.
!-------------------------------------------------------------------------------
  subroutine lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, fld)
    implicit none
    integer, intent(in)    :: XSTART, XEND, YSTART, YEND
    real,    intent(inout) :: fld(XSTART:XEND,YSTART:YEND)
    real    :: ssum(2*cpl_nlkgrp+2)
    logical :: land(XSTART:XEND,YSTART:YEND)
    integer :: k, n
    ! (collective: every rank takes part in the sum, lake cells or not)
    n = cpl_nlkgrp
    ! shoreline cells next to cells the land model does not cover can carry the
    ! coupler's missing value through the regridding; use valid values only
    land = cpl_ivgtyp /= cpl_iswater .and. cpl_ivgtyp /= cpl_islake .and. abs(fld) < 1.0e15
    do k = 1, n
       ssum(k)   = sum(fld, mask = land .and. cpl_shore_ix == k)
       ssum(n+k) = real(count(land .and. cpl_shore_ix == k))
    end do
    ssum(2*n+1) = sum(fld, mask = land)
    ssum(2*n+2) = real(count(land))
#ifdef MPP_LAND
    call sum_real1d(ssum, 2*n+2)
#endif
    do k = 1, n
       if (ssum(n+k) > 0.0) then
          where (cpl_lake_ix == k) fld = ssum(k) / ssum(n+k)
       else if (ssum(2*n+2) > 0.0) then
          where (cpl_lake_ix == k) fld = ssum(2*n+1) / ssum(2*n+2)
       end if
    end do
    if (ssum(2*n+2) > 0.0) then
       where (cpl_lake_ix == -1) fld = ssum(2*n+1) / ssum(2*n+2)
    end if
  end subroutine lake1d_fill_lake_cells

!-------------------------------------------------------------------------------
! One coupled step: lake-cell forcing from the shoreline, the forcing arrays as
! the HRLDAS driver prepares them, the shared lake step (lake1d_step), and the
! lake precipitation/evaporation hand-off to the level-pool reservoirs (as
! hrldas_cpl_HYDRO).  Forcing fields on the LSM tile (they are modified on lake
! cells): air temperature [K], specific humidity [kg/kg], wind speed [m/s],
! surface pressure [Pa], downward longwave and shortwave [W/m2], precipitation
! rate [kg/m2/s].  dt: LSM time step [s].
!-------------------------------------------------------------------------------
  subroutine lake1d_cpl_step(XSTART, XEND, YSTART, YEND, did, dt, t2, q2x, wind, sfcp, lwdown, short, prcp)
    implicit none
    integer, intent(in)    :: XSTART, XEND, YSTART, YEND, did
    real,    intent(in)    :: dt
    real, dimension(XSTART:XEND,YSTART:YEND), intent(inout) :: t2, q2x, wind, sfcp, lwdown, short, prcp
    integer :: k

    call lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, t2)
    call lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, q2x)
    call lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, wind)
    call lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, sfcp)
    call lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, lwdown)
    call lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, short)
    call lake1d_fill_lake_cells(XSTART, XEND, YSTART, YEND, prcp)

    do k = 1, 2
       cpl_t_phy(:,k,:) = t2
       cpl_qv(:,k,:)    = q2x / (1.0 - q2x)      ! mixing ratio, as the HRLDAS driver
       cpl_u(:,k,:)     = wind
       cpl_v(:,k,:)     = 0.0
       cpl_p8w(:,k,:)   = sfcp
    end do
    cpl_dz8w   = 2 * noah_lsm%zlvl              ! as the HRLDAS driver
    cpl_glw    = lwdown
    cpl_swdown = short
    cpl_rainbl = prcp * dt                      ! [mm per step]

    call lake1d_step(XSTART, XEND, YSTART, YEND, dt,                                         &
                     cpl_t_phy, cpl_p8w, cpl_dz8w, cpl_qv, cpl_u, cpl_v, cpl_glw, cpl_emiss,  &
                     cpl_rainbl, cpl_swdown, cpl_albedo, cpl_xlat, cpl_ivgtyp, cpl_terrain,  &
                     cpl_xland, cpl_iswater, cpl_xice, 0.5,                                  &
                     XSTART, XEND, YSTART, YEND, 1, 2, XSTART, XEND, YSTART, YEND, 1, 2,     &
                     XSTART, XEND, YSTART, YEND, 1, 2,                                       &
                     cpl_hfx, cpl_lh, cpl_grdflx, cpl_tsk, cpl_qfx)

    ! hand-off to the level-pool reservoirs, as hrldas_cpl_HYDRO
    if (nlakes > 0) then
       where (lake_grids > 0) lake_precip = lake_precip / lake_grids * 3600.0 / nlst(did)%dt
       where (lake_grids > 0) lake_evap   = lake_evap / lake_grids
       rt_domain(did)%lake_precip     = lake_precip
       rt_domain(did)%lake_evap       = lake_evap
       rt_domain(did)%lake_netoutflow = lake_netoutflow
    end if
  end subroutine lake1d_cpl_step

!-------------------------------------------------------------------------------
! After routing, as hrldas_cpl_HYDRO
!-------------------------------------------------------------------------------
  subroutine lake1d_cpl_post(did)
    implicit none
    integer, intent(in) :: did
    if (nlakes > 0) then
       lake_netoutflow = rt_domain(did)%lake_netoutflow
       lake_evap       = rt_domain(did)%lake_evap
    end if
  end subroutine lake1d_cpl_post

!-------------------------------------------------------------------------------
! Lake model restart (LAKE_RST.<date>_DOMAIN<n>), offline and coupled.
!
! The lake state (lake_rst_vars) has the variable names and dimensions of the
! lake state in the offline HRLDAS RESTART file of earlier versions, so the
! lake model can also restart from such a RESTART file.  LAKE_RST files also
! hold the lake surface fields that a coupled atmosphere receives over lake
! cells (lake_rst_sfc_vars), so that the first exchange after a restart is
! that of a continuous run; they are optional on input.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read the lake state on this tile.  Returns .false., reading nothing, if the
! file cannot be opened or lacks any lake state variable.
!-------------------------------------------------------------------------------
  logical function lake1d_read_restart(flnm, XSTART, XEND, YSTART, YEND) result(ok)
    implicit none
    character(len=*), intent(in) :: flnm
    integer,          intent(in) :: XSTART, XEND, YSTART, YEND
    integer :: ncid, varid, ierr, k

    ok = .false.
    ierr = nf90_open(flnm, NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) then
       write(6,*) 'WARNING: lake1d_read_restart: cannot open ', trim(flnm)
       return
    end if
    do k = 1, size(lake_rst_vars)
       if (nf90_inq_varid(ncid, trim(lake_rst_vars(k)), varid) /= NF90_NOERR) then
          write(6,*) 'WARNING: lake1d_read_restart: no ', trim(lake_rst_vars(k)), ' in ', trim(flnm)
          ierr = nf90_close(ncid)
          return
       end if
    end do
    call lake1d_get2d(ncid, flnm, 'LAKE_DEPTH', XSTART, XEND, YSTART, YEND, lakedepth2d)
    call lake1d_get2d(ncid, flnm, 'LAKE_SAVETKE1', XSTART, XEND, YSTART, YEND, savedtke12d)
    call lake1d_get2d(ncid, flnm, 'LAKE_SNOWP', XSTART, XEND, YSTART, YEND, snowdp2d)
    call lake1d_get2d(ncid, flnm, 'LAKE_H2OSNO', XSTART, XEND, YSTART, YEND, h2osno2d)
    call lake1d_get2d(ncid, flnm, 'LAKE_SNL', XSTART, XEND, YSTART, YEND, snl2d)
    call lake1d_get2d(ncid, flnm, 'LAKE_T_GRND', XSTART, XEND, YSTART, YEND, t_grnd2d)
    call lake1d_get3d(ncid, flnm, 'LAKE_T_SOILSNOW', XSTART, YSTART, t_soisno3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_H2OSOI_ICE', XSTART, YSTART, h2osoi_ice3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_H2OSOI_LIQ', XSTART, YSTART, h2osoi_liq3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_H2OSOI_VOL', XSTART, YSTART, h2osoi_vol3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_Z', XSTART, YSTART, z3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_DZ', XSTART, YSTART, dz3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_T_LAKE', XSTART, YSTART, t_lake3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_ICEFRAC', XSTART, YSTART, lake_icefrac3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_Z_LAKE', XSTART, YSTART, z_lake3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_DZ_LAKE', XSTART, YSTART, dz_lake3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_WATSAT', XSTART, YSTART, watsat3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_CSOL', XSTART, YSTART, csol3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_TKMG', XSTART, YSTART, tkmg3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_TKDRY', XSTART, YSTART, tkdry3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_TKSATU', XSTART, YSTART, tksatu3d)
    call lake1d_get3d(ncid, flnm, 'LAKE_ZI', XSTART, YSTART, zi3d)
    ierr = nf90_close(ncid)
    ok = .true.
  end function lake1d_read_restart

!-------------------------------------------------------------------------------
! Set the lake surface fields on lake cells (lakemask = 1) from a restart file,
! where present (LAKE_RST files; not in HRLDAS RESTART files).  Call after the
! lake model is initialized (lake1d_init sets the lake mask).
!-------------------------------------------------------------------------------
  subroutine lake1d_read_restart_sfc(flnm, XSTART, XEND, YSTART, YEND, &
                                     tsk, hfx, lh, qfx, grdflx, albedo, emiss)
    implicit none
    character(len=*), intent(in) :: flnm
    integer,          intent(in) :: XSTART, XEND, YSTART, YEND
    real, dimension(XSTART:XEND,YSTART:YEND), intent(inout) :: tsk, hfx, lh, qfx, grdflx, albedo, emiss
    real    :: a(XSTART:XEND,YSTART:YEND)
    integer :: ncid, varid, ierr, n
    ierr = nf90_open(flnm, NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) return
    n = 0
    call get('LAKE_TSK', tsk);       call get('LAKE_HFX', hfx);     call get('LAKE_LH', lh)
    call get('LAKE_QFX', qfx);       call get('LAKE_GRDFLX', grdflx)
    call get('LAKE_ALBEDO', albedo); call get('LAKE_EMISS', emiss)
    call get('LAKE_T2', t2);         call get('LAKE_TH2', th2);     call get('LAKE_Q2', q2)
    ierr = nf90_close(ncid)
    if (n > 0) write(6,*) 'lake1d_read_restart_sfc: ', n, ' lake surface fields from ', trim(flnm)
  contains
    subroutine get(name, fld)
      character(len=*), intent(in) :: name
      real, intent(inout) :: fld(XSTART:XEND,YSTART:YEND)
      if (nf90_inq_varid(ncid, name, varid) /= NF90_NOERR) return
      call lake1d_get2d(ncid, flnm, name, XSTART, XEND, YSTART, YEND, a)
      where (lakemask == 1) fld = a
      n = n + 1
    end subroutine get
  end subroutine lake1d_read_restart_sfc

  subroutine lake1d_get2d(ncid, flnm, name, XSTART, XEND, YSTART, YEND, a)
    implicit none                  ! file: (west_east, south_north, Time)
    integer,          intent(in)  :: ncid, XSTART, XEND, YSTART, YEND
    character(len=*), intent(in)  :: flnm, name
    real,             intent(out) :: a(XSTART:XEND,YSTART:YEND)
    integer :: varid, ierr
    ierr = nf90_inq_varid(ncid, name, varid)
    if (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, varid, a, start=(/XSTART,YSTART,1/), &
                                                count=(/XEND-XSTART+1,YEND-YSTART+1,1/))
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d restart: cannot read '//name//' from '//trim(flnm))
  end subroutine lake1d_get2d

  subroutine lake1d_get3d(ncid, flnm, name, XSTART, YSTART, a)
    implicit none                  ! file: (west_east, layers, south_north, Time)
    integer,          intent(in)  :: ncid, XSTART, YSTART
    character(len=*), intent(in)  :: flnm, name
    real,             intent(out) :: a(:,:,:)
    integer :: varid, ierr
    ierr = nf90_inq_varid(ncid, name, varid)
    if (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, varid, a, start=(/XSTART,1,YSTART,1/), &
                                                count=(/size(a,1),size(a,2),size(a,3),1/))
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d restart: cannot read '//name//' from '//trim(flnm))
  end subroutine lake1d_get3d

!-------------------------------------------------------------------------------
! Write LAKE_RST.<date>_DOMAIN<hgrid> (date: YYYY-MM-DD_HH:MM): the lake state
! and the lake surface fields (values on lake cells, 0 elsewhere).  Called when
! WRF-Hydro writes a HYDRO_RST, so the lake and reservoir restarts are paired.
! Collective over the WRF-Hydro ranks; the IO rank writes.
!-------------------------------------------------------------------------------
  subroutine lake1d_write_restart(date, XSTART, XEND, YSTART, YEND, &
                                  tsk, hfx, lh, qfx, grdflx, albedo, emiss)
    implicit none
    character(len=*), intent(in) :: date
    integer,          intent(in) :: XSTART, XEND, YSTART, YEND
    real, dimension(XSTART:XEND,YSTART:YEND), intent(in) :: tsk, hfx, lh, qfx, grdflx, albedo, emiss
    character(len=256) :: flnm
    integer :: ncid, ierr, nxg, nyg, k
    integer :: d_t, d_s, d_x, d_y, d_lk, d_lss, d_ls, d_lss2, vid
    logical :: io

    flnm = 'LAKE_RST.'//date(1:16)//'_DOMAIN'//trim(nlst(1)%hgrid)
#ifdef MPP_LAND
    nxg = global_nx; nyg = global_ny; io = (my_id == IO_id)
#else
    nxg = XEND-XSTART+1; nyg = YEND-YSTART+1; io = .true.
#endif
    if (io) then
       ierr = nf90_create(trim(flnm), ior(NF90_CLOBBER, NF90_64BIT_OFFSET), ncid)
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_write_restart: cannot create '//trim(flnm))
       ierr = nf90_def_dim(ncid, 'Time', NF90_UNLIMITED, d_t)
       ierr = nf90_def_dim(ncid, 'DateStrLen', 19, d_s)
       ierr = nf90_def_dim(ncid, 'west_east', nxg, d_x)
       ierr = nf90_def_dim(ncid, 'south_north', nyg, d_y)
       ierr = nf90_def_dim(ncid, 'lake_layers', nlevlake, d_lk)
       ierr = nf90_def_dim(ncid, 'lakesoilsnow_layers', nlevsnow+nlevsoil, d_lss)
       ierr = nf90_def_dim(ncid, 'lakesoil_layers', nlevsoil, d_ls)
       ierr = nf90_def_dim(ncid, 'lakesoilsnow_layers2', nlevsnow+nlevsoil+1, d_lss2)
       ierr = nf90_def_var(ncid, 'Times', NF90_CHAR, (/d_s,d_t/), vid)
       do k = 1, 6
          ierr = nf90_def_var(ncid, trim(lake_rst_vars(k)), NF90_FLOAT, (/d_x,d_y,d_t/), vid)
       end do
       do k = 7, 12
          ierr = nf90_def_var(ncid, trim(lake_rst_vars(k)), NF90_FLOAT, (/d_x,d_lss,d_y,d_t/), vid)
       end do
       do k = 13, 16
          ierr = nf90_def_var(ncid, trim(lake_rst_vars(k)), NF90_FLOAT, (/d_x,d_lk,d_y,d_t/), vid)
       end do
       do k = 17, 21
          ierr = nf90_def_var(ncid, trim(lake_rst_vars(k)), NF90_FLOAT, (/d_x,d_ls,d_y,d_t/), vid)
       end do
       ierr = nf90_def_var(ncid, 'LAKE_ZI', NF90_FLOAT, (/d_x,d_lss2,d_y,d_t/), vid)
       do k = 1, size(lake_rst_sfc_vars)
          ierr = nf90_def_var(ncid, trim(lake_rst_sfc_vars(k)), NF90_FLOAT, (/d_x,d_y,d_t/), vid)
       end do
       ierr = nf90_put_att(ncid, NF90_GLOBAL, 'TITLE', 'WRF-Hydro 1-D lake model restart')
       ierr = nf90_enddef(ncid)
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_write_restart: cannot define '//trim(flnm))
       ierr = nf90_inq_varid(ncid, 'Times', vid)
       ierr = nf90_put_var(ncid, vid, (/ date(1:16)//':00' /), start=(/1,1/), count=(/19,1/))
    end if

    call lake1d_put2d(ncid, io, 'LAKE_DEPTH', lakedepth2d)
    call lake1d_put2d(ncid, io, 'LAKE_SAVETKE1', savedtke12d)
    call lake1d_put2d(ncid, io, 'LAKE_SNOWP', snowdp2d)
    call lake1d_put2d(ncid, io, 'LAKE_H2OSNO', h2osno2d)
    call lake1d_put2d(ncid, io, 'LAKE_SNL', snl2d)
    call lake1d_put2d(ncid, io, 'LAKE_T_GRND', t_grnd2d)
    call lake1d_put3d(ncid, io, 'LAKE_T_SOILSNOW', t_soisno3d)
    call lake1d_put3d(ncid, io, 'LAKE_H2OSOI_ICE', h2osoi_ice3d)
    call lake1d_put3d(ncid, io, 'LAKE_H2OSOI_LIQ', h2osoi_liq3d)
    call lake1d_put3d(ncid, io, 'LAKE_H2OSOI_VOL', h2osoi_vol3d)
    call lake1d_put3d(ncid, io, 'LAKE_Z', z3d)
    call lake1d_put3d(ncid, io, 'LAKE_DZ', dz3d)
    call lake1d_put3d(ncid, io, 'LAKE_T_LAKE', t_lake3d)
    call lake1d_put3d(ncid, io, 'LAKE_ICEFRAC', lake_icefrac3d)
    call lake1d_put3d(ncid, io, 'LAKE_Z_LAKE', z_lake3d)
    call lake1d_put3d(ncid, io, 'LAKE_DZ_LAKE', dz_lake3d)
    call lake1d_put3d(ncid, io, 'LAKE_WATSAT', watsat3d)
    call lake1d_put3d(ncid, io, 'LAKE_CSOL', csol3d)
    call lake1d_put3d(ncid, io, 'LAKE_TKMG', tkmg3d)
    call lake1d_put3d(ncid, io, 'LAKE_TKDRY', tkdry3d)
    call lake1d_put3d(ncid, io, 'LAKE_TKSATU', tksatu3d)
    call lake1d_put3d(ncid, io, 'LAKE_ZI', zi3d)
    call put_sfc('LAKE_TSK', tsk);       call put_sfc('LAKE_HFX', hfx)
    call put_sfc('LAKE_LH', lh);         call put_sfc('LAKE_QFX', qfx)
    call put_sfc('LAKE_GRDFLX', grdflx); call put_sfc('LAKE_ALBEDO', albedo)
    call put_sfc('LAKE_EMISS', emiss);   call put_sfc('LAKE_T2', t2)
    call put_sfc('LAKE_TH2', th2);       call put_sfc('LAKE_Q2', q2)

    if (io) then
       ierr = nf90_close(ncid)
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_write_restart: cannot write '//trim(flnm))
       write(6,*) 'lake1d_write_restart: wrote ', trim(flnm)
    end if
#ifdef MPP_LAND
    call mpp_land_sync()
#endif
  contains
    subroutine put_sfc(name, fld)
      character(len=*), intent(in) :: name
      real, intent(in) :: fld(XSTART:XEND,YSTART:YEND)
      real :: a(XSTART:XEND,YSTART:YEND)
      a = 0.0
      where (lakemask == 1) a = fld
      call lake1d_put2d(ncid, io, name, a)
    end subroutine put_sfc
  end subroutine lake1d_write_restart

!-------------------------------------------------------------------------------
! Coupled (NUOPC cap): LAKE_RST with the coupled lake surface fields
!-------------------------------------------------------------------------------
  subroutine lake1d_cpl_write_restart(date, XSTART, XEND, YSTART, YEND)
    implicit none
    character(len=*), intent(in) :: date
    integer,          intent(in) :: XSTART, XEND, YSTART, YEND
    call lake1d_write_restart(date, XSTART, XEND, YSTART, YEND, &
                              cpl_tsk, cpl_hfx, cpl_lh, cpl_qfx, cpl_grdflx, cpl_albedo, cpl_emiss)
  end subroutine lake1d_cpl_write_restart

!-------------------------------------------------------------------------------
! Gather a tile field (XSTART:XEND, [levels,] YSTART:YEND) to the IO rank and
! write it to variable NAME of the open file NCID (IO rank only; the first
! record of a (west_east, [levels,] south_north, Time) variable).  Collective
! over the WRF-Hydro ranks.
!-------------------------------------------------------------------------------
  subroutine lake1d_put2d(ncid, io, name, a)
    implicit none
    integer,          intent(in) :: ncid
    logical,          intent(in) :: io
    character(len=*), intent(in) :: name
    real,             intent(in) :: a(:,:)
    real, allocatable :: g(:,:)
    integer :: varid, ierr
#ifdef MPP_LAND
    allocate(g(global_nx,global_ny))
    call write_IO_real(a, g)
#else
    allocate(g(size(a,1),size(a,2)))
    g = a
#endif
    if (io) then
       ierr = nf90_inq_varid(ncid, name, varid)
       if (ierr == NF90_NOERR) ierr = nf90_put_var(ncid, varid, g, start=(/1,1,1/), &
                                                   count=(/size(g,1),size(g,2),1/))
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_put2d: cannot write '//name)
    end if
    deallocate(g)
  end subroutine lake1d_put2d

  subroutine lake1d_put3d(ncid, io, name, a)
    implicit none
    integer,          intent(in) :: ncid
    logical,          intent(in) :: io
    character(len=*), intent(in) :: name
    real,             intent(in) :: a(:,:,:)
    real, allocatable :: g(:,:), g3(:,:,:)
    integer :: varid, ierr, k
#ifdef MPP_LAND
    allocate(g(global_nx,global_ny))
#else
    allocate(g(size(a,1),size(a,3)))
#endif
    allocate(g3(size(g,1),size(a,2),size(g,2)))
    do k = 1, size(a,2)
#ifdef MPP_LAND
       call write_IO_real(a(:,k,:), g)
#else
       g = a(:,k,:)
#endif
       g3(:,k,:) = g
    end do
    if (io) then
       ierr = nf90_inq_varid(ncid, name, varid)
       if (ierr == NF90_NOERR) ierr = nf90_put_var(ncid, varid, g3, start=(/1,1,1,1/), &
                                                   count=(/size(g3,1),size(g3,2),size(g3,3),1/))
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_put3d: cannot write '//name)
    end if
    deallocate(g, g3)
  end subroutine lake1d_put3d

!-------------------------------------------------------------------------------
! Lake model output of a coupled run, as the lake fields of the offline LDASOUT
! file: instantaneous values every OUTPUT_TIMESTEP (namelist.hrldas) in
! YYYYMMDDHHMM.LAKE1D_DOMAIN<hgrid>, -9999 on cells that are not lake cells.
! Called after each coupled step; DATE (YYYY-MM-DD_HH:MM) is the end of the step.
!-------------------------------------------------------------------------------
  subroutine lake1d_cpl_write_output(did, date, XSTART, XEND, YSTART, YEND)
    implicit none
    integer,          intent(in) :: did
    character(len=*), intent(in) :: date
    integer,          intent(in) :: XSTART, XEND, YSTART, YEND
    character(len=256) :: flnm
    integer :: ncid, ierr, nxg, nyg, d_t, d_s, d_x, d_y, d_lk, vid, k
    logical :: io
    real, allocatable :: w2(:,:), w3(:,:,:)
    real, parameter :: fillv = -9999.0

    if (noah_lsm%output_timestep <= 0) return
    cpl_out_seconds = cpl_out_seconds + nint(nlst(did)%dt)
    if (mod(cpl_out_seconds, noah_lsm%output_timestep) /= 0) return

    flnm = date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16)// &
           '.LAKE1D_DOMAIN'//trim(nlst(did)%hgrid)
#ifdef MPP_LAND
    nxg = global_nx; nyg = global_ny; io = (my_id == IO_id)
#else
    nxg = XEND-XSTART+1; nyg = YEND-YSTART+1; io = .true.
#endif
    if (io) then
       ierr = nf90_create(trim(flnm), ior(NF90_CLOBBER, NF90_64BIT_OFFSET), ncid)
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_cpl_write_output: cannot create '//trim(flnm))
       ierr = nf90_def_dim(ncid, 'Time', NF90_UNLIMITED, d_t)
       ierr = nf90_def_dim(ncid, 'DateStrLen', 19, d_s)
       ierr = nf90_def_dim(ncid, 'west_east', nxg, d_x)
       ierr = nf90_def_dim(ncid, 'south_north', nyg, d_y)
       ierr = nf90_def_dim(ncid, 'lake_layers', nlevlake, d_lk)
       ierr = nf90_def_var(ncid, 'Times', NF90_CHAR, (/d_s,d_t/), vid)
       call defv('LAKE_PRECIP2D', (/d_x,d_y,d_t/), 'mm h-1', 'precipitation over lakes')
       call defv('LAKE_EVAP2D',   (/d_x,d_y,d_t/), 'mm h-1', 'evaporation from lakes')
       call defv('LAKE_DEPTH2D',  (/d_x,d_y,d_t/), 'm', 'lake depth')
       call defv('LAKE_TSK2D',    (/d_x,d_y,d_t/), 'K', 'lake surface (skin) temperature')
       call defv('LAKE_T3D',      (/d_x,d_lk,d_y,d_t/), 'K', 'lake layer temperature')
       call defv('LAKE_ICE3D',    (/d_x,d_lk,d_y,d_t/), '-', 'lake layer ice fraction')
       call defv('LAKE_Z3D',      (/d_x,d_lk,d_y,d_t/), 'm', 'lake layer node depth')
       ierr = nf90_put_att(ncid, NF90_GLOBAL, 'TITLE', '1-D lake model output (WRF-Hydro NUOPC cap)')
       ierr = nf90_enddef(ncid)
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_cpl_write_output: cannot define '//trim(flnm))
       ierr = nf90_inq_varid(ncid, 'Times', vid)
       ierr = nf90_put_var(ncid, vid, (/ date(1:16)//':00' /), start=(/1,1/), count=(/19,1/))
    end if

    allocate(w2(XSTART:XEND,YSTART:YEND), w3(XSTART:XEND,nlevlake,YSTART:YEND))
    w2 = fillv; where (lakemask == 1) w2 = cpl_rainbl * 3600.0 / nlst(did)%dt
    call lake1d_put2d(ncid, io, 'LAKE_PRECIP2D', w2)
    w2 = fillv; where (lakemask == 1) w2 = cpl_qfx * 3600.0
    call lake1d_put2d(ncid, io, 'LAKE_EVAP2D', w2)
    w2 = fillv; where (lakemask == 1) w2 = lakedepth2d
    call lake1d_put2d(ncid, io, 'LAKE_DEPTH2D', w2)
    w2 = fillv; where (lakemask == 1) w2 = cpl_tsk
    call lake1d_put2d(ncid, io, 'LAKE_TSK2D', w2)
    do k = 1, nlevlake
       w3(:,k,:) = fillv; where (lakemask == 1) w3(:,k,:) = t_lake3d(:,k,:)
    end do
    call lake1d_put3d(ncid, io, 'LAKE_T3D', w3)
    do k = 1, nlevlake
       w3(:,k,:) = fillv; where (lakemask == 1) w3(:,k,:) = lake_icefrac3d(:,k,:)
    end do
    call lake1d_put3d(ncid, io, 'LAKE_ICE3D', w3)
    do k = 1, nlevlake
       w3(:,k,:) = fillv; where (lakemask == 1) w3(:,k,:) = z_lake3d(:,k,:)
    end do
    call lake1d_put3d(ncid, io, 'LAKE_Z3D', w3)
    deallocate(w2, w3)

    if (io) then
       ierr = nf90_close(ncid)
       if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_cpl_write_output: cannot write '//trim(flnm))
    end if
#ifdef MPP_LAND
    call mpp_land_sync()
#endif

  contains

    subroutine defv(name, dims, units, desc)
      character(len=*), intent(in) :: name, units, desc
      integer,          intent(in) :: dims(:)
      ierr = nf90_def_var(ncid, name, NF90_FLOAT, dims, vid)
      ierr = nf90_put_att(ncid, vid, 'units', units)
      ierr = nf90_put_att(ncid, vid, 'description', desc)
      ierr = nf90_put_att(ncid, vid, '_FillValue', fillv)
    end subroutine defv

  end subroutine lake1d_cpl_write_output

!-------------------------------------------------------------------------------
! A 2-D field (optionally with a Time dimension) of the setup file on
! XSTART:XEND, YSTART:YEND
!-------------------------------------------------------------------------------
  subroutine lake1d_read_field(flnm, name, XSTART, XEND, YSTART, YEND, fld)
    implicit none
    character(len=*), intent(in)  :: flnm, name
    integer,          intent(in)  :: XSTART, XEND, YSTART, YEND
    real,             intent(out) :: fld(XSTART:XEND,YSTART:YEND)
    integer :: ncid, varid, ndims, ierr
    integer :: start(3), cnt(3)
    ierr = nf90_open(trim(flnm), NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_read_field: cannot open '//trim(flnm))
    ierr = nf90_inq_varid(ncid, trim(name), varid)
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_read_field: no '//trim(name)//' in '//trim(flnm))
    ierr = nf90_inquire_variable(ncid, varid, ndims=ndims)
    start = (/ XSTART, YSTART, 1 /)
    cnt   = (/ XEND-XSTART+1, YEND-YSTART+1, 1 /)
    ierr = nf90_get_var(ncid, varid, fld, start=start(1:ndims), count=cnt(1:ndims))
    if (ierr /= NF90_NOERR) call lake1d_fatal('lake1d_read_field: cannot read '//trim(name)//' from '//trim(flnm))
    ierr = nf90_close(ncid)
  end subroutine lake1d_read_field

  subroutine lake1d_fatal(msg)
    implicit none
    character(len=*), intent(in) :: msg
    write(6,*) 'FATAL ERROR: '//trim(msg)
    call flush(6)
#ifdef MPP_LAND
    call fatal_error_stop('FATAL ERROR: '//trim(msg))
#else
    stop 1
#endif
  end subroutine lake1d_fatal
#endif

!-------------------------------------------------------------------------------
! LAKE_DEPTH and LAKE_ID from the setup (wrfinput) file, on the tile
!-------------------------------------------------------------------------------
  subroutine read_lake_depth(wrfinput_flnm,xstart,xend,ystart,yend,lake_depth)
    implicit none
    character(len=*),          intent(in)  :: wrfinput_flnm
    integer,                   intent(in)  :: xstart, xend, ystart, yend
    real,    dimension(xstart:xend,  ystart:yend), intent(out) :: lake_depth
    
    character(len=24)  :: name
    integer :: ierr,iret
    integer :: ncid, varid

    ! Open the NetCDF file.
    write(*,'("wrfinput_flnm: ''", A, "''")') trim(wrfinput_flnm)
    ierr = nf90_open(wrfinput_flnm, NF90_NOWRITE, ncid)
    if (ierr /= 0) then
       write(*,'("read_lake_depth:  Problem opening wrfinput file: ''", A, "''")') trim(wrfinput_flnm)
       stop
    endif

! Get lake depth (LAKE_DEPTH)
    name = "LAKE_DEPTH"
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret == 0) then
      ierr = nf90_get_var(ncid, varid, lake_depth, start=(/xstart,ystart/), count=(/xend-xstart+1,yend-ystart+1/))
    else
      stop "Problem retrieving LAKE_DEPTH with use_lake_depth = 1"
    endif

! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= 0) stop "MODULE_NOAHLSM_HRLDAS_INPUT: read_lake_depth :  NF90_CLOSE"


  end subroutine read_lake_depth

!---------------------------------------------------------------------------------------------------------
  subroutine read_lake_id(wrfinput_flnm,xstart,xend,ystart,yend,lake_id)
    implicit none
    character(len=*),          intent(in)  :: wrfinput_flnm
    integer,                   intent(in)  :: xstart, xend, ystart, yend
    integer,    dimension(xstart:xend,  ystart:yend), intent(out) :: lake_id
    
    character(len=24)  :: name
    integer :: ierr,iret
    integer :: ncid, varid

    ! Open the NetCDF file.
    write(*,'("wrfinput_flnm: ''", A, "''")') trim(wrfinput_flnm)
    ierr = nf90_open(wrfinput_flnm, NF90_NOWRITE, ncid)
    if (ierr /= 0) then
       write(*,'("read_lake_id:  Problem opening wrfinput file: ''", A, "''")') trim(wrfinput_flnm)
       stop
    endif

! Get lake ID (LAKE_ID)
    name = "LAKE_ID"
    iret = nf90_inq_varid(ncid,  name,  varid)
    if (iret == 0) then
      ierr = nf90_get_var(ncid, varid, lake_id, start=(/xstart,ystart/), count=(/xend-xstart+1,yend-ystart+1/))
    else
      stop "Problem retrieving LAKE_ID"
    endif

! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= 0) stop "MODULE_NOAHLSM_HRLDAS_INPUT: read_lake_id :  NF90_CLOSE"


  end subroutine read_lake_id

end module module_lake1d
