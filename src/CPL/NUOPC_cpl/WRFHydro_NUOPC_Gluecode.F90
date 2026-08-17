#include "WRFHydro_NUOPC_Macros.h"

module wrfhydro_nuopc_gluecode
! This module connects NUOPC initialize, advance,
! and finalize to WRFHYDRO.
  use ESMF
  use NUOPC
  use mpi
  use WRFHydro_ESMF_Extensions
  use wrfhydro_nuopc_utils, only: array_write_i
  use module_mpp_land, only: &
    HYDRO_COMM_WORLD, &
    numprocs, &
    global_nx, &
    global_ny, &
    startx, &
    starty, &
    local_nx_size, &
    local_ny_size, &
    log_map2d, &
    mpp_land_par_ini, &
    decompose_data_real, &
    write_io_real, my_id, &
    mpp_land_bcast_real1, &
    IO_id, my_id, &
    mpp_land_bcast_real, &
    mpp_land_bcast_int1, &
    MPP_LAND_INIT
  use module_noahmp_hrldas_driver, only: &
       land_driver_ini, land_driver_exe
  use module_HYDRO_drv, only: &
    HYDRO_ini, &
    HYDRO_exe
  use module_HYDRO_io, only: &
    get_file_dimension
  use module_CPL_LAND, only: &
    CPL_LAND_INIT, &
    cpl_outdate
  use module_rt_data, only: &
    rt_domain
  use module_lsm_forcing, only: &
    read_ldasout
  use config_base, only: &
    nlst, &
    init_namelist_rt_field
  use orchestrator_base
  use wrfhydro_nuopc_fields
  use wrfhydro_nuopc_flags
  use wrfhydro_nuopc_utils, only: check, check_nf, printa, create_dir_if_needed

  implicit none

  private

  public :: wrfhydro_nuopc_ini
  public :: wrfhydro_nuopc_run
  public :: wrfhydro_nuopc_fin
  public :: wrfhydro_open_mesh
  public :: wrfhydro_regrid_mesh
  public :: wrfhydro_GridCreate
  public :: wrfhydro_get_timestep
  public :: wrfhydro_set_timestep
  public :: wrfhydro_get_hgrid
  public :: wrfhydro_get_restart

  public :: wrfhydro_write_full_resolution_file
  public :: wrfhydro_grid_create_from_fulldom
  public :: regrid_import_mesh_to_grid
  public :: regrid_export_grid_to_mesh

  character(len=*), parameter :: full_resolution_file = 'hydro.fullres.nc'
  character(len=:), allocatable :: hires_file
  public :: full_resolution_file
  character(len=*), parameter :: vars_out_dir = "vars_out/"
  character(len=*), parameter :: vars_in_dir = "vars_in/"
  character(len=*), parameter :: weights_dir = "weights/"

  type(ESMF_RouteHandle) :: route_handle

  character(len=ESMF_MAXSTR), parameter :: file = __FILE__
  character(len=24), parameter :: filename = "WRFHydro_NUOPC_Gluecode"
  character(len=27), parameter :: modname = "WRFHydro_NUOPC_Gluecode.F90"

  ! Debugging Parameter
  logical, parameter :: debug = .false.
  logical :: show_import_once = .true.
  logical :: show_export_once = .true.
  logical :: print_import_statelog_once = .true.
  logical :: print_export_statelog_once = .true.
  logical :: debug_importexport_vars = .TRUE.
  character(len=:), allocatable :: fname
  character(len=32) :: scount

  ! Parameters
  character(len=ESMF_MAXSTR) :: indir = 'WRFHYDRO_FORCING'
  integer                    :: num_nests = UNINITIALIZED
  integer                    :: num_tiles
  integer                    :: nx_global(1)
  integer                    :: ny_global(1)
  integer                    :: x_start
  integer                    :: x_end
  integer                    :: y_start
  integer                    :: y_end
  integer                    :: nx_local
  integer                    :: ny_local
  integer                    :: sf_surface_physics = UNINITIALIZED

  ! added to consider the adaptive time step from driver.
  real                  :: dt0 = UNINITIALIZED
  real                  :: dtrt_ter0 = UNINITIALIZED
  real                  :: dtrt_ch0 = UNINITIALIZED
  integer               :: dt_factor0 = UNINITIALIZED
  integer               :: dt_factor = UNINITIALIZED
  ! added for check soil moisture and soiltype
  integer               :: checkSOIL_flag = UNINITIALIZED
  ! added to track the driver clock
  character(len=19)     :: startTimeStr = "0000-00-00_00:00:00"

  type(ESMF_DistGrid)   :: wrfhydro_distgrid ! One DistGrid created with ConfigFile dimensions
  logical :: wrfhydro_distgrid_initialized = .false.
  character(len=512) :: logMsg

  logical :: wrfhydro_open_mesh_initialized = .false.
  type(ESMF_Mesh) :: saved_wrfhydro_mesh
  logical :: io_rank = .false.

  interface read_mesh_var_and_regrid
     module procedure read_mesh_var_and_regrid_i
     module procedure read_mesh_var_and_regrid_r
  end interface read_mesh_var_and_regrid


  !-----------------------------------------------------------------------------
  ! Model Glue Code
  !-----------------------------------------------------------------------------
contains

#undef METHOD
#define METHOD "wrfhydro_nuopc_ini"
  subroutine wrfhydro_nuopc_ini(did,vm,clock,forcingDir,rc)
    integer, intent(in)                     :: did
    type(ESMF_VM),intent(in)                :: vm
    type(ESMF_Clock),intent(in)             :: clock
    character(len=*)                        :: forcingDir
    integer, intent(out)                    :: rc

    ! local variables
    integer                     :: localPet
    integer                     :: stat
    integer, allocatable        :: deBlockList(:,:,:)
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    integer                     :: i
    type(ESMF_Time)             :: startTime
    type(ESMF_TimeInterval)     :: timeStep
    real(ESMF_KIND_R8)          :: dt
    character(ESMF_MAXSTR)      :: logMsg

    if (debug) then
       call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
    end if

    rc = ESMF_SUCCESS

    ! Set mpiCommunicator for WRFHYDRO
    call ESMF_VMGet(vm, localPet=localPet, &
      mpiCommunicator=HYDRO_COMM_WORLD, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    ! Set focing directory
    indir=forcingDir

    ! Get the models timestep
    call ESMF_ClockGet(clock,timestep=timestep,startTime=startTime,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_TimeIntervalGet(timestep,s_r8=dt,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call WRFHYDRO_TimeToString(startTime,timestr=startTimeStr,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    call orchestrator%init()


    ! Set default namelist values
    read (startTimeStr(1:4),"(I4)")   nlst(did)%START_YEAR
    read (startTimeStr(6:7),"(I2)")   nlst(did)%START_MONTH
    read (startTimeStr(9:10),"(I2)")  nlst(did)%START_DAY
    read (startTimeStr(12:13),"(I2)") nlst(did)%START_HOUR
    read (startTimeStr(15:16),"(I2)") nlst(did)%START_MIN
    nlst(did)%startdate(1:19) = startTimeStr(1:19)
    nlst(did)%olddate(1:19)   = startTimeStr(1:19)
    nlst(did)%dt = dt
    cpl_outdate = startTimeStr(1:19)
    nlst(did)%nsoil=4
    allocate(nlst(did)%zsoil8(4),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=METHOD//': Allocation of model soil depths memory failed.', &
      file=FILENAME, rcToReturn=rc)) return ! bail out
    nlst(did)%zsoil8(1:4)=(/-0.1,-0.4,-1.0,-2.0/)
    nlst(did)%geo_static_flnm = full_resolution_file
    ! nlst(did)%geo_static_flnm = "fixed/frontrange.init.fixed.nc"
    ! geo_static_flnm is read in init_namelist_rt_field call

    nlst(did)%geo_finegrid_flnm = read_hires_filename_from_namelist()
    nlst(did)%sys_cpl = 2
    nlst(did)%sys_cpl = 5 ! added to couple to MPAS
    print *, "---FIX THIS TO NAMELIST OPTION OR SOMETHING ELSE---"

    nlst(did)%IGRID = did
    write(nlst(did)%hgrid,'(I1)') did

    if(nlst(did)%dt .le. 0) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg=METHOD//": Timestep less than 1 is not supported!", &
        file=FILENAME,rcToReturn=rc)
      return  ! bail out
    endif

     ! Read information from hydro.namelist config file
    call init_namelist_rt_field(did)
    nlst(did)%geo_static_flnm = full_resolution_file

    if (debug) then
       call WRFHYDRO_nlstLog(did,MODNAME,rc=rc)
       if(ESMF_STDERRORCHECK(rc)) return ! bail out
    end if

    if(nlst(did)%nsoil .gt. 4) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg=METHOD//": Maximum soil levels supported is 4.", &
        file=FILENAME,rcToReturn=rc)
      return  ! bail out
    endif

    print *, "FIX: getting dimension from ", trim(nlst(did)%geo_static_flnm)
    ! gets - west_east, south_north,
    call get_file_dimension(fileName=nlst(did)%geo_static_flnm, &
         ix=nx_global(1), jx=ny_global(1)) !, &
         ! x_dim_var="x", y_dim_var="y")
    if (debug) print *, "nx_global(1),ny_global(1) =", nx_global(1),ny_global(1)
    call MPP_LAND_INIT(nx_global(1), ny_global(1))


    if (debug) then
       write (logMsg,"(A,2(I0,A))") MODNAME//": Global Dimensions = (", &
            nx_global(1),",",ny_global(1),")"
       call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    end if

    call log_map2d()

    call ESMF_VMBroadcast(vm, nx_global, count=1, rootPet=IO_id, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_VMBroadcast(vm, ny_global, count=1, rootPet=IO_id, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    rt_domain(did)%ix = nx_global(1)
    rt_domain(did)%jx = ny_global(1)

    if (debug) print *, "DEBUGGING: rt_domain(did)%ix, rt_domain(did)%jx", &
         rt_domain(did)%ix, rt_domain(did)%jx
    if (debug) print*, "DEBUGGING: nlst(did)%AGGFACTRT =", nlst(did)%AGGFACTRT
    ! stop "DEBUGGING GLUE"


    call MPP_LAND_PAR_INI(1,rt_domain(did)%ix,rt_domain(did)%jx,&
         nlst(did)%AGGFACTRT)

    call ESMF_VMBroadcast(vm, startx, count=numprocs, rootPet=IO_id, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_VMBroadcast(vm, starty, count=numprocs, rootPet=IO_id, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_VMBroadcast(vm, local_nx_size, count=numprocs, rootPet=IO_id, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_VMBroadcast(vm, local_ny_size, count=numprocs, rootPet=IO_id, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    allocate(deBlockList(2,2,numprocs))
    do i = 1, numprocs
      deBlockList(:,1,i) = (/startx(i),starty(i)/)
      deBlockList(:,2,i) = (/startx(i)+local_nx_size(i)-1, &
                             starty(i)+local_ny_size(i)-1/)
!      write (logMsg,"(A,I0,A,4(I0,A))") MODNAME//": deBlockList ", i, " = (", &
!        deBlockList(1,1,i),":",deBlockList(1,2,i),",", &
!        deBlockList(2,1,i),":",deBlockList(2,2,i),")"
!      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    enddo

!    allocate(connectionList(1),stat=stat)
!    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
!      msg=METHOD//': Allocation of connection list memory failed.', &
!      file=FILENAME, rcToReturn=rc)) return ! bail out
!    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
!      tileIndexB=1, positionVector=(/nx_global(1), 0/), rc=rc)
!    if (ESMF_STDERRORCHECK(rc)) return  ! bail out


    ! Create DistGrid based on WRFHDYRO Config NX,NY
    wrfhydro_distgrid = ESMF_DistGridCreate( &
      minIndex=(/1,1/), maxIndex=(/nx_global(1),ny_global(1)/), &
!     indexflag = ESMF_INDEX_DELOCAL, &
     deBlockList=deBlockList, &
!     deLabelList=deLabelList, &
!     delayout=delayout, &
!     connectionList=connectionList, &
      rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    wrfhydro_distgrid_initialized = .true.
    deallocate(deBlockList)

!   deallocate(connectionList,stat=stat)
!   if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
!     msg=METHOD//': Deallocation of connection list memory failed.', &
!     file=FILENAME,rcToReturn=rc)) return ! bail out

    ! Get the Local Decomp Incides
    call set_local_indices(rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": Enter CPL_LAND_INIT", ESMF_LOGMSG_INFO)
#endif
    ! Initialize the internal Land <-> Hydro Coupling
    call CPL_LAND_INIT(x_start, x_end, y_start, y_end)
#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": Exit CPL_LAND_INIT", ESMF_LOGMSG_INFO)
#endif

    ! Routing timestep set in HYDRO_ini
#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": Enter HYDRO_ini", ESMF_LOGMSG_INFO)
#endif

    if(sf_surface_physics .eq. 5) then
      ! clm4
      ! Use wrfinput vegetation type and soil type
      call HYDRO_ini(ntime=1,did=did,ix0=1,jx0=1)
    else
      ! Use wrfinput vegetation type and soil type
      call HYDRO_ini(ntime=1,did=did,ix0=nx_local,jx0=ny_local)
    endif
#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": Exit HYDRO_ini", ESMF_LOGMSG_INFO)
    call WRFHYDRO_domainLog(did,MODNAME,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
#endif

    ! Override the clock configuration in hyro.namelist
    read (startTimeStr(1:4),"(I4)")   nlst(did)%START_YEAR
    read (startTimeStr(6:7),"(I2)")   nlst(did)%START_MONTH
    read (startTimeStr(9:10),"(I2)")  nlst(did)%START_DAY
    read (startTimeStr(12:13),"(I2)") nlst(did)%START_HOUR
    read (startTimeStr(15:16),"(I2)") nlst(did)%START_MIN
    nlst(did)%startdate(1:19) = startTimeStr(1:19)
    nlst(did)%olddate(1:19)   = startTimeStr(1:19)
    nlst(did)%dt = dt
    nlst(did)%nsoil=4
    cpl_outdate = startTimeStr(1:19)

    if(nlst(did)%dt .le. 0) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg=METHOD//": Timestep less than 1 is not supported!", &
        file=FILENAME,rcToReturn=rc)
      return  ! bail out
    endif

    ! Adjust the routing timestep and factor
    ! At this point the coupling driver timestep is unknown
    ! and uses WRFHYDRO Config as best guess
    if(nlst(did)%dtrt_ter .ge. nlst(did)%dt) then
       nlst(did)%dtrt_ter = nlst(did)%dt
       dt_factor0 = 1
    else
       dt_factor = nlst(did)%dt/nlst(did)%dtrt_ter
       if (dt_factor*nlst(did)%dtrt_ter .lt. nlst(did)%dt) &
         nlst(did)%dtrt_ter = nlst(did)%dt/dt_factor
       dt_factor0 = dt_factor
    endif

    if(nlst(did)%dtrt_ch .ge. nlst(did)%dt) then
      nlst(did)%dtrt_ch = nlst(did)%dt
      dt_factor0 = 1
    else
      dt_factor = nlst(did)%dt/nlst(did)%dtrt_ch
      if(dt_factor*nlst(did)%dtrt_ch .lt. nlst(did)%dt) &
        nlst(did)%dtrt_ch = nlst(did)%dt/dt_factor
      dt_factor0 = dt_factor
    endif

    dt0 = nlst(did)%dt
    dtrt_ter0 = nlst(did)%dtrt_ter
    dtrt_ch0 = nlst(did)%dtrt_ch

    RT_DOMAIN(did)%initialized = .true.

    num_nests = num_nests + 1

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

    ! THIS IS ZERO :D
    ! print *, "sfchead=",&
    !      rt_domain(did)%overland%control%surface_water_head_lsm(:,:)
    ! stop "SFCHEAD AT END OF WRFHYDRO_NUOPC_INI"
  end subroutine wrfhydro_nuopc_ini

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "wrfhydro_nuopc_run"

  subroutine wrfhydro_nuopc_run(did,lsm_forcings,clock,importState,&
  exportState,rc)
    integer, intent(in)                     :: did
    logical, intent(in)                     :: lsm_forcings
    type(ESMF_Clock),intent(in)             :: clock
    type(ESMF_State),intent(inout)          :: importState
    type(ESMF_State),intent(inout)          :: exportState
    integer, intent(out)                    :: rc

    ! local variables
    type(ESMF_TimeInterval)     :: timeStep

    ! debugging
    integer :: rank, ierr

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    if(.not. RT_DOMAIN(did)%initialized) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="WRHYDRO: Model has not been initialized!", &
        file=FILENAME,rcToReturn=rc)
      return  ! bail out
    endif

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    call WRFHYDRO_ClockToString(clock,timestr=cpl_outdate,rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return
    nlst(did)%olddate(1:19) = cpl_outdate(1:19) ! Current time is the

    nlst(did)%dt = WRFHYDRO_TimeIntervalGetReal(timeInterval=timeStep,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    if(nlst(did)%dt .le. 0) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg=METHOD//": Timestep less than 1 is not supported!", &
        file=FILENAME,rcToReturn=rc)
      return  ! bail out
    endif

    if((dt_factor0*nlst(did)%dtrt_ter) .ne. nlst(did)%dt) then   ! NUOPC driver time step changed.
      call ESMF_LogWrite(METHOD//": Driver timestep changed.",ESMF_LOGMSG_INFO)
      if(dtrt_ter0 .ge. nlst(did)%dt) then
        nlst(did)%dtrt_ter = nlst(did)%dt
        dt_factor0 = 1
      else
        dt_factor = nlst(did)%dt / dtrt_ter0
        if(dt_factor*dtrt_ter0 .lt. nlst(did)%dt) &
          nlst(did)%dtrt_ter = nlst(did)%dt / dt_factor
        dt_factor0 = dt_factor
      endif
    endif

    if((dt_factor0*nlst(did)%dtrt_ch) .ne. nlst(did)%dt) then   ! NUOPD driver time step changed.
      call ESMF_LogWrite(METHOD//": Driver timestep changed.",ESMF_LOGMSG_INFO)
      if(dtrt_ch0 .ge. nlst(did)%dt) then
        nlst(did)%dtrt_ch = nlst(did)%dt
        dt_factor0 = 1
      else
        dt_factor = nlst(did)%dt / dtrt_ch0
        if(dt_factor*dtrt_ch0 .lt. nlst(did)%dt) &
          nlst(did)%dtrt_ch = nlst(did)%dt / dt_factor
        dt_factor0 = dt_factor
      endif
    endif

    if(nlst(did)%SUBRTSWCRT .eq.0  .and. &
      nlst(did)%OVRTSWCRT .eq. 0 .and. &
      nlst(did)%GWBASESWCRT .eq. 0) then
       call ESMF_LogWrite(METHOD//": SUBRTSWCRT,OVRTSWCRT,GWBASESWCRT are zero!", &
            ESMF_LOGMSG_WARNING)
      !call ESMF_LogSetError(ESMF_FAILURE, &
      !  msg=METHOD//": SUBRTSWCRT,OVRTSWCRT,GWBASESWCRT are zero!", &
      !  file=FILENAME,rcToReturn=rc)
      !return  ! bail out
    endif

    if (debug) print *, "WRF-H GLUE: lsm_forcings =", lsm_forcings
    if((.not. RT_DOMAIN(did)%initialized) .and. (nlst(did)%rst_typ .eq. 1) ) then
      call ESMF_LogWrite(METHOD//": Restart initial data from offline file.", &
        ESMF_LOGMSG_INFO)
    else
      if (.not. lsm_forcings) then
        call read_ldasout(olddate=nlst(did)%olddate(1:19), &
          hgrid=nlst(did)%hgrid, &
          indir=trim(indir), dt=nlst(did)%dt, &
          ix=rt_domain(did)%ix,jx=rt_domain(did)%jx, &
          infxsrt=rt_domain(did)%infxsrt,soldrain=rt_domain(did)%soldrain)
      endif
    endif
    if (debug) print *, "DEBUGGING: COMPLETE"

    ! print *, "sfchead=",&
    !      rt_domain(1)%overland%control%surface_water_head_lsm(25:50,25:50)
    ! sfchead 0 here

    ! Call the WRF-HYDRO run routine
    call HYDRO_exe(did=did)

    ! print *, "sfchead=",&
    !      rt_domain(1)%overland%control%surface_water_head_lsm(40:50,40:50)
    ! stop "IT IS non-0 here??, but how big?!"

    ! rt_domain(did)%smc(:,:,:) = 0.333 * rank
    ! rt_domain(did)%overland%control%surface_water_head_lsm = 0.555 * rank

    ! rt_domain(did)%overland%control%surface_water_head_lsm = 0
    ! rt_domain(did)%overland%control%surface_water_head_lsm = 0.555 * rank
    ! print *, "CHANGING SFCHEAD TO TEST"
    ! rt_domain(did)%overland%control%surface_water_head_lsm = my_id / 10
    ! print *, "CHANGING SMC TO TEST"
    ! rt_domain(did)%smc(:,:,:) = 0.111

    ! provide groundwater soil flux to WRF for fully coupled simulations (FERSCH 09/2014)
    !if(nlst(did)%GWBASESWCRT .eq. 3 ) then
      !Wei Yu: comment the following two lines. Not ready
    !yw     qsgw(x_start(1):x_end(1),y_start(1):y_end(1)) = gw2d(did)%qsgw
    !yw     config_flags%gwsoilcpl = nlst(did)%gwsoilcpl
    !end if

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine wrfhydro_nuopc_run

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "wrfhydro_nuopc_fin"

  subroutine wrfhydro_nuopc_fin(did,rc)
    ! ARGUMENTES
    integer, intent(inout)      :: did
    integer, intent(out)        :: rc

    ! LOCAL VARIABLES
    integer                     :: stat

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    ! WRF-Hydro finish routine cannot be called because it stops MPI

!    DCR - Turned off to let the model deallocate memory
!    deallocate(nlst(did)%zsoil8)
!    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
!      msg=METHOD//': Deallocation of model soil depth memory failed.', &
!      file=FILENAME,rcToReturn=rc)) return ! bail out

    RT_DOMAIN(did)%initialized = .false.

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine


  subroutine write_netcdf_full_resolution_file(dx, dy, hgt, soil_cat, lat, &
       lon, lu_index, landmask)
    integer, intent(in) :: dx, dy
    real, dimension(:,:), intent(in) :: hgt, lat, lon
    integer, dimension(:,:), intent(in) :: soil_cat, lu_index, landmask

    ! Constants
    character(len=*), parameter :: times_value = '0000-00-00_00:00:00'
    integer, parameter :: DateStrLen = len(times_value)

    ! Locals
    integer :: stat, ncid
    integer :: dim_time, dim_we, dim_sn, dim_dsl
    integer :: dim_soil_cat, dim_land_cat
    integer :: var_times, var_hgt, var_soil_cat, var_lat, var_lon, &
         var_lu_index, var_landmask

    ! Local parameters
    character(len=*), parameter :: timestr = "0000-00-00_00:00:00"
    integer, parameter :: soil_cat_n = 16
    integer, parameter :: land_cat_n = 24

    if (io_rank) print *, "=== entering write_netcdf_full_resolution_file ==="
    ! Begin work
    ! Create/clobber file
    call check_nf(nf90_create(full_resolution_file, NF90_CLOBBER, ncid))

    ! Define Dimensions
    call check_nf(nf90_def_dim(ncid, 'Time', NF90_UNLIMITED, dim_time))
    call check_nf(nf90_def_dim(ncid, 'west_east', dx, dim_we))
    call check_nf(nf90_def_dim(ncid, 'south_north', dy, dim_sn))
    call check_nf(nf90_def_dim(ncid, 'soil_cat', soil_cat_n, dim_soil_cat))
    call check_nf(nf90_def_dim(ncid, 'land_cat', land_cat_n, dim_land_cat))
    call check_nf(nf90_def_dim(ncid, 'DateStrLen', DateStrLen, dim_dsl))

    ! Define Variables
    call check_nf(nf90_def_var(ncid, 'Times', NF90_CHAR, [dim_dsl, dim_time], var_times))
    call check_nf(nf90_def_var(ncid, 'HGT_M', NF90_FLOAT, [dim_we, dim_sn, dim_time], var_hgt))
    call check_nf(nf90_def_var(ncid, 'SCT_DOM', NF90_INT, [dim_we, dim_sn, dim_time], var_soil_cat))
    call check_nf(nf90_def_var(ncid, 'XLAT_M', NF90_FLOAT, [dim_we, dim_sn, dim_time], var_lat))
    call check_nf(nf90_def_var(ncid, 'XLONG_M', NF90_FLOAT, [dim_we, dim_sn, dim_time], var_lon))
    ! print *, "TOOK OUT VEGTYP, IS THIS CORRECT?"
    ! call check_nf(nf90_def_var(ncid, 'VEGTYP', NF90_INT, [dim_we, dim_sn, dim_time], var_veg))
    call check_nf(nf90_def_var(ncid, 'LU_INDEX', NF90_INT, [dim_we, dim_sn, dim_time], var_lu_index))
    call check_nf(nf90_def_var(ncid, 'LANDMASK', NF90_INT, [dim_we, dim_sn, dim_time], var_landmask))

    ! Attributes
    call check_nf(nf90_put_att(ncid, NF90_GLOBAL, 'ISWATER', 16))
    call check_nf(nf90_put_att(ncid, NF90_GLOBAL, 'ISLAKE', -1))
    call check_nf(nf90_put_att(ncid, NF90_GLOBAL, 'ISICE', 24))
    call check_nf(nf90_put_att(ncid, NF90_GLOBAL, 'ISURBAN', 1))
    call check_nf(nf90_put_att(ncid, NF90_GLOBAL, 'ISOILWATER', 14))

    ! call check_nf(nf90_put_att(ncid, var_hgt, 'description', 'topography height'))

    ! End define mode before writing data
    call check_nf(nf90_enddef(ncid))

    ! Write Variables
    call check_nf(nf90_put_var(ncid, var_times, timestr, start=[1,1], count=[DateStrLen,1]))
    call check_nf(nf90_put_var(ncid, var_hgt, hgt, start=[1,1,1], count=[dx,dy,1]))
    call check_nf(nf90_put_var(ncid, var_soil_cat, soil_cat, start=[1,1,1], count=[dx,dy,1]))
    call check_nf(nf90_put_var(ncid, var_lat, lat, start=[1,1,1], count=[dx,dy,1]))
    call check_nf(nf90_put_var(ncid, var_lon, lon, start=[1,1,1], count=[dx,dy,1]))
    ! call check_nf(nf90_put_var(ncid, var_veg, veg, start=[1,1,1], count=[dx,dy,1]))
    call check_nf(nf90_put_var(ncid, var_lu_index, lu_index, start=[1,1,1], count=[dx,dy,1]))
    call check_nf(nf90_put_var(ncid, var_landmask, landmask, start=[1,1,1], count=[dx,dy,1]))
    ! call check_nf(nf90_put_var(ncid, var_, , start=[1,1,1], count=[dx,dy,1]))



    call check_nf(nf90_close(ncid))
    print *, "=== wrote to ", full_resolution_file, " ==="
    print *, "=== exiting write_netcdf_full_resolution_file ==="
  end subroutine write_netcdf_full_resolution_file

  subroutine read_mesh_var_and_regrid_r(var_name, outvar, &
       ncid, nCells, &
       regrid_handle, srcPtr, dstPtr, f_src, f_dst)
    character(len=*), intent(in) :: var_name
    real, allocatable, intent(inout) :: outvar(:,:)
    integer, intent(in) :: ncid, nCells
    real(ESMF_KIND_R8), pointer, intent(inout) :: srcPtr(:), dstPtr(:,:)
    type(ESMF_Field) :: f_src, f_dst
    type(ESMF_RouteHandle), intent(inout) :: regrid_handle

    ! local
    real, allocatable :: var(:)
    integer :: rc, stat, varid
    print *, "Regridding ", trim(var_name)

    ! read in NetCDF Vars
    stat = nf90_inq_varid(ncid, var_name, varid)
    call check_nf(stat)
    allocate(var(nCells))
    stat = nf90_get_var(ncid, varid, var)
    call check_nf(stat)

    ! print *, "ncells =", ncells
    ! print *, "var_r(1:10) =", var(1:10)

    if (size(srcPtr) /= size(var)) then
       rc = ESMF_RC_ARG_SIZE
       stop 'var_r length does not match local mesh ELEMENT count.'
    end if
    ! if (size(dstPtr,1) /= size(outvar,1) .or. &
    !     size(dstPtr,2) /= size(outvar,2)) then
    !    rc = ESMF_RC_ARG_SIZE
    !     stop 'outvar shape does not match local grid CENTER size.'
    ! end if

    srcPtr = real(var, kind=ESMF_KIND_R8)
    dstPtr = 0.0_ESMF_KIND_R8

    call ESMF_FieldRegrid(srcField=f_src, dstField=f_dst, &
         routehandle=regrid_handle, &
         rc=rc)
    call check(rc, __LINE__, file)

    outvar = dstPtr

    print *, "  outvar_r(1:2,1:2) =", outvar(1:2,1:2)
    print *, "  outvar_r shape =", shape(outvar)
  end subroutine read_mesh_var_and_regrid_r


  subroutine read_mesh_var_and_regrid_i(var_name, outvar, &
       ncid, nCells, &
       regrid_handle, srcPtr, dstPtr, f_src, f_dst)
    character(len=*), intent(in) :: var_name
    integer, allocatable, intent(inout) :: outvar(:,:)
    integer, intent(in) :: ncid, nCells
    real(ESMF_KIND_R8), pointer, intent(inout) :: srcPtr(:), dstPtr(:,:)
    type(ESMF_Field) :: f_src, f_dst
    type(ESMF_RouteHandle), intent(inout) :: regrid_handle

    ! local
    integer, allocatable :: var_i(:)
    integer :: rc, stat, varid
    print *, "Regridding ", trim(var_name)

    ! read in NetCDF Vars
    stat = nf90_inq_varid(ncid, var_name, varid)
    call check_nf(stat)
    allocate(var_i(nCells))
    stat = nf90_get_var(ncid, varid, var_i)
    call check_nf(stat)

    if (size(srcPtr) /= size(var_i)) then
       rc = ESMF_RC_ARG_SIZE
       print *, "size(srcPtr) =", size(srcPtr), "size(var_i)) =", size(var_i)
       error stop 'var_i length does not match local mesh ELEMENT count.'
    end if
    ! if (size(dstPtr,1) /= size(outvar,1) .or. &
    !     size(dstPtr,2) /= size(outvar,2)) then
    !    rc = ESMF_RC_ARG_SIZE
    !     stop 'outvar shape does not match local grid CENTER size.'
    ! end if

    srcPtr = real(var_i, kind=ESMF_KIND_R8)
    dstPtr = 0.0_ESMF_KIND_R8

    call ESMF_FieldRegrid(srcField=f_src, dstField=f_dst, &
         routehandle=regrid_handle, &
         rc=rc)
    call check(rc, __LINE__, file)

    ! print*, trim(var_name), "dstPtr =", dstPtr
    outvar = nint(dstPtr)

    print *, "  outvar_i(1:2,1:2) =", outvar(1:2,1:2)
    print *, "  outvar_i shape =", shape(outvar)
  end subroutine read_mesh_var_and_regrid_i

  subroutine wrfhydro_write_full_resolution_file(wrfhydro_grid, wrfhydro_mesh, &
       regrid_handle_bl, regrid_handle_nn_stod, regrid_handle_nn_dtos)
    use netcdf
    type(ESMF_Grid), intent(in) :: wrfhydro_grid
    type(ESMF_Mesh), intent(in) :: wrfhydro_mesh
    type(ESMF_RouteHandle), intent(inout) :: regrid_handle_bl
    type(ESMF_RouteHandle), intent(inout) :: regrid_handle_nn_stod
    type(ESMF_RouteHandle), intent(inout) :: regrid_handle_nn_dtos
    ! type(ESMF_Mesh)            :: wrfhydro_mesh
    character(:), allocatable :: mpas_file

    ! NetCDF Variables
    integer :: stat, ncid, dimid, nCells, varid
    integer :: ncid_out
    integer, allocatable :: var_i(:)
    real, allocatable :: var_r(:)

    ! ESMF Variables
    type(ESMF_Field) :: f_src, f_dst
    real(ESMF_KIND_R8), pointer :: srcPtr(:) => null()
    real(ESMF_KIND_R8), pointer :: dstPtr(:,:) => null()
    integer :: rc

    ! vars to write
    real, dimension(:,:), allocatable :: hgt, lat, lon
    integer, dimension(:,:), allocatable :: soil_cat, lu_index, landmask

    ! locals
    integer :: nx, ny
    ! type(ESMF_RouteHandle) :: regrid_handle_tmp

    print *, "=== entering wrfhydro_write_full_resolution_file ==="
    ! initialize values
    ! regrid_handle_tmp = regrid_handle

    ! setup MPAS NetCDF variables
    mpas_file = "frontrange.static.nc"
    print *, "TODO: mpas_file for geo is hardcoded to: ", mpas_file

    stat = nf90_open(trim(mpas_file), NF90_NOWRITE, ncid)
    call check_nf(stat)

    stat = nf90_inq_dimid(ncid, 'nCells', dimid)
    call check_nf(stat)

    stat = nf90_inquire_dimension(ncid, dimid, len=nCells)
    call check_nf(stat)

    ! setup MPAS regridding
    !   regrid MPAS mesh to WRF-Hydro hi-res grid
    !   Create source and destination field on grid (CENTER)
    f_src = ESMF_FieldCreate(mesh=wrfhydro_mesh, &
         typekind=ESMF_TYPEKIND_R8, &
         meshloc=ESMF_MESHLOC_ELEMENT, &
         name='var_src', rc=rc)
    call check(rc, __LINE__, file)
    f_dst = ESMF_FieldCreate(grid=wrfhydro_grid, &
         typekind=ESMF_TYPEKIND_R8, &
         indexflag=ESMF_INDEX_DELOCAL, &
         name='var_dst', rc=rc)
    call check(rc, __LINE__, file)

    call ESMF_FieldGet(f_src, farrayPtr=srcPtr, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_FieldGet(f_dst, farrayPtr=dstPtr, rc=rc)
    call check(rc, __LINE__, file)

    ! print *, "prencells =", ncells
    ! print *, "prelat(1:2,1:2) =", lat(1:2,1:2)
    ! print *, "prelon(1:2,1:2) =", lon(1:2,1:2)
    ! print *, "prehgt(1:2,1:2) =", hgt(1:2,1:2)
    ! print *, "presoil_cat(1:2,1:2) =", soil_cat(1:2,1:2)
    ! print *, "preveg(1:2,1:2) =", veg(1:2,1:2)
    ! print *, "prelu_index(1:2,1:2) =", lu_index(1:2,1:2)
    ! print *, "prelandmask(1:2,1:2) =", landmask(1:2,1:2)
    ! print *, "=================="


    ! nearest neighbor regridding
    call read_mesh_var_and_regrid('isltyp', soil_cat, ncid, nCells, &
         regrid_handle_nn_stod, srcPtr, dstPtr, f_src, f_dst)
    call read_mesh_var_and_regrid('ivgtyp', lu_index, ncid, nCells, &
         regrid_handle_nn_stod, srcPtr, dstPtr, f_src, f_dst)
    call read_mesh_var_and_regrid('landmask', landmask, ncid, nCells, &
         regrid_handle_nn_stod, srcPtr, dstPtr, f_src, f_dst)
    ! bilinear regridding
    call read_mesh_var_and_regrid('ter', hgt, ncid, nCells, &
         regrid_handle_bl, srcPtr, dstPtr, f_src, f_dst)
    call read_mesh_var_and_regrid('latCell', lat, ncid, nCells, &
         regrid_handle_bl, srcPtr, dstPtr, f_src, f_dst)
    call read_mesh_var_and_regrid('lonCell', lon, ncid, nCells, &
         regrid_handle_bl, srcPtr, dstPtr, f_src, f_dst)

    stat = nf90_close(ncid)
    call check_nf(stat)
    call ESMF_FieldDestroy(f_src, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_FieldDestroy(f_dst, rc=rc)
    call check(rc, __LINE__, file)


    nx = size(hgt, dim=1)
    ny = size(hgt, dim=2)
    call write_netcdf_full_resolution_file(nx, ny, hgt, soil_cat, lat, lon, &
         lu_index, landmask)

    print *, "=== exiting wrfhydro_write_full_resolution_file ==="
  end subroutine wrfhydro_write_full_resolution_file

  function wrfhydro_open_mesh(vm, rc) result(mesh)
    type(ESMF_VM), intent(inout) :: vm
    integer, intent(out) :: rc
    type(ESMF_Mesh) :: mesh
    type(ESMF_DistGrid) :: distgrid
    integer, allocatable :: gindex(:)
    character(:), allocatable :: mpas_grid_file, scrip_mesh_file
    integer :: unit, iostat, irank, localCount, tmp, ierr
    integer :: idx, inode
    logical :: exists

    rc = ESMF_SUCCESS

    if (wrfhydro_open_mesh_initialized) then
       mesh = saved_wrfhydro_mesh
       return
    end if


    mpas_grid_file = get_mpas_grid_filename()
    scrip_mesh_file = mpas_to_scrip_filename(mpas_grid_file)

    print *, "mpas_grid_file=", mpas_grid_file
    print *, "scrip_mesh_file=", scrip_mesh_file

    inquire(file=trim(scrip_mesh_file), exist=exists)
    if (.not. exists) then
       call mpas_to_scrip_mesh(mpas_grid_file, scrip_mesh_file)
    end if

    distgrid  = get_mpas_dist_grid(vm, rc)
    mesh = ESMF_MeshCreate(filename=scrip_mesh_file, &
         elementDistgrid=distgrid, &
         ! convertToDual=.false., & ! didn't do anything
         fileformat=ESMF_FILEFORMAT_SCRIP, rc=rc)
    call check(rc, __LINE__, file)

    ! save initialized mesh for future calls
    saved_wrfhydro_mesh = mesh
    wrfhydro_open_mesh_initialized = .true.
  end function wrfhydro_open_mesh

  function wrfhydro_regrid_mesh(input_grid, mesh, regrid_method, &
       did, importState, vm, rc) result(handle)
    type(ESMF_Grid), intent(in) :: input_grid
    type(ESMF_Mesh), intent(in) :: mesh
    type(ESMF_RegridMethod_Flag), intent(in) :: regrid_method
    integer, intent(in) :: did
    type(ESMF_State), intent(inout) :: importState
    type(ESMF_VM), intent(inout) :: vm
    integer, intent(out) :: rc
    type(ESMF_RouteHandle) :: handle
    type(ESMF_Grid) :: grid

    type(ESMF_Field) :: import_field, new_field
    character(:), allocatable :: file

    character(:), allocatable :: mpas_grid_file, scrip_mesh_file, hires_scrip_file

    ! testing variables
    logical :: realizeImport, connected
    character(:), allocatable :: st_name



    ! error stop "WRFHYDRO_REGRID_MESH IS ENTERED?"
    file = __FILE__
    rc = ESMF_SUCCESS
    print *, "enter: wrfhydro_regrid_mesh"
    call ESMF_LogWrite("WRFH: enter wrfhydro_regrid_mesh", &
         ESMF_LOGMSG_INFO, rc=rc)

    if (regrid_method == ESMF_REGRIDMETHOD_BILINEAR) then
       st_name = "bilinear"
    else if (regrid_method == ESMF_REGRIDMETHOD_CONSERVE) then
       st_name = "conserve"
    else if (regrid_method == ESMF_REGRIDMETHOD_NEAREST_STOD) then
       st_name = "nn_stod"
    else if (regrid_method == ESMF_REGRIDMETHOD_NEAREST_DTOS) then
       st_name = "nn_dtos"
    else
       error stop "regridding method not supported"
    end if
    ! check for import field

    if (debug) then
       print *, "========", trim(st_name), "========"
       call ESMF_StateLog(importState, logMsgFlag=ESMF_LOGMSG_INFO, rc=rc)
       if (ESMF_STDERRORCHECK(rc)) return
    end if

    import_field = ESMF_FieldCreate(name=st_name, &
         mesh=mesh, typekind=ESMF_TYPEKIND_R8, &
         meshloc=ESMF_MESHLOC_ELEMENT, & ! ESMF_MESHLOC_ELEMENT for cell-center vars
         rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return
    new_field = ESMF_FieldCreate(name=st_name, &
         staggerloc=ESMF_STAGGERLOC_CENTER,         &
         grid=input_grid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return

    ! THIS DOESN'T WORK BECAUSE IT IS A REDISTRIBUTION, WHEN AN INTERPOLATION
    !   IS NEEDED
    ! call ESMF_FieldRedistStore(srcField=import_field, dstField=new_field, &
    !      routehandle=route_handle, rc=rc)

    ! Build the interpolation operator once
    if (debug) print *, "Initialized route_handle"
    call ESMF_LogWrite("Initialized route_handle", &
         ESMF_LOGMSG_INFO, rc=rc)

    ! stop "SHOULD THIS BE REACHED?"
    ! generate regrid weights file and read in to handle
    if (.not. allocated(hires_file)) &
         error stop "hires_file variable not defined"
    call create_dir_if_needed(weights_dir)

    mpas_grid_file = get_mpas_grid_filename()
    scrip_mesh_file = mpas_to_scrip_filename(mpas_grid_file)
    hires_scrip_file = fulldom_to_scrip_filename(hires_file)
    print *, "mpas_grid_file =", trim(mpas_grid_file)
    print *, "scrip_grid_file =", trim(scrip_mesh_file)
    print *, "hires_scrip_file =", trim(hires_scrip_file)
    call ESMF_RegridWeightGen(&
         srcFile=scrip_mesh_file, &
         dstFile=hires_scrip_file, &
         weightFile=weights_dir//'setup_'//st_name//'.nc', &
         regridmethod=regrid_method, &
         srcFileType=ESMF_FILEFORMAT_SCRIP, &
         dstFileType=ESMF_FILEFORMAT_SCRIP, &
         srcRegionalFlag=.true., &
         dstRegionalFlag=.true., &
         normType=ESMF_NORMTYPE_FRACAREA, &
         rc=rc)
    call check(rc, __LINE__, file)
    ! Precompute Field sparse matrix multiplication with local factors
    call ESMF_FieldSMMStore(srcField=import_field, dstField=new_field, &
         filename=weights_dir//'setup_'//st_name//'.nc', &
         routehandle=route_handle, rc=rc)
    call check(rc, __LINE__, file)

    ! generate regrid weights in to handle
    ! call ESMF_FieldRegridStore(srcfield=import_field, dstfield=new_field, &
    !      routehandle=route_handle, rc=rc, &
    !      regridmethod=ESMF_REGRIDMETHOD_BILINEAR, & ! or CONSERVE / PATCH
    !      lineType=ESMF_LINETYPE_GREAT_CIRCLE, &
    !      unmappedaction=ESMF_UNMAPPEDACTION_IGNORE)
    ! write handle to binary file
    ! call ESMF_RouteHandleWrite(route_handle, fileName="mpas_hydro.RH", rc=rc)
    ! if(ESMF_STDERRORCHECK(rc)) return
    if(ESMF_STDERRORCHECK(rc)) return


    ! each step (after src has values):
    print *, "Regrid field, applying weights"
    call ESMF_LogWrite("Regrid field, applying weights", &
         ESMF_LOGMSG_INFO, rc=rc)
    ! Compute a regridding operation
    call ESMF_FieldRegrid(import_field, new_field, route_handle, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return



    ! if (ESMF_STDERRORCHECK(rc)) return
    ! call NUOPC_Realize(importState, field=field_import, rc=rc)
    ! if (ESMF_STDERRORCHECK(rc)) return
    ! connected = NUOPC_IsConnected(importState, &
    !      fieldName=trim(st_name),rc=rc)


    ! ! get import field
    ! ! field_import=field_create(fld_name=fieldList(n)%st_name, &
    ! !       mesh=mesh, did=did, memflg=memr_import, rc=rc)

    ! ! field = ESMF_FieldCreate(name="sss", mesh=meshIn, &
    ! !   typekind=ESMF_TYPEKIND_R8, rc=rc)
    ! ! call NUOPC_Realize(importState, field=field, rc=rc)

    ! memory cleanup, only used to create route handle
    call ESMF_FieldDestroy(import_field, rc=rc)
    call check(rc, __LINE__, file)

    call ESMF_FieldDestroy(new_field, rc=rc)
    call check(rc, __LINE__, file)


    handle = route_handle ! return variable
    call ESMF_LogWrite("WRFH: exit wrfhydro_regrid_mesh", &
         ESMF_LOGMSG_INFO, rc=rc)
    print *, "WRFH: exit wrfhydro_regrid_mesh"
    ! error stop "where is this being called?"
  end function wrfhydro_regrid_mesh

  subroutine regrid_import_mesh_to_grid(grid, mesh, state, &
       did, memflg)
    type(ESMF_Grid), intent(in) :: grid
    type(ESMF_Mesh), intent(in) :: mesh
    type(ESMF_State),intent(in) :: state
    integer, intent(in) :: did
    type(memory_flag), intent(in) :: memflg
    type(ESMF_Field) :: meshField !, gridField
    integer :: n, rc, itemCount
    character(len=64), allocatable :: itemNameList(:)
    logical :: imported
    character(:), allocatable :: mpas_grid_file, scrip_mesh_file, hires_scrip_file

    ! debugging
    integer :: unit, i
    logical :: exists
    integer :: rank, ierr
    character(len=3) :: rank_s

    if (print_import_statelog_once) then
       if (io_rank) print *, 'WRFH: enter regrid_import_mesh_to_grid'
       call ESMF_LogWrite('WRFH: enter regrid_import_mesh_to_grid', ESMF_LOGMSG_INFO)
       call ESMF_LogWrite('--- WRFH: IMPORT STATE DEBUG: ', ESMF_LOGMSG_INFO)
       call ESMF_StateLog(state, rc=rc)
       call ESMF_LogWrite('--- WRFH: POST STATE PRINT: ', ESMF_LOGMSG_INFO)
       print_import_statelog_once = .false.
    end if


    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)
    allocate(itemNameList(itemCount), stat=rc)

    call ESMF_StateGet(state,itemNameList=itemNameList, rc=rc)
    call check(rc, __LINE__, file)
    if (debug) then
       do i=1,itemCount
          print *, "WRFH: import itemnamelist: ", itemNameList(i)
          call ESMF_LogWrite("WRFH: import itemnamelist: "// itemNameList(i), ESMF_LOGMSG_INFO, rc=rc)
       end do
    end if

    do n=lbound(cap_fld_list,1),ubound(cap_fld_list,1)
       if (cap_fld_list(n)%ad_import) then
          imported = NUOPC_IsConnected(state, &
            fieldName=trim(cap_fld_list(n)%st_name), rc=rc)
          call check(rc, __LINE__, file)

          if (debug) print *, "WRFH: imported = ", imported,&
               ", name: ", cap_fld_list(n)%st_name
          if (debug) call ESMF_LogWrite("WRFH: imported "//cap_fld_list(n)%st_name, &
               ESMF_LOGMSG_INFO, rc=rc)

          call ESMF_StateGet(state, itemName=trim(cap_fld_list(n)%st_name), &
               field=meshField, rc=rc)
          call check(rc, __LINE__, file)

          ! gridField = field_create(cap_fld_list(n)%st_name, grid, did, &
          !      memflg, rc)
          call ESMF_GridValidate(grid, rc=rc)
          call check(rc, __LINE__, file)


          ! new way to handle gridField
          if (cap_fld_list(n)%import_field_init .eqv. .false.) then
             cap_fld_list(n)%import_field = field_create(cap_fld_list(n)%st_name, grid, did, &
                  memflg, rc)
             cap_fld_list(n)%import_field_init = .true.
          end if

          ! CREAT LOCAL COPY, NOT CONNECTED TO WRF-H MODEL MEMORY
          ! gridField = ESMF_FieldCreate(grid=grid, &
          !      typekind=ESMF_TYPEKIND_R8, &
          !      indexflag=ESMF_INDEX_DELOCAL, &
          !      name=cap_fld_list(n)%st_name, rc=rc)
          call check(rc, __LINE__, file)


          call check(rc, __LINE__, file)

          ! --- Debugging: write import mesh before regrid
          ! if (debug) then
          if (show_import_once) then
                 call create_dir_if_needed(vars_in_dir)
             call ESMF_FieldWrite(cap_fld_list(n)%import_field, &
                  vars_in_dir//trim(cap_fld_list(n)%st_name)//&
                  "_grid_pre.nc", &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call ESMF_FieldWrite(meshField, &
                  vars_in_dir//trim(cap_fld_list(n)%st_name)//&
                  "_mesh_pre.nc", &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call check(rc, __LINE__, file)
          end if

          if (cap_fld_list(n)%import_handle_init .eqv. .false.) then
             ! call ESMF_FieldRegridStore(srcField=meshField, &
             call check(rc, __LINE__, file)

             if (.not. allocated(hires_file)) &
                  hires_file = read_hires_filename_from_namelist()
             hires_scrip_file = fulldom_to_scrip_filename(hires_file)
             mpas_grid_file = get_mpas_grid_filename()
             scrip_mesh_file = mpas_to_scrip_filename(mpas_grid_file)

             call create_dir_if_needed(weights_dir)
             ! center, dims and imask
             call ESMF_RegridWeightGen(&
                  srcFile=scrip_mesh_file, &
                  dstFile=hires_scrip_file, &
                  weightFile=weights_dir//trim(cap_fld_list(n)%st_name)//'_import.nc', &
                  regridmethod=cap_fld_list(n)%regrid_method, &
                  srcFileType=ESMF_FILEFORMAT_SCRIP, &
                  dstFileType=ESMF_FILEFORMAT_SCRIP, &
                  srcRegionalFlag=.true., &
                  dstRegionalFlag=.true., &
                  normType=ESMF_NORMTYPE_FRACAREA, &
                  rc=rc)
             call check(rc, __LINE__, file)
             ! Precompute Field sparse matrix multiplication with local factors
             call ESMF_FieldSMMStore(srcField=meshfield, dstField=cap_fld_list(n)%import_field, &
                  filename=weights_dir//trim(cap_fld_list(n)%st_name)//'_import.nc', &
                  routehandle=cap_fld_list(n)%import_handle, &
                  ! transposeRoutehandle=cap_fld_list(n)%export_handle, &
                  rc=rc)
             call check(rc, __LINE__, file)

             cap_fld_list(n)%import_handle_init = .true.
          end if


          if (debug) print*, "DEBUGGING: using var's handle to regrid field ", &
               trim(cap_fld_list(n)%st_name)
          call ESMF_FieldRegrid(meshField, cap_fld_list(n)%import_field, &
               cap_fld_list(n)%import_handle, rc=rc)
          call check(rc, __LINE__, file)

          ! --- Debugging: write import mesh after regrid
          if (show_import_once) then
             call create_dir_if_needed(vars_in_dir)
             call ESMF_FieldWrite(meshField, &
                  vars_in_dir//trim(cap_fld_list(n)%st_name)// &
                  "_mesh_post.nc", &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call check(rc, __LINE__, file)
             call ESMF_FieldWrite(cap_fld_list(n)%import_field, &
                  vars_in_dir//trim(cap_fld_list(n)%st_name)// &
                  "_grid_post.nc", &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call check(rc, __LINE__, file)
          end if
          ! --- end debug ---
          ! call ESMF_FieldDestroy(gridField, rc=rc)
          call check(rc, __LINE__, file)
       end if
    end do
    show_import_once = .false.
    ! stop "look at import meshes"
  end subroutine regrid_import_mesh_to_grid

  subroutine regrid_export_grid_to_mesh(grid, mesh, state, &
       saved_import_state, did, memflg)
    use MPI
    type(ESMF_Grid), intent(in) :: grid
    type(ESMF_Mesh), intent(in) :: mesh
    type(ESMF_State),intent(in) :: state
    type(ESMF_State),intent(in) :: saved_import_state
    integer, intent(in) :: did
    type(memory_flag), intent(in) :: memflg
    type(ESMF_Field) :: meshField!, gridField
    type(ESMF_Field) :: importMeshField
    integer :: n, rc, itemCount
    character(len=64), allocatable :: itemNameList(:)
    logical :: exported

    real(ESMF_KIND_R8), pointer :: gridPtr(:,:) => null()
    type(ESMF_Array) :: gridarray
    real(ESMF_KIND_R8), pointer :: r82d(:,:) => null()
    real(ESMF_KIND_R4), pointer :: r42d(:,:) => null()
    real(ESMF_KIND_FIELD), pointer :: farrayPtr2d(:,:)


    character(:), allocatable :: mpas_grid_file, scrip_mesh_file, hires_scrip_file

    ! debugging
    integer :: unit
    logical :: exists

    integer :: rank, localDECount, localDE
    type(ESMF_TypeKind_Flag) :: tk

    real(ESMF_KIND_R8), pointer :: import_mesh(:)
    real(ESMF_KIND_R8), pointer :: export_mesh(:)


    integer :: ierr

    if (print_export_statelog_once) then
       if (io_rank) print *, 'WRFH: enter regrid_export_grid_to_mesh'
       if (io_rank) call ESMF_LogWrite('--- WRFH: EXPORT STATE DEBUG: ', ESMF_LOGMSG_INFO)
       if (debug) call ESMF_LogWrite('___ EXPORT STATE DEBUG: ', ESMF_LOGMSG_INFO)
       if (debug) call ESMF_StatePrint(state, rc=rc)
       print_export_statelog_once = .false.
    end if

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)
    allocate(itemNameList(itemCount), stat=rc)

    call ESMF_StateGet(state,itemNameList=itemNameList, rc=rc)
    call check(rc, __LINE__, file)
    if (debug) print *, "============ Export itemnamelist: ", itemNameList

    do n=lbound(cap_fld_list,1),ubound(cap_fld_list,1)
       if (cap_fld_list(n)%ad_export) then
          exported = NUOPC_IsConnected(state, &
               fieldName=trim(cap_fld_list(n)%st_name), rc=rc)
          if (.not. exported) then
             error stop "Variable " // trim(cap_fld_list(n)%st_name) // " should be exported"
          end if
          call check(rc, __LINE__, file)

          if (debug) print *, "exported = ", exported,", name: ", cap_fld_list(n)%st_name
          if (debug) call ESMF_LogWrite("WRFH: exported "//cap_fld_list(n)%st_name, &
               ESMF_LOGMSG_INFO, rc=rc)

          ! FOOBAR LOOK HERE
          ! meshField = field_create(cap_fld_list(n)%st_name, mesh, did, &
          !      memflg, rc)
          ! call check(rc, __LINE__, file)
          ! will this work better??
          ! THIS LOOKS BAD FOR NP=2

          ! sfhead is export only, so doen't have an import mesh to copy
          if (cap_fld_list(n)%ad_import) then
             call ESMF_StateGet(saved_import_state, &
                  itemName=trim(cap_fld_list(n)%st_name), &
                  field=importMeshField, rc=rc)
             call check(rc, __LINE__, file)
             call ESMF_FieldGet(importMeshField,   farrayPtr=import_mesh, rc=rc)
             call check(rc, __LINE__, file)

             call ESMF_StateGet(state, &
                  itemName=trim(cap_fld_list(n)%st_name), &
                  field=meshField, rc=rc)
             call check(rc, __LINE__, file)
             call ESMF_FieldGet(meshField,   farrayPtr=export_mesh, rc=rc)
             call check(rc, __LINE__, file)

             export_mesh(:) = import_mesh(:)
          else
             call ESMF_StateGet(state, &
                  itemName=trim(cap_fld_list(n)%st_name), &
                  field=meshField, rc=rc)
             call check(rc, __LINE__, file)
          end if

          ! new way to handle gridField
          if (cap_fld_list(n)%export_field_init .eqv. .false.) then
             cap_fld_list(n)%export_field = field_create(cap_fld_list(n)%st_name, grid, did, &
                  memflg, rc)
             call check(rc, __LINE__, file)
             cap_fld_list(n)%export_field_init = .true.
          end if

          ! Debugging: write export mesh before regrid
          ! if (debug_importexport_vars) then
          if (show_export_once) then
             call create_dir_if_needed(vars_out_dir)

             fname = vars_out_dir &
                  //trim(cap_fld_list(n)%st_name) &
                  // "_mesh_pre" &
                  ! //trim(scount) &
                  //".nc"

             call ESMF_FieldWrite(meshField, &
                  fname, &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call check(rc, __LINE__, file)
             fname = vars_out_dir &
                  //trim(cap_fld_list(n)%st_name) &
                  // "_grid_pre" &
                  ! //trim(scount) &
                  //".nc"
             call ESMF_FieldWrite(cap_fld_list(n)%export_field, &
                  fname, &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call check(rc, __LINE__, file)
          end if


          if (cap_fld_list(n)%export_handle_init .eqv. .false.) then
             call printa("Regridding export variables")
             ! call ESMF_FieldRegridStore(srcField=gridField, &
             ! call ESMF_RegridWeightGen(&
             ! ! Precompute Field sparse matrix multiplication with local factors
             ! call ESMF_FieldSMMStore(srcField=gridField, dstField=meshField, &


             ! ESMF_UNMAPPEDACTION_IGNORE:
             ! Destination points which do not have corresponding source
             ! points are ignored and zeros are used for the entries of the
             ! sparse matrix that is generated.
             ! ESMF_REGION_EMPTY:
             ! the elements in dstField will not be modified prior to the
             ! sparse matrix multiplication and results will be added to the
             ! incoming element values
             ! ESMF_REGION_SELECT:
             ! only zero out those elements in the destination Field that
             ! will be updated by the sparse matrix multiplication
             if (.not. allocated(hires_file)) &
                  error stop "hires_file variable not defined"
             hires_scrip_file = fulldom_to_scrip_filename(hires_file)

             ! hires_scrip_file needs grid corner, grid
             ! center, dims and imask

             mpas_grid_file = get_mpas_grid_filename()
             scrip_mesh_file = mpas_to_scrip_filename(mpas_grid_file)

             call ESMF_RegridWeightGen(&
                  srcFile=hires_scrip_file, & ! grid
                  dstFile=scrip_mesh_file, &            ! to mesh
                  weightFile=weights_dir//trim(cap_fld_list(n)%st_name)//'_export.nc', &
                  regridmethod=cap_fld_list(n)%regrid_method, &
                  unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                  srcFileType=ESMF_FILEFORMAT_SCRIP, &
                  dstFileType=ESMF_FILEFORMAT_SCRIP, &
                  srcRegionalFlag=.true., &
                  dstRegionalFlag=.true., &
                  normType=ESMF_NORMTYPE_FRACAREA, &
                  rc=rc)
             call check(rc, __LINE__, file)
             ! Precompute Field sparse matrix multiplication with local factors
             call ESMF_FieldSMMStore(srcField=cap_fld_list(n)%export_field, dstField=meshfield, &
                  filename=weights_dir//trim(cap_fld_list(n)%st_name)//'_export.nc', &
                  routehandle=cap_fld_list(n)%export_handle, &
                  rc=rc)
             call check(rc, __LINE__, file)

             cap_fld_list(n)%export_handle_init = .true.
          end if

          if (debug) print*, "DEBUGGING: using var's handle to regrid field ", &
               trim(cap_fld_list(n)%st_name)
          call ESMF_FieldRegrid(cap_fld_list(n)%export_field, meshField, &
               cap_fld_list(n)%export_handle, &
               zeroregion=ESMF_REGION_SELECT, rc=rc)
          call check(rc, __LINE__, file)


          if (show_export_once) then
             fname = vars_out_dir &
                  //trim(cap_fld_list(n)%st_name) &
                  // "_mesh_post" &
                  ! //trim(scount) &
                  //".nc"

             call ESMF_FieldWrite(meshField, &
                  fname, &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call check(rc, __LINE__, file)

             fname = vars_out_dir &
                  //trim(cap_fld_list(n)%st_name) &
                  // "_grid_post" &
                  ! //trim(scount) &
                  //".nc"
             call ESMF_FieldWrite(cap_fld_list(n)%export_field, &
                  fname, &
                  overwrite=.true., &
                  variableName=trim(cap_fld_list(n)%st_name), rc=rc)
             call check(rc, __LINE__, file)
          end if

          ! call ESMF_FieldDestroy(gridField, rc=rc)
          call check(rc, __LINE__, file)
       end if
    end do

    ! stop "test"
    show_export_once = .false.
  end subroutine regrid_export_grid_to_mesh


  subroutine dump_state(state, label)
    use ESMF
    implicit none
    type(ESMF_State), intent(in) :: state
    character(len=*), intent(in) :: label

    integer :: rc, n, i, pet
    type(ESMF_VM) :: vm
    type(ESMF_StateItem_Flag), allocatable :: itemType(:)
    character(len=ESMF_MAXSTR), allocatable :: itemName(:)
    type(ESMF_StateItem_Flag) :: case
    type(ESMF_Field) :: f
    character(len=ESMF_MAXSTR) :: nm, stdn
    type(ESMF_State) :: sub

    call ESMF_VMGetCurrent(vm, rc=rc)
    call ESMF_VMGet(vm, localPet=pet, rc=rc)
    if (pet /= 0) return   ! print once

    call ESMF_StateGet(state, itemCount=n, rc=rc)
    if (n <= 0) then
       call ESMF_LogWrite(trim(label)//" state: (empty)", ESMF_LOGMSG_INFO)
       return
    endif

    allocate(itemType(n), itemName(n))
    call ESMF_StateGet(state, itemTypeList=itemType, itemNameList=itemName, rc=rc)

    do i = 1, n
       print*, "itemType(i) =", itemType(i)
       case = itemType(i)
       if (case == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(state, itemName=itemName(i), field=f, rc=rc)
          call ESMF_FieldGet(f, name=nm, rc=rc)
          stdn = ""
          call ESMF_AttributeGet(f, name="standard_name", value=stdn, defaultValue="", rc=rc)
          call ESMF_LogWrite("Import Field: "//trim(nm)// &
               merge(" ["//trim(stdn)//"]","",len_trim(stdn)>0), ESMF_LOGMSG_INFO)

          ! case (ESMF_STATEITEM_FIELDBUNDLE)
          !   type(ESMF_FieldBundle) :: fb
          !   integer :: k, nf
          !   type(ESMF_Field), allocatable :: flist(:)
          !   character(len=ESMF_MAXSTR) :: bname, fname
          !   call ESMF_StateGet(state, itemName=itemName(i), fieldbundle=fb, rc=rc)
          !   call ESMF_FieldBundleGet(fb, name=bname, fieldCount=nf, rc=rc)
          !   allocate(flist(nf))
          !   call ESMF_FieldBundleGet(fb, fieldList=flist, rc=rc)
          !   call ESMF_LogWrite("Import FieldBundle: "//trim(bname), ESMF_LOGMSG_INFO)
          !   do k=1,nf
          !     call ESMF_FieldGet(flist(k), name=fname, rc=rc)
          !     call ESMF_LogWrite("  - "//trim(fname), ESMF_LOGMSG_INFO)
          !   enddo
          !   deallocate(flist)

       ! else if (case == ESMF_STATEITEM_STATE) then
          ! call ESMF_StateGet(state, itemName=itemName(i), state=sub, rc=rc)
          ! call ESMF_LogWrite("Nested State: "//trim(itemName(i)), ESMF_LOGMSG_INFO)
          ! call dump_state(sub, "  "//trim(label))  ! recurse
       else
          call ESMF_LogWrite("Unknown item: "//trim(itemName(i)), ESMF_LOGMSG_INFO)
       end if
    enddo
  end subroutine dump_state


  function read_hires_filename_from_namelist() result(geo_finegrid_flnm)
    implicit none
    ! integer :: sys_cpl
    ! character(len=256) :: geo_finegrid_flnm, geo_static_flnm
    ! namelist /hydro_nlist/ sys_cpl, geo_static_flnm, geo_finegrid_flnm
    character(len=256) :: geo_finegrid_flnm

    integer  :: NSOIL
    real :: ZSOIL8(8)
    type(crocus_options) :: crocus_opts

    logical            :: dir_e
    character(len=1024) :: reservoir_obs_dir
    integer:: RT_OPTION, CHANRTSWCRT, channel_option, &
         SUBRTSWCRT,OVRTSWCRT,AGGFACTRT, &
         GWBASESWCRT,  GW_RESTART,RSTRT_SWC,TERADJ_SOLAR, &
         sys_cpl, rst_typ, rst_bi_in, rst_bi_out, &
         gwChanCondSw, GwPreCycles, GwSpinCycles, GwPreDiagInterval, &
         gwsoilcpl, UDMP_OPT, io_form_outputs, bucket_loss, imperv_adj
    real:: DTRT_TER,DTRT_CH,dxrt, gwChanCondConstIn, gwChanCondConstOut, &
         gwIhShift
    character(len=256) :: route_topo_f="", route_chan_f="", route_link_f=""
    integer            :: channel_loss_option = 0, lake_option
    character(len=256) :: route_lake_f=""
    character(len=256) :: diversions_file=""
    character(len=256) :: reservoir_parameter_file=""
    character(len=256) :: reservoir_usgs_timeslice_path=""
    character(len=256) :: reservoir_usace_timeslice_path=""
    integer :: reservoir_observation_lookback_hours = 24, &
         reservoir_observation_update_time_interval_seconds = 86400, &
         reservoir_rfc_forecasts_lookback_hours = 28
    character(len=256) :: reservoir_rfc_forecasts_time_series_path=""
    character(len=256) :: route_direction_f="", route_order_f="", &
         gwbasmskfil ="", gwstrmfil ="", udmap_file ="", GWBUCKPARM_file = ""
    integer :: SOLVEG_INITSWC, reservoir_data_ingest
    real*8 :: out_dt, rst_dt
    character(len=256)  :: RESTART_FILE = "", hydrotbl_f   = ""
    logical            :: GwPreDiag, GwSpinUp, compound_channel, &
         reservoir_persistence_usgs, reservoir_persistence_usace, &
         reservoir_rfc_forecasts, reservoir_type_specified
    integer :: split_output_count, order_to_write, igrid, io_config_outputs, &
         t0OutputFlag, output_channelBucket_influx
    character(len=256) :: geo_static_flnm = "", land_spatial_meta_flnm = ""
    integer  :: DEEPGWSPIN, CHRTOUT_DOMAIN, CHRTOUT_GRID, LSMOUT_DOMAIN, &
         RTOUT_DOMAIN, output_gw, outlake, frxst_pts_out, CHANOBS_DOMAIN

    namelist /hydro_nlist/ nsoil, zsoil8,&
         restart_file,split_output_count,igrid,&
         geo_static_flnm, &
         land_spatial_meta_flnm, &
         out_dt, rst_dt, &
         deepgwspin, solveg_initswc, &
         rt_option, chanrtswcrt, channel_option, &
         subrtswcrt,ovrtswcrt,aggfactrt, dtrt_ter,dtrt_ch,dxrt,&
         gwspincycles, gwprecycles, gwspinup, gwprediag, gwprediaginterval, &
         gwihshift, &
         gwbaseswcrt, gwchancondsw, gwchancondconstin, gwchancondconstout , &
         route_topo_f,route_chan_f,route_link_f, compound_channel, &
         channel_loss_option, &
         lake_option, route_lake_f,  diversions_file, &
         route_direction_f,route_order_f,gwbasmskfil, &
         geo_finegrid_flnm, gwstrmfil,gw_restart,rstrt_swc,teradj_solar, sys_cpl, &
         order_to_write , rst_typ, rst_bi_in, rst_bi_out, gwsoilcpl, &
         chrtout_domain,chanobs_domain,chrtout_grid,lsmout_domain,&
         rtout_domain, output_gw, outlake, &
         frxst_pts_out, udmap_file, udmp_opt, gwbuckparm_file, bucket_loss, &
         io_config_outputs, io_form_outputs, hydrotbl_f, t0outputflag, &
         output_channelbucket_influx, imperv_adj


    integer :: unit, ios
    character(len=512) :: iomsg

    geo_finegrid_flnm = ""

    open(newunit=unit, file="hydro.namelist", status="old", action="read", &
         iostat=ios)
    if (ios /= 0) error stop "Error opening hydro.namelist"

    read(unit, nml=hydro_nlist, iostat=ios, iomsg=iomsg)
    close(unit)

    if (ios /= 0) then
       ! print *, "iomsg = ",
       error stop trim(iomsg)
    end if
  end function read_hires_filename_from_namelist

  function wrfhydro_grid_create_from_fulldom(rc) &
       result(grid)
    implicit none
    type(ESMF_Grid) :: grid
    integer, intent(out) :: rc
    integer :: nx, ny
    real :: lon0, lat0, dlon, dlat
    real, allocatable :: lats(:), lons(:)
    real :: lat_start, lat_end, lon_start, lon_end
    real(ESMF_KIND_R8), pointer :: lon(:,:), lat(:,:)
    integer :: i,j

    integer :: ncid, dim_xid, dim_yid
    integer :: var_lat_id, var_lon_id
    integer :: stat
    real, allocatable :: latitude(:,:), longitude(:,:)

    ! read from hydro.namelist
    hires_file = read_hires_filename_from_namelist()

    ! get lat/lon and nx/ny from hires_file
    print *, "WRFH: Opening high-resolution file ", hires_file
    stat = nf90_open(trim(hires_file), NF90_NOWRITE, ncid)
    call check_nf(stat)


    ! Dimension IDs and sizes
    stat = nf90_inq_dimid(ncid, 'x', dim_xid)
    call check_nf(stat)
    stat = nf90_inq_dimid(ncid, 'y', dim_yid)
    call check_nf(stat)
    stat = nf90_inquire_dimension(ncid, dim_xid, len=nx)
    call check_nf(stat)
    stat = nf90_inquire_dimension(ncid, dim_yid, len=ny)
    call check_nf(stat)

    if (nx <= 0 .or. ny <= 0) then
       stat = NF90_EBADID ! why: surface a clear failure if dims are invalid
       print *, "Error: ", trim(hires_file), " nx <=0 or ny <= 0"
       error stop "High resolution file's nx <=0 or ny <= 0"
    end if

    print *, trim(hires_file), " nx, ny:", nx, ny

    ! Variable IDs
    stat = nf90_inq_varid(ncid, 'LATITUDE', var_lat_id)
    call check_nf(stat)
    stat = nf90_inq_varid(ncid, 'LONGITUDE', var_lon_id)
    call check_nf(stat)

    ! Read variables (library converts to real64 as needed)
    allocate(latitude(nx, ny), longitude(nx, ny))
    stat = nf90_get_var(ncid, var_lat_id, latitude)
    call check_nf(stat)
    stat = nf90_get_var(ncid, var_lon_id, longitude)
    call check_nf(stat)

    ! print *, "lat shape =", shape(latitude)
    ! print *, "lat(1:2,1:2) =", latitude(1:2,1:2)

    grid = ESMF_GridCreate(&
         regDecomp=[1, 1], &
         decompflag=[ESMF_DECOMP_BALANCED, ESMF_DECOMP_BALANCED], &
         maxIndex=[nx,ny], &
         rc=rc)
    call check(rc, __LINE__, file)

    ! add center coordinates to grid
    call ESMF_GridAddCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check(rc, __LINE__, file)

    ! Get writable pointers to the coordinate arrays and fill a regular lat/lon
    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
         coordDim=1, farrayPtr=lon, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
         coordDim=2, farrayPtr=lat, rc=rc)
    call check(rc, __LINE__, file)

    lat(:,:) = latitude(:,:)
    lon(:,:) = longitude(:,:)
  end function wrfhydro_grid_create_from_fulldom

#undef METHOD
#define METHOD "wrfhydro_GridCreate"
  function wrfhydro_GridCreate(did,rc) result(wrfhydro_grid)
    type(ESMF_Grid) :: wrfhydro_grid
    ! ARGUMENTS
    integer, intent(in)                     :: did
    integer, intent(out)                    :: rc
    ! LOCAL VARIABLES
    integer                     :: stat
    real                        :: min_lat, max_lat, min_lon, max_lon
    real, allocatable           :: latitude(:,:), longitude(:,:)
    integer, allocatable        :: mask(:,:)
    integer                     :: lbnd(2),ubnd(2)
    real(ESMF_KIND_COORD), pointer :: coordXcenter(:,:)
    real(ESMF_KIND_COORD), pointer :: coordYcenter(:,:)
    real(ESMF_KIND_COORD), pointer :: coordXcorner(:,:)
    real(ESMF_KIND_COORD), pointer :: coordYcorner(:,:)
    integer(ESMF_KIND_I4), pointer :: gridmask(:,:)
    integer                     :: i,j, i1,j1
    character(len=16)           :: xlat_corner_name, xlon_corner_name
#ifdef DEBUG
    character(ESMF_MAXSTR)      :: logMsg
#endif
    character(:), allocatable :: grid_file

    integer :: ierr, rank


#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    grid_file = 'WRFHYDRO_Grid_'//trim(nlst(did)%hgrid)
    ! grid_file = 'Fulldom_hires.nc'
    print *, "move wrfhydro grid file to namelist: ", grid_file

    ! print *, "wrfhydro_distGrid =", wrfhydro_distGrid
    ! stop "pink panther"
    if (wrfhydro_distgrid_initialized .eqv. .false.) then
       error stop "Make sure wrfhydro_distgrid_initalized is initialized"
    end if
    wrfhydro_grid = ESMF_GridCreate(name=grid_file, &
         distgrid=wrfhydro_distGrid, &
         coordSys=ESMF_COORDSYS_SPH_DEG, &
         coordTypeKind=ESMF_TYPEKIND_COORD, &
         ! gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
         rc = rc)
    call check(rc, __LINE__, file)

    ! CENTERS

    ! Get Local Latitude (lat)
    allocate(latitude(nx_local,ny_local),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=METHOD//': Allocation of latitude memory failed.', &
      file=FILENAME, rcToReturn=rc)) return ! bail out
    call WRFHYDRO_ESMF_NetcdfReadIXJX("XLAT_M",nlst(did)%geo_static_flnm, &
      (/x_start,y_start/),latitude,rc=rc)
    call check(rc, __LINE__, file)

    ! Get Local Longitude (lon)
    allocate(longitude(nx_local,ny_local),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=METHOD//': Allocation of longitude memory failed.', &
      file=FILENAME, rcToReturn=rc)) return ! bail out
    call WRFHYDRO_ESMF_NetcdfReadIXJX("XLONG_M",nlst(did)%geo_static_flnm, &
      (/x_start,y_start/),longitude,rc=rc)
    call check(rc, __LINE__, file)

! #ifdef DEBUG
    ! Print Local Lat Lon Lower Left / Upper Right Centers
    write(logMsg,"(A,4(F0.3,A))") MODNAME//": Center Coordinates = (", &
      longitude(1,1),":",longitude(nx_local,ny_local),",", &
      latitude(1,1),":",latitude(nx_local,ny_local),")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
! #endif

    ! Add Center Coordinates to Grid
    call ESMF_GridAddCoord(wrfhydro_grid, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check(rc, __LINE__, file)

    call ESMF_GridGetCoord(wrfhydro_grid, coordDim=1, localDE=0, &
      staggerloc=ESMF_STAGGERLOC_CENTER, &
      computationalLBound=lbnd, computationalUBound=ubnd, &
      farrayPtr=coordXcenter, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_GridGetCoord(wrfhydro_grid, coordDim=2, localDE=0, &
      staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordYcenter, rc=rc)
    call check(rc, __LINE__, file)

    do j = lbnd(2),ubnd(2)
    do i = lbnd(1),ubnd(1)
      coordXcenter(i,j) = longitude(i,j)
      coordYcenter(i,j) = latitude(i,j)
    enddo
    enddo

    deallocate(latitude,longitude,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=METHOD//': Deallocation of longitude and latitude memory failed.', &
      file=FILENAME,rcToReturn=rc)) return ! bail out

    ! Get Local Mask
    allocate(mask(nx_local,ny_local),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=METHOD//': Allocation of mask memory failed.', &
      file=FILENAME, rcToReturn=rc)) return ! bail out
    call WRFHYDRO_ESMF_NetcdfReadIXJX("LANDMASK",nlst(did)%geo_static_flnm, &
      (/x_start,y_start/),mask,rc=rc)
    call check(rc, __LINE__, file)

    ! Add Grid Mask
    call ESMF_GridAddItem(wrfhydro_grid, itemFlag=ESMF_GRIDITEM_MASK, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check(rc, __LINE__, file)
    ! Get pointer to Grid Mask array
    call ESMF_GridGetItem(wrfhydro_grid, itemflag=ESMF_GRIDITEM_MASK, &
      localDE=0, &
      staggerloc=ESMF_STAGGERLOC_CENTER, &
      farrayPtr=gridmask, rc=rc)
    call check(rc, __LINE__, file)

    do j = lbnd(2),ubnd(2)
    do i = lbnd(1),ubnd(1)
      gridmask(i,j) = mask(i,j)
      gridmask(i,j) = mask(i,j)
    enddo
    enddo

    deallocate(mask,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=METHOD//': Deallocation of mask memory failed.', &
      file=FILENAME,rcToReturn=rc)) return ! bail out

    ! CORNERS
    ! The original WPS implementation used the _CORNER names
    ! but it was then changes to the _C names.  Support both
    ! options.
    if (WRFHYDRO_ESMF_NetcdfIsPresent("XLAT_CORNER",nlst(did)%geo_static_flnm) .AND. &
         WRFHYDRO_ESMF_NetcdfIsPresent("XLONG_CORNER",nlst(did)%geo_static_flnm)) then
       xlat_corner_name = "XLAT_CORNER"
       xlon_corner_name = "XLONG_CORNER"
    else if (WRFHYDRO_ESMF_NetcdfIsPresent("XLAT_C",nlst(did)%geo_static_flnm) .AND. &
         WRFHYDRO_ESMF_NetcdfIsPresent("XLONG_C",nlst(did)%geo_static_flnm)) then
       xlat_corner_name = "XLAT_C"
       xlon_corner_name = "XLONG_C"
    else
       xlat_corner_name = ""
       xlon_corner_name = ""
    endif

    if (trim(xlat_corner_name) /= "") then
      ! Get Local Latitude (lat)
      allocate(latitude(nx_local+1,ny_local+1),stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg=METHOD//': Allocation of corner latitude memory failed.', &
        file=FILENAME, rcToReturn=rc)) return ! bail out
      call WRFHYDRO_ESMF_NetcdfReadIXJX(trim(xlat_corner_name),nlst(did)%geo_static_flnm, &
        (/x_start,y_start/),latitude,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out

      ! Get Local Longitude (lon)
      allocate(longitude(nx_local+1,ny_local+1),stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
       msg=METHOD//': Allocation of corner longitude memory failed.', &
       file=FILENAME, rcToReturn=rc)) return ! bail out
      call WRFHYDRO_ESMF_NetcdfReadIXJX(trim(xlon_corner_name),nlst(did)%geo_static_flnm, &
        (/x_start,y_start/),longitude,rc=rc)
      if(ESMF_STDERRORCHECK(rc)) return ! bail out

#ifdef DEBUG
      ! Print Local Lat Lon Lower Left / Upper Right Corners
      write(logMsg,"(A,4(F0.3,A))") MODNAME//": Corner Coordinates = (", &
        longitude(1,1),":",longitude(nx_local+1,ny_local+1),",", &
        latitude(1,1),":",latitude(nx_local+1,ny_local+1),")"
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
#endif

      ! Add Corner Coordinates to Grid
      call ESMF_GridAddCoord(wrfhydro_grid, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc)
      call check(rc, __LINE__, file)

      call ESMF_GridGetCoord(wrfhydro_grid, coordDim=1, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CORNER, &
        computationalLBound=lbnd, computationalUBound=ubnd, &
        farrayPtr=coordXcorner, rc=rc)
      call check(rc, __LINE__, file)
      call ESMF_GridGetCoord(wrfhydro_grid, coordDim=2, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=coordYcorner, rc=rc)
      call check(rc, __LINE__, file)

      do j = lbnd(2),ubnd(2)
      do i = lbnd(1),ubnd(1)
        coordXcorner(i,j) = longitude(i,j)
        coordYcorner(i,j) = latitude(i,j)
      enddo
      enddo

      deallocate(latitude,longitude,stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg=METHOD//': Deallocation of corner longitude and latitude memory failed.', &
        file=FILENAME,rcToReturn=rc)) return ! bail out

      call add_area(wrfhydro_grid, rc=rc)
      call check(rc, __LINE__, file)

    else
#ifdef DEBUG
      ! Warning no corners in domain file
      call ESMF_LogWrite(MODNAME//": No Corner Coordinates.", ESMF_LOGMSG_WARNING)
#endif
    endif

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif



  end function wrfhydro_GridCreate


  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "add_area"

  subroutine add_area(grid,rc)
    type(ESMF_Grid), intent(inout)          :: grid
    integer, intent(out)                    :: rc

    ! Local Variables
    integer(ESMF_KIND_I4), PARAMETER :: R = 6376000 ! metres
    type(ESMF_Field)                 :: fieldArea
    type(ESMF_Array)                 :: areaArray
    integer                          :: i,j
    integer                          :: lbnd(2),ubnd(2)
    real(ESMF_KIND_R8), pointer      :: radianarea(:,:)
    real(ESMF_KIND_R8), pointer      :: gridarea(:,:)

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    fieldArea = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldRegridGetArea(fieldArea, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_FieldGet(fieldArea, localDE=0, &
      farrayPtr=radianarea, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_GridAddItem(grid, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
      staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

    call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, localDE=0, &
      staggerloc=ESMF_STAGGERLOC_CENTER, &
      computationalLBound=lbnd, computationalUBound=ubnd, &
      farrayPtr=gridarea, rc=rc)
    if (ESMF_STDERRORCHECK(rc)) return

     do j = lbnd(2),ubnd(2)
     do i = lbnd(1),ubnd(1)
       gridarea(i,j) = radianarea(i,j) * R * R
     enddo
     enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "set_local_indices"

  subroutine set_local_indices(rc)
    ! ARGUMENTS
    integer, intent(out)                    :: rc

    ! LOCAL VARIABLES
    integer                     :: stat
    type(ESMF_VM)               :: currentVM
    integer                     :: localPet
    integer                     :: petCount
    type(ESMF_DELayout)         :: delayout
    integer, allocatable        :: dimExtent(:,:)
    integer, allocatable        :: iIndexList(:), jIndexList(:)
#ifdef DEBUG
    character(ESMF_MAXSTR)      :: logMsg
#endif

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    !! Get VM Info to see if this will give me the PET info I need
    call ESMF_VMGetCurrent(currentVM, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    call ESMF_VMGet(currentVM, localPet=localPet, petCount=petCount, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    !! Get the grid distribution for this pet
    allocate(dimExtent(2, 0:(petCount - 1)),stat=stat) ! (dimCount, deCount)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=METHOD//': Allocation of indexCountPDe memory failed.', &
      file=FILENAME, rcToReturn=rc)) return ! bail out
    call ESMF_DistGridGet(wrfhydro_distGrid, delayout=delayout, &
      indexCountPDe=dimExtent, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    allocate(iIndexList(dimExtent(1, localPet)),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=METHOD//': Allocation of iIndexList memory failed.', &
      file=FILENAME, rcToReturn=rc)) return ! bail out
    call ESMF_DistGridGet(wrfhydro_distGrid, localDe=0, dim=1, &
      indexList=iIndexList, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    allocate(jIndexList(dimExtent(2, localPet)),stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg=METHOD//': Allocation of jIndexList memory failed.', &
      file=FILENAME, rcToReturn=rc)) return ! bail out
    call ESMF_DistGridGet(wrfhydro_distGrid, localDe=0, dim=2, &
      indexList=jIndexList, rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    x_start = minVal(iIndexList)
    x_end   = maxVal(iIndexList)
    y_start = minVal(jIndexList)
    y_end   = maxVal(jIndexList)

    nx_local = x_end - x_start + 1
    ny_local = y_end - y_start + 1

#ifdef DEBUG
    write (logMsg,"(A,6(I0,A))") MODNAME//": Local Indices = (", &
      x_start,":",x_end,",",y_start,":",y_end,") Local Size = (", &
      nx_local,"x",ny_local,")"
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
#endif

    deallocate(iIndexList,jIndexList,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=METHOD//': Deallocation of IndexList memory failed.', &
      file=FILENAME,rcToReturn=rc)) return ! bail out

    deallocate(dimExtent,stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg=METHOD//': Deallocation of indexCountPDeo memory failed.', &
      file=FILENAME,rcToReturn=rc)) return ! bail out

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_get_timestep"

  function WRFHYDRO_get_timestep(did,rc)
    ! RETURN VALUE
    real :: WRFHYDRO_get_timestep
    ! ARGUMENTS
    integer, intent(in)         :: did
    integer, intent(out)        :: rc

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    WRFHYDRO_get_timestep = nlst(did)%dt

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end function

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_set_timestep"

  subroutine WRFHYDRO_set_timestep(did,dt,rc)
    ! ARGUMENTS
    integer, intent(in)           :: did
    real                          :: dt
    integer, intent(out)          :: rc

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    nlst(did)%dt = dt

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_get_hgrid"

  subroutine WRFHYDRO_get_hgrid(did,hgrid,rc)
    ! ARGUMENTS
    integer, intent(in)         :: did
    character, intent(out)      :: hgrid
    integer, intent(out)        :: rc

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    hgrid = nlst(did)%hgrid

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_get_restart"

  subroutine WRFHYDRO_get_restart(did,restart,rc)
    ! ARGUMENTS
    integer, intent(in)         :: did
    logical, intent(out)        :: restart
    integer, intent(out)        :: rc

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    if (nlst(did)%rst_typ .eq. 0) then
      restart = .FALSE.
    else
      restart = .TRUE.
    endif

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------
  ! Conversion Utilities
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_ClockToString"

  subroutine WRFHYDRO_ClockToString(clock, timestr, rc)
    ! ARGUMENTS
    type(ESMF_Clock)                :: clock
    integer, intent(out),optional   :: rc
    character (len=*), intent(out)  :: timestr

    ! LOCAL VARIABLES
    type(ESMF_Time)            :: currTime

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(rc)) rc = ESMF_SUCCESS  ! Initialize

    ! Get the current time from the clock
    call ESMF_ClockGet(clock=clock,currTime=currTime,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    call WRFHYDRO_TimeToString(currTime,timestr,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

!-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_TimeToString"

  subroutine WRFHYDRO_TimeToString(time, timestr, rc)
    ! ARGUMENTS
    type(ESMF_Time)                 :: time
    integer, intent(out),optional   :: rc
    character (len=*), intent(out)  :: timestr

    ! LOCAL VARIABLES
    character (len=256)        :: tmpstr = ''
    integer                    :: strlen

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(rc)) rc = ESMF_SUCCESS  ! Initialize

    timestr = '' ! clear string

    if (len(timestr) < 19) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg=METHOD//": Time string is too short!", &
        file=FILENAME,rcToReturn=rc)
      return  ! bail out
    endif

    CALL ESMF_TimeGet(time,timeString=tmpstr,rc=rc )
    if(ESMF_STDERRORCHECK(rc)) return ! bail out

    strlen = min(len(timestr),len_trim(tmpstr))
    timestr(1:strlen) = tmpstr(1:strlen)
    timestr(11:11) = '_'

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_TimeIntervalGetReal"

  function WRFHYDRO_TimeIntervalGetReal(timeInterval,rc)
    ! RETURN VALUE:
    real                                :: WRFHYDRO_TimeIntervalGetReal
    ! ARGUMENTS
    type(ESMF_TimeInterval),intent(in)  :: timeInterval
    integer, intent(out), optional      :: rc

    ! LOCAL VARIABLES
    real(ESMF_KIND_R8)                  :: s_r8

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    if(present(rc)) rc = ESMF_SUCCESS

    WRFHYDRO_TimeIntervalGetReal = -9999

    call ESMF_TimeIntervalGet(timeInterval,s_r8=s_r8,rc=rc)
    if(ESMF_STDERRORCHECK(rc)) return ! bail out
    WRFHYDRO_TimeIntervalGetReal = s_r8

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end function

  !-----------------------------------------------------------------------------
  ! Log Utilities
  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_nlstLog"

  subroutine WRFHYDRO_nlstLog(did,label,rc)
    ! ARGUMENTS
    integer,intent(in)                   :: did
    character(len=*),intent(in),optional :: label
    integer,intent(out)                  :: rc

    ! LOCAL VARIABLES
    integer                     :: layerIndex
    character(len=64)           :: l_label

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    if (present(label)) then
      l_label = label
    else
      l_label = METHOD
    endif

    write (logMsg,"(A,I0)") ": Domain ID      = ",did
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)

    write (logMsg,"(A,5(I0,A))") ": Start Date     = ", &
      nlst(did)%START_YEAR,"-",nlst(did)%START_MONTH,"-", &
      nlst(did)%START_DAY,"_",nlst(did)%START_HOUR,":", &
      nlst(did)%START_MIN
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,F0.3)") ": Timestep       = ",nlst(did)%dt
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,F0.3)") ": Output Step    = ",nlst(did)%out_dt
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,F0.3)") ": Restart Step   = ",nlst(did)%rst_dt
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,F0.3)") ": Ter Routing Step   = ",nlst(did)%dtrt_ter
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,F0.3)") ": Ch Routing Step   = ",nlst(did)%dtrt_ch
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)

    write (logMsg,"(A,I0)") ": Grid ID        = ",nlst(did)%igrid
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A)") ": Hydro Grid     = ",nlst(did)%hgrid
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A)") ": Geo Grid File  = ",nlst(did)%geo_static_flnm
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A)") ": Fine Grid File = ",nlst(did)%geo_finegrid_flnm
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A)") ": GW Basin File  = ",nlst(did)%gwbasmskfil
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)

    write (logMsg,"(A,I0)") ": Restart Type   = ",nlst(did)%rst_typ
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,A)") ": Restart file   = ",nlst(did)%restart_file
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Coupling       = ",nlst(did)%sys_cpl
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)

    write (logMsg,"(A,I0)") ": Channel RT     = ",nlst(did)%CHANRTSWCRT
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Subsurface RT  = ",nlst(did)%SUBRTSWCRT
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Overland RT    = ",nlst(did)%OVRTSWCRT
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": GW Baseflow RT = ",nlst(did)%GWBASESWCRT
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Routing Option = ",nlst(did)%RT_OPTION
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Channel Option = ",nlst(did)%channel_option
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Aggr Factor    = ",nlst(did)%AGGFACTRT
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": GW Restart     = ",nlst(did)%GW_RESTART
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": SWC Restart    = ",nlst(did)%RSTRT_SWC
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)

    write (logMsg,"(A,I0)") ": Soil Layers    = ",nlst(did)%nsoil
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    do layerIndex=1,nlst(did)%nsoil
      write (logMsg,"(A,I0,A,F0.3)") ": Soil layer depth (", &
        layerIndex,") = ",nlst(did)%ZSOIL8(layerIndex)
      call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    enddo

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine

  !-----------------------------------------------------------------------------

#undef METHOD
#define METHOD "WRFHYDRO_domainLog"

  subroutine WRFHYDRO_domainLog(did,label,rc)
    ! ARGUMENTS
    integer,intent(in)                   :: did
    character(len=*),intent(in),optional :: label
    integer,intent(out)                  :: rc

    ! LOCAL VARIABLES
    character(len=64)           :: l_label

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": entered "//METHOD, ESMF_LOGMSG_INFO)
#endif

    rc = ESMF_SUCCESS

    if (present(label)) then
      l_label = label
    else
      l_label = METHOD
    endif

    write (logMsg,"(A,I0)") ": Domain ID      = ",did
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)

    write (logMsg,"(A,L1)") ": Domain Init    = ",rt_domain(did)%initialized
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Domain IX      = ",rt_domain(did)%IX
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Domain JX      = ",rt_domain(did)%JX
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Domain IXRT    = ",rt_domain(did)%IXRT
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Domain JXRT    = ",rt_domain(did)%JXRT
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Domain Forc    = ",rt_domain(did)%FORC_TYP
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Max Links      = ",rt_domain(did)%NLINKS
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Num Lakes      = ",rt_domain(did)%NLAKES
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)
    write (logMsg,"(A,I0)") ": Num Basins     = ",rt_domain(did)%numbasns
    call ESMF_LogWrite(trim(l_label)//logMsg,ESMF_LOGMSG_INFO)

#ifdef DEBUG
    call ESMF_LogWrite(MODNAME//": leaving "//METHOD, ESMF_LOGMSG_INFO)
#endif

  end subroutine


  ! TODO: remove this function
  function init_distgrid(vm) result(distGrid)
    type(ESMF_VM), intent(in) :: vm
    type(ESMF_DistGrid) :: distGrid

    integer :: localPet, rc
    call ESMF_VMGet(vm, localPet=localPet, &
      mpiCommunicator=HYDRO_COMM_WORLD, rc=rc)
!     if(ESMF_STDERRORCHECK(rc)) return ! bail out

!     call ESMF_VMBroadcast(vm, startx, count=numprocs, rootPet=IO_id, rc=rc)
!     if(ESMF_STDERRORCHECK(rc)) return ! bail out
!     call ESMF_VMBroadcast(vm, starty, count=numprocs, rootPet=IO_id, rc=rc)
!     if(ESMF_STDERRORCHECK(rc)) return ! bail out
!     call ESMF_VMBroadcast(vm, local_nx_size, count=numprocs, rootPet=IO_id, rc=rc)
!     if(ESMF_STDERRORCHECK(rc)) return ! bail out
!     call ESMF_VMBroadcast(vm, local_ny_size, count=numprocs, rootPet=IO_id, rc=rc)
!     if(ESMF_STDERRORCHECK(rc)) return ! bail out

!     allocate(deBlockList(2,2,numprocs))
!     do i = 1, numprocs
!       deBlockList(:,1,i) = (/startx(i),starty(i)/)
!       deBlockList(:,2,i) = (/startx(i)+local_nx_size(i)-1, &
!                              starty(i)+local_ny_size(i)-1/)
! !      write (logMsg,"(A,I0,A,4(I0,A))") MODNAME//": deBlockList ", i, " = (", &
! !        deBlockList(1,1,i),":",deBlockList(1,2,i),",", &
! !        deBlockList(2,1,i),":",deBlockList(2,2,i),")"
! !      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
!     enddo

! !    allocate(connectionList(1),stat=stat)
! !    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
! !      msg=METHOD//': Allocation of connection list memory failed.', &
! !      file=FILENAME, rcToReturn=rc)) return ! bail out
! !    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
! !      tileIndexB=1, positionVector=(/nx_global(1), 0/), rc=rc)
! !    if (ESMF_STDERRORCHECK(rc)) return  ! bail out


!     ! Create DistGrid based on WRFHDYRO Config NX,NY
!     WRFHYDRO_distgrid = ESMF_DistGridCreate( &
!       minIndex=(/1,1/), maxIndex=(/nx_global(1),ny_global(1)/), &
! !     indexflag = ESMF_INDEX_DELOCAL, &
!      deBlockList=deBlockList, &
! !     deLabelList=deLabelList, &
! !     delayout=delayout, &
! !     connectionList=connectionList, &
!       rc=rc)
!     if(ESMF_STDERRORCHECK(rc)) return ! bail out

!     deallocate(deBlockList)

! !   deallocate(connectionList,stat=stat)
! !   if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
! !     msg=METHOD//': Deallocation of connection list memory failed.', &
! !     file=FILENAME,rcToReturn=rc)) return ! bail out

!     ! Get the Local Decomp Incides
!     call set_local_indices(rc)
!     if(ESMF_STDERRORCHECK(rc)) return ! bail out
  end function init_distgrid

  function get_config_block_decomp_file_prefix() &
       result(decomp_filename)
    character(len=:), allocatable :: decomp_filename
    ! copied from mpas/src/core_atmosphere/inc/namelist_defines.inc
    character (len=512) :: config_block_decomp_file_prefix = &
         'x1.40962.graph.info.part.'
    integer :: config_number_of_blocks = 0
    logical :: config_explicit_proc_decomp = .false.
    character (len=512) :: config_proc_decomp_file_prefix = &
         'graph.info.part.'
    integer :: unit, stat
    namelist /decomposition/ &
         config_block_decomp_file_prefix, &
         config_number_of_blocks, &
         config_explicit_proc_decomp, &
         config_proc_decomp_file_prefix

    open(newunit=unit, file='namelist.atmosphere', status='old', &
         action='read', iostat=stat)
    if (stat /= 0) stop 'Could not open namelist.atmosphere'
    read(unit, decomposition, iostat=stat)
    if (stat /= 0) stop 'Error reading namelist.atmosphere'
    decomp_filename = config_block_decomp_file_prefix
  end function get_config_block_decomp_file_prefix

  function get_mpas_grid_filename() &
       result(mpas_grid_file)
    character(len=:), allocatable :: mpas_grid_file, &
         config_block_decomp_file_prefix
    integer :: i

    config_block_decomp_file_prefix = get_config_block_decomp_file_prefix()

    i = index(trim(config_block_decomp_file_prefix), 'graph.info.part', &
         back=.true.)
    if (i == 0) then
       print *, "config_block_decomp_file_prefix =", &
            config_block_decomp_file_prefix
       error stop &
            "Error: graph.info.part not in config_block_decomp_file_prefix"
    end if
    mpas_grid_file = config_block_decomp_file_prefix(:i-1) // "grid.nc"
  end function get_mpas_grid_filename

  function get_mpas_dist_grid(vm, rc) &
       result(dist_grid)
    type(ESMF_VM), intent(inout) :: vm
    integer, intent(inout) :: rc
    type(ESMF_DistGrid) :: dist_grid
    character(len=:), allocatable :: mpas_grid_file, iomsg
    character(len=:), allocatable :: config_block_decomp_file_prefix
    character(len=512) :: mpas_graph_file
    integer :: unit, iostat, irank, localCount, ierr, tmp
    integer :: i, rank, np, idx, inode
    integer, allocatable :: gindex(:)
    rc = ESMF_SUCCESS

    call ESMF_VMGetGlobal(vm, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_VMGet(vm, localPet=rank, petCount=np, rc=rc)
    call check(rc, __LINE__, file)

    config_block_decomp_file_prefix = get_config_block_decomp_file_prefix()

    print *, "np =", np
    print *, "config_block_decomp_file_prefix=", config_block_decomp_file_prefix
    ! open .graph.info when np==1 or .graph.info.np when np>1
    if (np == 1) then
       i = index(trim(config_block_decomp_file_prefix), '.part', back=.true.)
       mpas_graph_file = config_block_decomp_file_prefix(:i-1)
       print *, "MPAS: opening mpas_graph_file =", mpas_graph_file
       open(newunit=unit, file=mpas_graph_file, status="old", &
            action="read", iostat=ierr)
       if (ierr /= 0) error stop "Failed to open frontrange.graph.info"

       read(unit, *, iostat=ierr) localCount, tmp   ! reads: first_int second_int
       if (ierr /= 0) error stop &
            "Failed to read first line of frontrange.graph.info"
       close(unit)
    else
       write(mpas_graph_file, '(A,I0)') &
            trim(config_block_decomp_file_prefix), np
       print *, "MPAS: opening mpas_graph_file =", mpas_graph_file
       open(newunit=unit, file=mpas_graph_file, &
            status='old', action='read', iostat=iostat, iomsg=iomsg)
       if (iostat /= 0) then
          print *, trim(iomsg)
          stop "Error opening [casename].graph.info.part.[np]"
       end if
       localCount = 0
       do
          read(unit, *, iostat=iostat) irank
          if (iostat /= 0) exit
          if (irank == rank) localCount = localCount + 1
       end do
    end if

    allocate(gindex(localCount))
    print *, rank, "/", np, ": with localCount =", localCount

    ! setup grid distribution
    if (np == 1) then
       inode = 1
       do inode = 1, localCount
          gindex(inode) = inode
       end do
    else
       rewind(unit)
       idx   = 0
       inode = 0

       do
          read(unit, *, iostat=iostat) irank
          if (iostat /= 0) exit

          inode = inode + 1 ! inode = line number
          if (irank == rank) then
             idx = idx + 1
             gindex(idx) = inode ! seqIndex = global node id
          end if
       end do
    end if
    close(unit)

    print *, rank,":gindex=", gindex(1:10)
    ! stop "here"
    dist_grid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
  end function get_mpas_dist_grid


  function mpas_to_scrip_filename(mpas_grid_file) result(scrip_mesh_file)
    character(len=*), intent(in) :: mpas_grid_file
    character(len=:), allocatable :: scrip_mesh_file
    integer :: i
    i = index(trim(mpas_grid_file), '.grid.nc', back=.true.)
    if (i == 0) then
       print *, "Error: mpas_grid_file ", trim(mpas_grid_file), &
            " does not have .grid.nc suffix"
       error stop "Error: mpas_grid_file in incorrect format"
    end if
    scrip_mesh_file = mpas_grid_file(:i-1) // ".tmp.scrip.nc"
  end function mpas_to_scrip_filename


  function fulldom_to_scrip_filename(fulldom_file) result(scrip_file)
    character(len=*), intent(in) :: fulldom_file
    character(len=:), allocatable :: scrip_file
    integer :: i
    i = index(trim(fulldom_file), '.nc', back=.true.)
    if (i == 0) then
       print *, "Error: fulldom_file ", trim(fulldom_file), &
            " does not have .nc suffix"
       error stop "Error: fulldom_file in incorrect format"
    end if
    scrip_file = fulldom_file(:i-1) // ".tmp.scrip.nc"
  end function fulldom_to_scrip_filename


  ! Convert MPAS mesh to scrip format. Based on
  ! mpas-dev.github.io/MPAS-Tools/0.24.0/_modules/mpas_tools/scrip/from_mpas.html
  subroutine mpas_to_scrip_mesh(mpasFile, scripFile, useLandIceMask)
    use netcdf
    use iso_fortran_env, only : real64, int32
    implicit none

    character(len=*), intent(in) :: mpasFile
    character(len=*), intent(in) :: scripFile
    logical, intent(in), optional :: useLandIceMask

    logical :: doLandIceMask
    integer :: stat
    integer :: fin, fout
    integer :: dimid, nCells, nVertices, maxVertices
    integer :: varid

    real(ESMF_KIND_R8), parameter :: SHR_CONST_REARTH = 6.37122e6

    ! input arrays
    real(ESMF_KIND_R8), allocatable :: latCell(:), lonCell(:)
    real(ESMF_KIND_R8), allocatable :: latVertex(:), lonVertex(:)
    real(ESMF_KIND_R8), allocatable :: areaCell(:)
    integer(ESMF_KIND_I4), allocatable :: verticesOnCell(:,:)
    integer(ESMF_KIND_I4), allocatable :: nEdgesOnCell(:)
    integer(ESMF_KIND_I4), allocatable :: landIceMask1d(:)
    integer(ESMF_KIND_I4), allocatable :: landIceMask2d(:,:)

    ! output arrays
    real(ESMF_KIND_R8), allocatable :: grid_corner_lat(:,:), grid_corner_lon(:,:)
    real(ESMF_KIND_R8), allocatable :: grid_area(:)
    integer(ESMF_KIND_I4), allocatable :: grid_imask(:)
    integer(ESMF_KIND_I4) :: grid_dims(1)

    ! output variable ids
    integer :: dim_grid_size, dim_grid_corners, dim_grid_rank
    integer :: var_grid_center_lat, var_grid_center_lon
    integer :: var_grid_corner_lat, var_grid_corner_lon
    integer :: var_grid_area, var_grid_imask, var_grid_dims

    ! helpers
    real(ESMF_KIND_R8) :: sphereRadius
    real(ESMF_KIND_R8) :: pi
    character(len=:), allocatable :: on_a_sphere
    integer :: attlen
    integer :: iCell, iVertex, lastValidVertex
    integer :: ndims_landIceMask, dimids_landIceMask(NF90_MAX_VAR_DIMS)

    doLandIceMask = .false.
    if (present(useLandIceMask)) doLandIceMask = useLandIceMask

    if (doLandIceMask) then
       write(*,'(A)') ' -- Landice Masks are enabled'
    else
       write(*,'(A)') ' -- Landice Masks are disabled'
    end if
    write(*,*)

    pi = acos(-1.0)

    stat = nf90_open(trim(mpasFile), nf90_nowrite, fin)
    call check_nc(stat, 'nf90_open('//trim(mpasFile)//')')

    stat = nf90_inq_dimid(fin, 'nCells', dimid)
    call check_nc(stat, 'nf90_inq_dimid(nCells)')
    stat = nf90_inquire_dimension(fin, dimid, len=nCells)
    call check_nc(stat, 'nf90_inquire_dimension(nCells)')

    stat = nf90_inq_dimid(fin, 'nVertices', dimid)
    call check_nc(stat, 'nf90_inq_dimid(nVertices)')
    stat = nf90_inquire_dimension(fin, dimid, len=nVertices)
    call check_nc(stat, 'nf90_inquire_dimension(nVertices)')

    stat = nf90_inq_dimid(fin, 'maxEdges', dimid)
    call check_nc(stat, 'nf90_inq_dimid(maxEdges)')
    stat = nf90_inquire_dimension(fin, dimid, len=maxVertices)
    call check_nc(stat, 'nf90_inquire_dimension(maxEdges)')

    allocate(latCell(nCells), lonCell(nCells))
    allocate(latVertex(nVertices), lonVertex(nVertices))
    allocate(areaCell(nCells))
    allocate(verticesOnCell(maxVertices, nCells))
    allocate(nEdgesOnCell(nCells))

    stat = nf90_inq_varid(fin, 'latCell', varid)
    call check_nc(stat, 'nf90_inq_varid(latCell)')
    stat = nf90_get_var(fin, varid, latCell)
    call check_nc(stat, 'nf90_get_var(latCell)')

    stat = nf90_inq_varid(fin, 'lonCell', varid)
    call check_nc(stat, 'nf90_inq_varid(lonCell)')
    stat = nf90_get_var(fin, varid, lonCell)
    call check_nc(stat, 'nf90_get_var(lonCell)')

    stat = nf90_inq_varid(fin, 'latVertex', varid)
    call check_nc(stat, 'nf90_inq_varid(latVertex)')
    stat = nf90_get_var(fin, varid, latVertex)
    call check_nc(stat, 'nf90_get_var(latVertex)')

    stat = nf90_inq_varid(fin, 'lonVertex', varid)
    call check_nc(stat, 'nf90_inq_varid(lonVertex)')
    stat = nf90_get_var(fin, varid, lonVertex)
    call check_nc(stat, 'nf90_get_var(lonVertex)')

    stat = nf90_inq_varid(fin, 'verticesOnCell', varid)
    call check_nc(stat, 'nf90_inq_varid(verticesOnCell)')
    stat = nf90_get_var(fin, varid, verticesOnCell)
    call check_nc(stat, 'nf90_get_var(verticesOnCell)')

    stat = nf90_inq_varid(fin, 'nEdgesOnCell', varid)
    call check_nc(stat, 'nf90_inq_varid(nEdgesOnCell)')
    stat = nf90_get_var(fin, varid, nEdgesOnCell)
    call check_nc(stat, 'nf90_get_var(nEdgesOnCell)')

    stat = nf90_inq_varid(fin, 'areaCell', varid)
    call check_nc(stat, 'nf90_inq_varid(areaCell)')
    stat = nf90_get_var(fin, varid, areaCell)
    call check_nc(stat, 'nf90_get_var(areaCell)')

    stat = nf90_get_att(fin, nf90_global, 'sphere_radius', sphereRadius)
    call check_nc(stat, 'nf90_get_att(sphere_radius)')

    stat = nf90_inquire_attribute(fin, nf90_global, 'on_a_sphere', len=attlen)
    call check_nc(stat, 'nf90_inquire_attribute(on_a_sphere)')
    allocate(character(len=attlen) :: on_a_sphere)
    stat = nf90_get_att(fin, nf90_global, 'on_a_sphere', on_a_sphere)
    call check_nc(stat, 'nf90_get_att(on_a_sphere)')

    if (any(lonCell < 0.0 .or. lonCell > 2.0*pi)) then
       error stop 'lonCell is not in the desired range (0, 2pi)'
    end if

    if (any(lonVertex < 0.0 .or. lonVertex > 2.0*pi)) then
       error stop 'lonVertex is not in the desired range (0, 2pi)'
    end if

    if (sphereRadius <= 0.0) then
       sphereRadius = SHR_CONST_REARTH
       write(*,'(A,ES24.16)') ' -- WARNING: sphereRadius<=0 so setting sphereRadius = ', &
            SHR_CONST_REARTH
    end if

    if (trim(on_a_sphere) == 'NO') then
       write(*,'(A)') " -- WARNING: 'on_a_sphere' attribute is 'NO', which means there may be some disagreement regarding area between the planar (source) and spherical (target) mesh"
    end if

    if (doLandIceMask) then
       stat = nf90_inq_varid(fin, 'landIceMask', varid)
       call check_nc(stat, 'nf90_inq_varid(landIceMask)')

       stat = nf90_inquire_variable(fin, varid, ndims=ndims_landIceMask, dimids=dimids_landIceMask)
       call check_nc(stat, 'nf90_inquire_variable(landIceMask)')

       if (ndims_landIceMask == 1) then
          allocate(landIceMask1d(nCells))
          stat = nf90_get_var(fin, varid, landIceMask1d)
          call check_nc(stat, 'nf90_get_var(landIceMask 1D)')
       else if (ndims_landIceMask == 2) then
          allocate(landIceMask2d(1, nCells))
          stat = nf90_get_var(fin, varid, landIceMask2d, start=(/1,1/), count=(/1,nCells/))
          call check_nc(stat, 'nf90_get_var(landIceMask 2D first slice)')
       else
          error stop 'landIceMask has unsupported rank'
       end if
    end if

    ! allocate(grid_corner_lat(nCells, maxVertices))
    ! allocate(grid_corner_lon(nCells, maxVertices))
    allocate(grid_area(nCells))
    allocate(grid_imask(nCells))

    ! grid_corner_lat = 0.0
    ! grid_corner_lon = 0.0
    grid_area = areaCell / (sphereRadius*sphereRadius)
    grid_dims(1) = nCells


    allocate(grid_corner_lat(maxVertices, nCells))
    allocate(grid_corner_lon(maxVertices, nCells))

    grid_corner_lat = 0.0
    grid_corner_lon = 0.0

    do iCell = 1, nCells
       lastValidVertex = verticesOnCell(nEdgesOnCell(iCell), iCell)

       do iVertex = 1, maxVertices
          if (iVertex <= nEdgesOnCell(iCell)) then
             grid_corner_lat(iVertex, iCell) = latVertex(verticesOnCell(iVertex, iCell))
             grid_corner_lon(iVertex, iCell) = lonVertex(verticesOnCell(iVertex, iCell))
          else
             grid_corner_lat(iVertex, iCell) = latVertex(lastValidVertex)
             grid_corner_lon(iVertex, iCell) = lonVertex(lastValidVertex)
          end if
       end do
    end do


    ! do iCell = 1, nCells
    !   lastValidVertex = verticesOnCell(nEdgesOnCell(iCell), iCell)

    !   do iVertex = 1, maxVertices
    !     if (iVertex <= nEdgesOnCell(iCell)) then
    !       grid_corner_lat(iCell, iVertex) = latVertex(verticesOnCell(iVertex, iCell))
    !       grid_corner_lon(iCell, iVertex) = lonVertex(verticesOnCell(iVertex, iCell))
    !     else
    !       grid_corner_lat(iCell, iVertex) = latVertex(lastValidVertex)
    !       grid_corner_lon(iCell, iVertex) = lonVertex(lastValidVertex)
    !     end if
    !   end do
    ! end do

    if (doLandIceMask) then
       if (allocated(landIceMask1d)) then
          grid_imask = 1 - landIceMask1d
       else
          grid_imask = 1 - landIceMask2d(1,:)
       end if
    else
       grid_imask = 1
    end if

    stat = nf90_create(trim(scripFile), ior(nf90_clobber, nf90_netcdf4), fout)
    call check_nc(stat, 'nf90_create('//trim(scripFile)//', NETCDF4)')

    stat = nf90_def_dim(fout, 'grid_size',    nCells,      dim_grid_size)
    call check_nc(stat, 'nf90_def_dim(grid_size)')
    stat = nf90_def_dim(fout, 'grid_corners', maxVertices, dim_grid_corners)
    call check_nc(stat, 'nf90_def_dim(grid_corners)')
    stat = nf90_def_dim(fout, 'grid_rank',    1,           dim_grid_rank)
    call check_nc(stat, 'nf90_def_dim(grid_rank)')

    stat = nf90_def_var(fout, 'grid_center_lat', nf90_double, (/dim_grid_size/), var_grid_center_lat)
    call check_nc(stat, 'nf90_def_var(grid_center_lat)')
    stat = nf90_put_att(fout, var_grid_center_lat, 'units', 'radians')
    call check_nc(stat, 'nf90_put_att(grid_center_lat:units)')

    stat = nf90_def_var(fout, 'grid_center_lon', nf90_double, (/dim_grid_size/), var_grid_center_lon)
    call check_nc(stat, 'nf90_def_var(grid_center_lon)')
    stat = nf90_put_att(fout, var_grid_center_lon, 'units', 'radians')
    call check_nc(stat, 'nf90_put_att(grid_center_lon:units)')

    stat = nf90_def_var(fout, 'grid_corner_lat', nf90_double, &
         &(/dim_grid_corners,dim_grid_size /), var_grid_corner_lat)
    call check_nc(stat, 'nf90_def_var(grid_corner_lat)')
    stat = nf90_put_att(fout, var_grid_corner_lat, 'units', 'radians')
    call check_nc(stat, 'nf90_put_att(grid_corner_lat:units)')

    stat = nf90_def_var(fout, 'grid_corner_lon', nf90_double, &
         &(/dim_grid_corners,dim_grid_size/), var_grid_corner_lon)
    call check_nc(stat, 'nf90_def_var(grid_corner_lon)')
    stat = nf90_put_att(fout, var_grid_corner_lon, 'units', 'radians')
    call check_nc(stat, 'nf90_put_att(grid_corner_lon:units)')

    stat = nf90_def_var(fout, 'grid_area', nf90_double, (/dim_grid_size/), var_grid_area)
    call check_nc(stat, 'nf90_def_var(grid_area)')
    stat = nf90_put_att(fout, var_grid_area, 'units', 'radian^2')
    call check_nc(stat, 'nf90_put_att(grid_area:units)')

    stat = nf90_def_var(fout, 'grid_imask', nf90_int, (/dim_grid_size/), var_grid_imask)
    call check_nc(stat, 'nf90_def_var(grid_imask)')
    stat = nf90_put_att(fout, var_grid_imask, 'units', 'unitless')
    call check_nc(stat, 'nf90_put_att(grid_imask:units)')

    stat = nf90_def_var(fout, 'grid_dims', nf90_int, (/dim_grid_rank/), var_grid_dims)
    call check_nc(stat, 'nf90_def_var(grid_dims)')

    stat = nf90_enddef(fout)
    call check_nc(stat, 'nf90_enddef')

    stat = nf90_put_var(fout, var_grid_center_lat, latCell)
    call check_nc(stat, 'nf90_put_var(grid_center_lat)')
    stat = nf90_put_var(fout, var_grid_center_lon, lonCell)
    call check_nc(stat, 'nf90_put_var(grid_center_lon)')
    stat = nf90_put_var(fout, var_grid_corner_lat, grid_corner_lat)
    call check_nc(stat, 'nf90_put_var(grid_corner_lat)')
    stat = nf90_put_var(fout, var_grid_corner_lon, grid_corner_lon)
    call check_nc(stat, 'nf90_put_var(grid_corner_lon)')
    stat = nf90_put_var(fout, var_grid_area, grid_area)
    call check_nc(stat, 'nf90_put_var(grid_area)')
    stat = nf90_put_var(fout, var_grid_imask, grid_imask)
    call check_nc(stat, 'nf90_put_var(grid_imask)')
    stat = nf90_put_var(fout, var_grid_dims, grid_dims)
    call check_nc(stat, 'nf90_put_var(grid_dims)')

    stat = nf90_close(fin)
    call check_nc(stat, 'nf90_close(input)')
    stat = nf90_close(fout)
    call check_nc(stat, 'nf90_close(output)')

  contains

    subroutine check_nc(status, where)
      integer, intent(in) :: status
      character(len=*), intent(in) :: where
      if (status /= nf90_noerr) then
         write(*,'(A)') 'NetCDF error in '//trim(where)//': '//trim(nf90_strerror(status))
         error stop 1
      end if
    end subroutine check_nc

  end subroutine mpas_to_scrip_mesh

!-----------------------------------------------------------------------------

end module
