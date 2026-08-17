#include "WRFHydro_NUOPC_Macros.h"

module wrfhydro_nuopc_fields
! !MODULE: wrfhydro_nuopc_fields
!
! !DESCRIPTION:
!   This module connects NUOPC field information for WRFHYDRO
!
! !REVISION HISTORY:
!  21Jul23    Dan Rosen  Initial Specification
!
! !USES:
  use ESMF
  use NUOPC
  use WRFHydro_ESMF_Extensions
  use WRFHydro_NUOPC_Flags
  use wrfhydro_nuopc_utils, only: printa, check
  use config_base,      only: nlst
  use module_rt_data,   only: rt_domain
  use overland_data,    only: overland_struct
  use overland_control, only: overland_control_struct

  implicit none

  private

  logical, parameter :: debug_esmf = .false.

  type cap_fld_type
    sequence
    character(len=64)      :: sd_name   = "dummy" ! standard name
    character(len=64)      :: st_name   = "dummy" ! state name
    character(len=64)      :: units     = "-"     ! units
    type(ESMF_RegridMethod_Flag) :: regrid_method = ESMF_REGRIDMETHOD_BILINEAR
    logical                :: ad_import = .false. ! advertise import
    logical                :: ad_export = .false. ! advertise export
    real(ESMF_KIND_R8)     :: vl_fillv  = ESMF_MISSING_VALUE ! default
    integer                :: num_dims  = 2
    logical                :: rl_import = .false. ! realize import
    logical                :: rl_export = .false. ! realize export
    type(ESMF_RouteHandle) :: import_handle
    logical                :: import_handle_init = .false.
    type(ESMF_RouteHandle) :: export_handle
    logical                :: export_handle_init = .false.
    type(ESMF_Field)       :: import_field
    logical                :: import_field_init = .false.
    type(ESMF_Field)       :: export_field
    logical                :: export_field_init = .false.
  end type cap_fld_type

  character(len=ESMF_MAXSTR), parameter :: file = __FILE__
  logical, parameter :: IMPORT_T = .true.
  logical, parameter :: IMPORT_F = .false.
  logical, parameter :: EXPORT_T = .true.
  logical, parameter :: EXPORT_F = .false.
  logical, parameter :: TMP_EXPORT_T = .false.
  logical, parameter :: TMP_IMPORT_T = .false.

  logical, parameter :: EXPORT_SF_HEAD = .true.
  ! logical, parameter :: EXPORT_SF_HEAD = .false.

  logical, parameter :: EXPORT_SMC = .true.
  ! logical, parameter :: EXPORT_SMC = .false.

  logical, parameter :: EXPORT_SH20 = .true.
  ! logical, parameter :: EXPORT_SH20 = .false.


  type(cap_fld_type),target,dimension(21) :: cap_fld_list = (/          &
    !cap_fld_type("inst_total_soil_moisture_ratio","smcratio", &
    !             "-", ESMF_REGRIDMETHOD_CONSERVE, &
    !             IMPORT_F, EXPORT_F, 1.0d0, 3),         &
    !cap_fld_type("inst_total_soil_moisture_content_returned","smcret", &
    !             "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
    !             IMPORT_F, EXPORT_F, 0.20d0, 3),         &
    !cap_fld_type("inst_soil_moisture_content","slc", &
    !             "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
    !             TMP_IMPORT_T, TMP_EXPORT_T, 0.20d0, 3),         &
    !cap_fld_type("inst_soil_temperature","stc", &
    !             "K     ", ESMF_REGRIDMETHOD_CONSERVE, &
    !             TMP_IMPORT_T, EXPORT_F, 288.d0, 3),             &
    cap_fld_type("liquid_fraction_of_soil_moisture_layer_1","sh2ox1", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_SH20, 0.20d0),         &
    cap_fld_type("liquid_fraction_of_soil_moisture_layer_2","sh2ox2", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_SH20, 0.20d0),         &
    cap_fld_type("liquid_fraction_of_soil_moisture_layer_3","sh2ox3", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_SH20, 0.20d0),         &
    cap_fld_type("liquid_fraction_of_soil_moisture_layer_4","sh2ox4", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_SH20, 0.20d0),         &
    ! MPAS -> hydro leg: import-only, carries the current/previous smois
    ! ratio computed on the MPAS side.
    cap_fld_type("soil_moisture_ratio_layer_1","smcratio1", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "-", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 1.0d0),         &
    cap_fld_type("soil_moisture_ratio_layer_2","smcratio2", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "-", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 1.0d0),         &
    cap_fld_type("soil_moisture_ratio_layer_3","smcratio3", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "-", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 1.0d0),         &
    cap_fld_type("soil_moisture_ratio_layer_4","smcratio4", &
                 ! "m3 m-3", ESMF_REGRIDMETHOD_BILINEAR, &
                 "-", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 1.0d0),         &
    ! hydro -> MPAS leg: export-only, carries WRF-Hydro's actual (real, not
    ! ratio) soil moisture. Distinct name/standard name from smc1-4 so the
    ! two directions never share a field or its backing memory.
    cap_fld_type("soil_moisture_fraction_layer_1_returned","smcret1", &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_F, EXPORT_SMC, 0.20d0),         &
    cap_fld_type("soil_moisture_fraction_layer_2_returned","smcret2", &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_F, EXPORT_SMC, 0.20d0),         &
    cap_fld_type("soil_moisture_fraction_layer_3_returned","smcret3", &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_F, EXPORT_SMC, 0.20d0),         &
    cap_fld_type("soil_moisture_fraction_layer_4_returned","smcret4", &
                 "m3 m-3", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_F, EXPORT_SMC, 0.20d0),         &
    cap_fld_type("soil_temperature_layer_1","stc1", &
                 "K", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 288.d0),                 &
    cap_fld_type("soil_temperature_layer_2","stc2", &
                 "K     ", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 288.d0),             &
    cap_fld_type("soil_temperature_layer_3","stc3", &
                 "K     ", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 288.d0),             &
    cap_fld_type("soil_temperature_layer_4","stc4", &
                 "K     ", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_T, EXPORT_F, 288.d0),             &
    cap_fld_type("soil_porosity","smcmax1", &
                 "1     ", ESMF_REGRIDMETHOD_CONSERVE, &
                 IMPORT_F, EXPORT_F, 0.45d0),                 &
    cap_fld_type("vegetation_type","vegtyp", &
                 ! "1     ", ESMF_REGRIDMETHOD_NEAREST_DTOS, &
                 "1     ", ESMF_REGRIDMETHOD_NEAREST_STOD, &
                 IMPORT_F, EXPORT_F, 16.0d0),                 & !FIX
    cap_fld_type("surface_water_depth","sfchead", &
                 "mm    ", ESMF_REGRIDMETHOD_CONSERVE, &
                 ! "mm    ", ESMF_REGRIDMETHOD_BILINEAR, &
                 IMPORT_F, EXPORT_SF_HEAD, 0.00d0),             &
    cap_fld_type("time_step_infiltration_excess","infxsrt", &
                 "mm    ", ESMF_REGRIDMETHOD_CONSERVE, &
                 ! "mm    ", ESMF_REGRIDMETHOD_BILINEAR, &
                 IMPORT_T, EXPORT_F, 0.00d0),             &
    cap_fld_type("soil_column_drainage","soldrain", &
                 "mm    ", ESMF_REGRIDMETHOD_CONSERVE, &
                 ! "mm    ", ESMF_REGRIDMETHOD_BILINEAR, &
                 IMPORT_T, EXPORT_F, 0.00d0)             &
    ! FOOBAR: double check this is ok
    ! these two accumulated variables break during runtime
    ! it could be they are pointing to the same variable
    ! as infxsrt and soldrain on the MPAS side
    !cap_fld_type("surface_runoff_accumulated","sfcrunoff", &
    !             "mm    ", ESMF_REGRIDMETHOD_CONSERVE, &
    !             IMPORT_F, EXPORT_F, 0.00d0),             &
    !cap_fld_type("subsurface_runoff_accumulated","udrunoff", &
    !             "mm    ", ESMF_REGRIDMETHOD_CONSERVE, &
    !             ! IMPORT_F, EXPORT_SMC, 0.00d0)              &
    !             IMPORT_F, EXPORT_F, 0.00d0)              &
    /)

  public cap_fld_list
  public field_dictionary_add
  public field_create
  public field_realize
  public field_advertise
  public check_lsm_forcings
  public field_advertise_log
  public field_realize_log
  public read_impexp_config_flnm
  public field_find_standardname
  public field_find_statename
  public state_fill_uniform
  public state_fill_prescribe
  public state_fill_file
  public state_copy_tohyd
  public state_copy_frhyd
  public state_check_missing
  public state_prescribe_missing
  public model_debug

  interface field_create
     module procedure field_create_grid
     module procedure field_create_mesh
  end interface field_create

  interface field_realize
     module procedure field_realize_grid
     module procedure field_realize_mesh
  end interface field_realize

  character(*), parameter :: filename = "WRFHydro_NUOPC_Fields.F90"

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------

  subroutine field_dictionary_add(fieldList, rc)
    type(cap_fld_type), intent(in) :: fieldList(:)
    integer, intent(out)           :: rc
    ! local variables
    integer :: n
    logical :: isPresent
    character(*), parameter :: method ="field_dictionary_add"


    rc = ESMF_SUCCESS

    do n=lbound(fieldList,1),ubound(fieldList,1)
      isPresent = NUOPC_FieldDictionaryHasEntry( &
        fieldList(n)%sd_name, rc=rc)
      call check(rc, __LINE__, file)
      if (.not.isPresent) then
        call NUOPC_FieldDictionaryAddEntry( &
          StandardName=trim(fieldList(n)%sd_name), &
          canonicalUnits=trim(fieldList(n)%units), &
          rc=rc)
        call check(rc, __LINE__, file)
      end if
    end do

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine field_realize_mesh(fieldList, importState, exportState, mesh, &
  did, realizeAllImport, realizeAllExport, memr_import, memr_export, rc)
    type(cap_fld_type), intent(inout) :: fieldList(:)
    type(ESMF_State), intent(inout)   :: importState
    type(ESMF_State), intent(inout)   :: exportState
    type(ESMF_Mesh), intent(in)       :: mesh
    integer, intent(in)               :: did
    logical, intent(in)               :: realizeAllImport
    logical, intent(in)               :: realizeAllExport
    type(memory_flag)                 :: memr_import
    type(memory_flag)                 :: memr_export
    integer, intent(out)              :: rc
    ! local variables
    integer :: n
    logical :: realizeImport
    logical :: realizeExport
    type(ESMF_Field) :: field_import
    type(ESMF_Field) :: field_export
    type(ESMF_Field) :: field
    character(*), parameter :: method = "field_realize_mpas"

    rc = ESMF_SUCCESS

   if (debug_esmf) call ESMF_LogWrite("WRFH: enter field_realize_mesh", &
               ESMF_LOGMSG_INFO, rc=rc)
    call check(rc, __LINE__, file)

    ! call ESMF_LogWrite("WRFH: realize mesh, validating states", &
    !            ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_StateValidate(importState, rc=rc)
    call check(rc, __LINE__, file)
    if (rc /= ESMF_SUCCESS) then
       error stop "WRFH: import state not valid"
    end if
    if (debug_esmf) call ESMF_LogWrite("WRFH: PRINTING IMPORTSTATE: ", &
            ESMF_LOGMSG_INFO, rc=rc)
    if (debug_esmf) call ESMF_StateLog(importState, logMsgFlag=ESMF_LOGMSG_INFO, rc=rc)
    ! call ESMF_StateLog(importState, logMsgFlag=ESMF_LOGMSG_INFO, rc=rc)
    ! call ESMF_StateValidate(exportState, rc=rc)
    ! call check(rc, __LINE__)
    ! if (rc /= ESMF_SUCCESS) then
    !    error stop "WRFH: export state not valid"
    ! end if
    ! call ESMF_StateLog(exportState, logMsgFlag=ESMF_LOGMSG_INFO, rc=rc)
    ! stop "Validating states"

    do n=lbound(fieldList,1),ubound(fieldList,1)
       ! fieldList(n)%ad_import = .true.
       ! print *, "WRFH: after statelog: ",trim(fieldList(n)%st_name)," and ",fieldList(n)%ad_import
       ! call ESMF_LogWrite("WRFH: after statelog: "//trim(fieldList(n)%st_name), &
       !      ESMF_LOGMSG_INFO, rc=rc)

      ! this looks good
      ! print *, trim(fieldList(n)%st_name), ": ad_import=", fieldList(n)%ad_import
      ! stop "debug"
      ! check realize import
      if (fieldList(n)%ad_import) then
        if (realizeAllImport) then
          realizeImport = .true.
       else
          realizeImport = NUOPC_IsConnected(importState, &
            fieldName=trim(fieldList(n)%st_name), rc=rc)
          call check(rc, __LINE__, file)
          if (realizeImport .eqv. .true. ) then
             call printa("realize_import True, name = "// &
                  fieldList(n)%st_name// "|")
          end if
          if (realizeImport .eqv. .false. ) then
             call printa("realize_import False, name = "// &
                  fieldList(n)%st_name// "|")
          end if

          ! call ESMF_StateLog(importState, logMsgFlag=ESMF_LOGMSG_INFO, rc=rc)
          ! call check(rc, __LINE__, file)

          ! call ESMF_Stateget(importState, fieldList(n)%st_name, field, rc=rc)
          ! call check(rc, __LINE__, file)

          ! call ESMF_FieldValidate(field, rc=rc)
          ! call check(rc, __LINE__, file)

          ! call ESMF_FieldWrite(field, "foo.txt", rc=rc)
          ! call check(rc, __LINE__, file)
          ! call ESMF_FieldGet(f, name=nm, rank=rank,  rc=rc)
          ! call check(rc, __LINE__)

          ! call dump_sig(fieldList(n)%st_name, field, importState)
          ! call model_debug(importState, did, memr_import, 'wrfh_debug_', rc&
          !      &=rc)

          ! call check(rc, __LINE__, file)

          ! stop "DEBUG THIS"
        end if
      else
        realizeImport = .false.
      end if
      ! create import field
      if ( realizeImport ) then
        field_import=field_create(fld_name=fieldList(n)%st_name, &
             mesh=mesh, did=did, memflg=memr_import, rc=rc)
        call check(rc, __LINE__, file)
        call NUOPC_Realize(importState, field=field_import, rc=rc)
        call check(rc, __LINE__, file)
        if (debug_esmf) call printa("realize import "// trim(fieldList(n)%st_name))
        fieldList(n)%rl_import = .true.
     else
        if (debug_esmf) call printa("realize remove import "//&
             trim(fieldList(n)%st_name)// &
             " ---------------------------")
        call ESMF_StateRemove(importState, (/fieldList(n)%st_name/), &
             relaxedflag=.true., rc=rc)
        call check(rc, __LINE__, file)
        fieldList(n)%rl_import = .false.
      end if

      ! check realize export
      if (fieldList(n)%ad_export) then
        if (realizeAllExport) then
          realizeExport = .true.
       else
          ! call ESMF_LogWrite("WRFH: printing export state realize", &
          !      ESMF_LOGMSG_INFO, rc=rc)
          ! call ESMF_StateLog(exportState, logMsgFlag=ESMF_LOGMSG_INFO, rc=rc)
          realizeExport = NUOPC_IsConnected(exportState, &
               fieldName=trim(fieldList(n)%st_name), rc=rc)
          call check(rc, __LINE__, file)
          if (realizeExport .eqv. .true.) then
             if (debug_esmf) call printa("export, true realizeexport of " &
                  // trim(fieldList(n)%st_name))
          else
             if (debug_esmf) call printa("export, false realizeexport of " &
                  // trim(fieldList(n)%st_name))
          end if
          ! call ESMF_LogWrite("WRFH: this realize needs to be true", &
          !      ESMF_LOGMSG_INFO, rc=rc)
        end if
      else
          realizeExport = .false.
      end if
      ! create export field
      if( realizeExport ) then
        field_export = field_create(fld_name=fieldList(n)%st_name, &
             mesh=mesh, did=did, memflg=memr_export, rc=rc)
        ! stop "investigating export, one should stop"
        call check(rc, __LINE__, file)
        call NUOPC_Realize(exportState, field=field_export, rc=rc)
        call check(rc, __LINE__, file)
        fieldList(n)%rl_export = .true.
        ! print *, "WRFH: created export ", trim(fieldList(n)%st_name)
        if (debug_esmf) call printa("realize export "// trim(fieldList(n)%st_name))
      else
        if (debug_esmf) call printa("realize remove export "//&
             trim(fieldList(n)%st_name)// &
             " ---------------------------")
        call ESMF_StateRemove(exportState, (/fieldList(n)%st_name/), &
             relaxedflag=.true., rc=rc)
        call check(rc, __LINE__, file)
        fieldList(n)%rl_export = .false.
      end if
   end do
   if (debug_esmf) call printa("exit field_realize_mesh")
   ! stop "debugging field_realize_mesh"
  end subroutine field_realize_mesh


  subroutine field_realize_grid(fieldList, importState, exportState, grid, &
  did, realizeAllImport, realizeAllExport, memr_import, memr_export, rc)
    type(cap_fld_type), intent(inout) :: fieldList(:)
    type(ESMF_State), intent(inout)   :: importState
    type(ESMF_State), intent(inout)   :: exportState
    type(ESMF_Grid), intent(in)       :: grid
    integer, intent(in)               :: did
    logical, intent(in)               :: realizeAllImport
    logical, intent(in)               :: realizeAllExport
    type(memory_flag)                 :: memr_import
    type(memory_flag)                 :: memr_export
    integer, intent(out)              :: rc
    ! local variables
    integer :: n
    logical :: realizeImport
    logical :: realizeExport
    type(ESMF_Field) :: field_import
    type(ESMF_Field) :: field_export
    character(*), parameter :: method = "field_realize"

    rc = ESMF_SUCCESS

    do n=lbound(fieldList,1),ubound(fieldList,1)
      ! check realize import
      if (fieldList(n)%ad_import) then
        if (realizeAllImport) then
          realizeImport = .true.
        else
          realizeImport = NUOPC_IsConnected(importState, &
            fieldName=trim(fieldList(n)%st_name), rc=rc)
          call check(rc, __LINE__, file)
        end if
      else
        realizeImport = .false.
      end if
      ! create import field
      if ( realizeImport ) then
        field_import=field_create(fld_name=fieldList(n)%st_name, &
          grid=grid, did=did, memflg=memr_import, rc=rc)
        call check(rc, __LINE__, file)
        call NUOPC_Realize(importState, field=field_import, rc=rc)
        call check(rc, __LINE__, file)
        fieldList(n)%rl_import = .true.
      else
        call ESMF_StateRemove(importState, (/fieldList(n)%st_name/), &
          relaxedflag=.true., rc=rc)
        call check(rc, __LINE__, file)
        fieldList(n)%rl_import = .false.
      end if

      ! check realize export
      if (fieldList(n)%ad_export) then
        if (realizeAllExport) then
          realizeExport = .true.
        else
          realizeExport = NUOPC_IsConnected(exportState, &
            fieldName=trim(fieldList(n)%st_name), rc=rc)
          call check(rc, __LINE__, file)
        end if
      else
        realizeExport = .false.
      end if
      ! create export field
      if( realizeExport ) then
        field_export=field_create(fld_name=fieldList(n)%st_name, &
          grid=grid, did=did, memflg=memr_export, rc=rc)
        call check(rc, __LINE__, file)
        call NUOPC_Realize(exportState, field=field_export, rc=rc)
        call check(rc, __LINE__, file)
        fieldList(n)%rl_export = .true.
      else
        call ESMF_StateRemove(exportState, (/fieldList(n)%st_name/), &
          relaxedflag=.true., rc=rc)
        call check(rc, __LINE__, file)
        fieldList(n)%rl_export = .false.
      end if
    end do
  end subroutine

  !-----------------------------------------------------------------------------

  function check_lsm_forcings(importState,rc)
    ! RETURN
    logical :: check_lsm_forcings
    ! ARGUMENTS
    type(ESMF_State), intent(in) :: importState
    integer, intent(out)         :: rc
    ! LOCAL VARIABLES
    integer                    :: fieldIndex
    type(ESMF_StateItem_Flag)  :: itemType
    integer                    :: s_smc, s_smc1, s_smc2, s_smc3, s_smc4
    integer                    :: s_slc, s_slc1, s_slc2, s_slc3, s_slc4
    integer                    :: s_stc, s_stc1, s_stc2, s_stc3, s_stc4
    integer                    :: s_infxsrt
    integer                    :: s_soldrain
    logical                    :: c_smc
    logical                    :: c_slc
    logical                    :: c_stc
    logical                    :: c_infxsrt
    logical                    :: c_soldrain
    character(*), parameter :: method = "lsm_forcings"

    ! total soil moisture content
    !call ESMF_StateGet(importState,itemSearch="smcratio", itemCount=s_smc, rc=rc)
    !call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="smcratio1",itemCount=s_smc1, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="smcratio2",itemCount=s_smc2, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="smcratio3",itemCount=s_smc3, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="smcratio4",itemCount=s_smc4, rc=rc)
    call check(rc, __LINE__, file)
    !if (s_smc.gt.0) then
    !  c_smc = NUOPC_IsConnected(importState, fieldName="smcratio")
    !elseif ((s_smc1.gt.0) .and. (s_smc2.gt.0) .and. &
    !        (s_smc3.gt.0) .and. (s_smc4.gt.0)) then
    if ((s_smc1.gt.0) .and. (s_smc2.gt.0) .and. &
            (s_smc3.gt.0) .and. (s_smc4.gt.0)) then
      c_smc = (NUOPC_IsConnected(importState, fieldName="smcratio1") .and. &
               NUOPC_IsConnected(importState, fieldName="smcratio2") .and. &
               NUOPC_IsConnected(importState, fieldName="smcratio3") .and. &
               NUOPC_IsConnected(importState, fieldName="smcratio4"))
    else
      c_smc = .false.
    end if

    ! liquid soil moisture content
    !call ESMF_StateGet(importState,itemSearch="slc", itemCount=s_slc, rc=rc)
    !call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="sh2ox1",itemCount=s_slc1, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="sh2ox2",itemCount=s_slc2, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="sh2ox3",itemCount=s_slc3, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="sh2ox4",itemCount=s_slc4, rc=rc)
    call check(rc, __LINE__, file)
    !if (s_slc.gt.0) then
    !  c_slc = NUOPC_IsConnected(importState, fieldName="slc")
    !elseif ((s_slc1.gt.0) .and. (s_slc2.gt.0) .and. &
    !        (s_slc3.gt.0) .and. (s_slc4.gt.0)) then
    if ((s_slc1.gt.0) .and. (s_slc2.gt.0) .and. &
            (s_slc3.gt.0) .and. (s_slc4.gt.0)) then
      c_slc = (NUOPC_IsConnected(importState, fieldName="sh2ox1") .and. &
               NUOPC_IsConnected(importState, fieldName="sh2ox2") .and. &
               NUOPC_IsConnected(importState, fieldName="sh2ox3") .and. &
               NUOPC_IsConnected(importState, fieldName="sh2ox4"))
    else
      c_slc = .false.
    end if

    ! soil temperature
    !call ESMF_StateGet(importState,itemSearch="stc", itemCount=s_stc, rc=rc)
    !call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="stc1",itemCount=s_stc1, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="stc2",itemCount=s_stc2, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="stc3",itemCount=s_stc3, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_StateGet(importState,itemSearch="stc4",itemCount=s_stc4, rc=rc)
    call check(rc, __LINE__, file)
    !if (s_stc.gt.0) then
    !  c_stc = NUOPC_IsConnected(importState, fieldName="stc")
    !elseif ((s_stc1.gt.0) .and. (s_stc2.gt.0) .and. &
    !        (s_stc3.gt.0) .and. (s_stc4.gt.0)) then
    if ((s_stc1.gt.0) .and. (s_stc2.gt.0) .and. &
            (s_stc3.gt.0) .and. (s_stc4.gt.0)) then
      c_stc = (NUOPC_IsConnected(importState, fieldName="stc1") .and. &
               NUOPC_IsConnected(importState, fieldName="stc2") .and. &
               NUOPC_IsConnected(importState, fieldName="stc3") .and. &
               NUOPC_IsConnected(importState, fieldName="stc4"))
    else
      c_stc = .false.
    end if

    ! infiltration excess
    call ESMF_StateGet(importState,itemSearch="infxsrt",itemCount=s_infxsrt, rc=rc)
    call check(rc, __LINE__, file)
    if (s_infxsrt.gt.0) then
      c_infxsrt = NUOPC_IsConnected(importState, fieldName="infxsrt")
    else
      c_infxsrt = .false.
    end if

    ! soil drainage
    call ESMF_StateGet(importState,itemSearch="soldrain",itemCount=s_soldrain, rc=rc)
    call check(rc, __LINE__, file)
    if (s_soldrain.gt.0) then
      c_soldrain = NUOPC_IsConnected(importState, fieldName="soldrain")
    else
      c_soldrain = .false.
    end if

    check_lsm_forcings = c_smc .and. c_slc .and. c_stc .and. &
                         c_infxsrt .and. c_soldrain

    ! print *, "WRFH NUOPC TODO: fix check_lsm_forcings"
    ! ! check_lsm_forcings = .false.
    ! check_lsm_forcings = .true.
  end function

  !-----------------------------------------------------------------------------

  subroutine field_advertise(fieldList, importState, exportState, &
  transferOffer, rc)
    type(cap_fld_type), intent(in)    :: fieldList(:)
    type(ESMF_State), intent(inout)   :: importState
    type(ESMF_State), intent(inout)   :: exportState
    character(*), intent(in),optional :: transferOffer
    integer, intent(out)              :: rc
    ! local variables
    integer :: n
    character(*), parameter :: method = "field_advertise"

    rc = ESMF_SUCCESS

    do n=lbound(fieldList,1),ubound(fieldList,1)
      if (fieldList(n)%ad_import) then
        call NUOPC_Advertise(importState, &
          StandardName=fieldList(n)%sd_name, &
          Units=fieldList(n)%units, &
          TransferOfferGeomObject=transferOffer, &
          name=fieldList(n)%st_name, &
          rc=rc)
        call check(rc, __LINE__, file)
        if (debug_esmf) then
           print *, "WRFH: advertise import ", trim(fieldList(n)%sd_name), &
             trim(fieldList(n)%st_name), &
             trim(fieldList(n)%units), &
             trim(transferOffer)
        end if

      end if
      if (fieldList(n)%ad_export) then
        call NUOPC_Advertise(exportState, &
          StandardName=fieldList(n)%sd_name, &
          Units=fieldList(n)%units, &
          TransferOfferGeomObject=transferOffer, &
          name=fieldList(n)%st_name, &
          rc=rc)
        call check(rc, __LINE__, file)
        if (debug_esmf) print *, "WRFH: advertise export ", fieldList(n)%st_name
      end if
    end do

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine field_advertise_log(fieldList, cname, rc)
    type(cap_fld_type), intent(in) :: fieldList(:)
    character(*), intent(in)       :: cname
    integer, intent(out)           :: rc
    ! local variables
    integer                    :: cntImp
    integer                    :: cntExp
    integer                    :: n
    character(32)              :: label
    character(ESMF_MAXSTR)     :: logMsg
    character(*), parameter :: method = "field_advertise_log"

    rc = ESMF_SUCCESS

    label = trim(cname)

    ! count advertised import and export fields
    cntImp = 0
    cntExp = 0
    do n = lbound(fieldList,1), ubound(fieldList,1)
      if (fieldList(n)%ad_import) cntImp = cntImp + 1
      if (fieldList(n)%ad_export) cntExp = cntExp + 1
    enddo

    ! log advertised import fields
    write(logMsg,'(a,a,i0,a)') trim(label)//': ', &
      'List of advertised import fields(',cntImp,'):'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') trim(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    cntImp = 0
    do n=lbound(fieldList,1), ubound(fieldList,1)
      if (.NOT.fieldList(n)%ad_import) cycle
      cntImp = cntImp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') trim(label)//': ', &
        cntImp,' ',trim(fieldList(n)%st_name), &
        ' ',trim(fieldList(n)%sd_name)
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    enddo

    ! log advertised export fields
    write(logMsg,'(a,a,i0,a)') trim(label)//': ', &
      'List of advertised export fields(',cntExp,'):'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') trim(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    cntExp = 0
    do n=lbound(fieldList,1), ubound(fieldList,1)
      if (.NOT.fieldList(n)%ad_export) cycle
      cntExp = cntExp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') trim(label)//': ', &
        cntExp,' ',trim(fieldList(n)%st_name), &
        ' ',trim(fieldList(n)%sd_name)
      call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    enddo

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine field_realize_log(fieldList, cname, rc)
    type(cap_fld_type), intent(in) :: fieldList(:)
    character(*), intent(in)       :: cname
    integer, intent(out)           :: rc
    ! local variables
    integer                    :: cntImp
    integer                    :: cntExp
    integer                    :: n
    character(32)              :: label
    character(ESMF_MAXSTR)     :: logMsg
    character(*), parameter :: method = "field_realize_log"

    rc = ESMF_SUCCESS

    label = trim(cname)

    ! count realized import and export fields
    cntImp = 0
    cntExp = 0
    do n = lbound(fieldList,1), ubound(fieldList,1)
      if (fieldList(n)%rl_import) cntImp = cntImp + 1
      if (fieldList(n)%rl_export) cntExp = cntExp + 1
    enddo

    ! log realized import fields
    write(logMsg,'(a,a,i0,a)') trim(label)//': ', &
      'List of realized import fields(',cntImp,'):'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') trim(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    cntImp = 0
    do n=lbound(fieldList,1), ubound(fieldList,1)
      if (.NOT.fieldList(n)%rl_import) cycle
      cntImp = cntImp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') trim(label)//': ', &
        cntImp,' ',trim(fieldList(n)%st_name), &
        ' ',trim(fieldList(n)%sd_name)
      call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
    enddo

    ! log realized export fields
    write(logMsg,'(a,a,i0,a)') trim(label)//': ', &
      'List of realized export fields(',cntExp,'):'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    write(logMsg,'(a,a5,a,a16,a,a)') trim(label)//': ', &
      'index',' ','name',' ','standardName'
    call ESMF_LogWrite(trim(logMsg), ESMF_LOGMSG_INFO)
    cntExp = 0
    do n=lbound(fieldList,1), ubound(fieldList,1)
      if (.NOT.fieldList(n)%rl_export) cycle
      cntExp = cntExp + 1
      write(logMsg,'(a,i5,a,a16,a,a)') trim(label)//': ', &
        cntExp,' ',trim(fieldList(n)%st_name), &
        ' ',trim(fieldList(n)%sd_name)
      call ESMF_LogWrite(trim(LogMsg), ESMF_LOGMSG_INFO)
    enddo

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine read_impexp_config_flnm(fname, fieldList, rc)
    character(len=30),intent(in)     :: fname
    type(cap_fld_type),intent(inout) :: fieldList(:)
    integer,intent(out)              :: rc

    ! local variables
    type(ESMF_Config)                  :: fieldsConfig
    type(NUOPC_FreeFormat)             :: attrFF
    integer                            :: lineCount
    integer                            :: tokenCount
    character(len=NUOPC_FreeFormatLen),allocatable :: tokenList(:)
    integer                            :: i,j
    integer                            :: stat
    character(*), parameter :: method = "read_impexp_config_flnm"

    rc = ESMF_SUCCESS

!   load fname into fieldsConfig
    fieldsConfig = ESMF_ConfigCreate(rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_ConfigLoadFile(fieldsConfig, fname, rc=rc)
    call check(rc, __LINE__, file)

!   read export fields from config
    attrFF = NUOPC_FreeFormatCreate(fieldsConfig, &
      label="hyd_fields", rc=rc)
    call check(rc, __LINE__, file)
    call NUOPC_FreeFormatGet(attrFF, lineCount=lineCount, rc=rc)
    call check(rc, __LINE__, file)
    do i=1, lineCount
      call NUOPC_FreeFormatGetLine(attrFF, line=i, &
        tokenCount=tokenCount, rc=rc)
      call check(rc, __LINE__, file)
      if (.not.((tokenCount.eq.5).or.(tokenCount.eq.6))) then
        call ESMF_LogSetError(ESMF_FAILURE, &
          msg="Malformed ocn_export_fields item FORMAT="// &
            "'STATE_NAME' 'STANDARD_NAME' 'UNITS' 'IMPORT' 'EXPORT' "// &
!            "['FILLVAL'] "// &
            "in file: "//trim(fname), &
          CONTEXT, rcToReturn=rc)
        return
      end if
      allocate(tokenList(tokenCount))
      call NUOPC_FreeFormatGetLine(attrFF, line=i, tokenList=tokenList, rc=rc)
      call check(rc, __LINE__, file)
      call field_find_statename(fieldList, tokenList(1), location=j, rc=rc)
      call check(rc, __LINE__, file)
      fieldList(j)%st_name=tokenList(1)
      fieldList(j)%sd_name=tokenList(2)
      fieldList(j)%units=tokenList(3)
      tokenList(4) = ESMF_UtilStringUpperCase(tokenList(4), rc=rc)
      call check(rc, __LINE__, file)
      fieldList(j)%ad_import=((tokenList(4).eq.".TRUE.") .or. &
                         (tokenList(4).eq."TRUE"))
      print*, " int(fieldList(j)%ad_import)=", fieldList(j)%ad_import
      error stop "ADIMPORT?? DEBUGGING THIS TO SEE WHY AD_IMPORT WOULD BE FALSE?"
      tokenList(5) = ESMF_UtilStringUpperCase(tokenList(5), rc=rc)
      call check(rc, __LINE__, file)
      fieldList(j)%ad_export=((tokenList(5).eq.".TRUE.") .or. &
                         (tokenList(5).eq."TRUE"))
      if (tokenCount.eq.6) then
        fieldList(j)%vl_fillv = ESMF_UtilString2Real(tokenList(6), rc=rc)
        call check(rc, __LINE__, file)
      end if
      deallocate(tokenList)
    enddo

!   cleanup
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    call check(rc, __LINE__, file)
    call ESMF_ConfigDestroy(fieldsConfig, rc=rc)
    call check(rc, __LINE__, file)

  end subroutine read_impexp_config_flnm

  !-----------------------------------------------------------------------------

  subroutine field_find_standardname(fieldList, standardName, location, &
  fillValue, rc)
    type(cap_fld_type), intent(in)          :: fieldList(:)
    character(len=64), intent(in)           :: standardName
    integer, intent(out), optional          :: location
    real(ESMF_KIND_R8),intent(out),optional :: fillValue
    integer, intent(out)                    :: rc
    ! local variables
    integer :: n
    character(*), parameter :: method = "field_find_standardname"

    rc = ESMF_RC_NOT_FOUND

    if (present(location)) location = lbound(fieldList,1) - 1
    if (present(fillValue)) fillValue = ESMF_MISSING_VALUE

    do n=lbound(fieldList,1),ubound(fieldList,1)
      if (fieldList(n)%sd_name .eq. standardName) then
        if (present(location)) location = n
        if (present(fillValue)) fillValue = fieldList(n)%vl_fillv
        rc = ESMF_SUCCESS
        return
      end if
    end do

    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="Field not found in fieldList "//trim(standardName), &
      line=__LINE__, &
      file=__FILE__)) &
      return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine field_find_statename(fieldList, stateName, location, &
  fillValue, rc)
    type(cap_fld_type), intent(in)          :: fieldList(:)
    character(len=64), intent(in)           :: stateName
    integer, intent(out), optional          :: location
    real(ESMF_KIND_R8),intent(out),optional :: fillValue
    integer, intent(out)                    :: rc
    ! local variables
    integer :: n
    character(*), parameter :: method = "field_find_statename"

    rc = ESMF_RC_NOT_FOUND

    if (present(location)) location = lbound(fieldList,1) - 1
    if (present(fillValue)) fillValue = ESMF_MISSING_VALUE

    do n=lbound(fieldList,1),ubound(fieldList,1)
      if (fieldList(n)%st_name .eq. stateName) then
        if (present(location)) location = n
        if (present(fillValue)) fillValue = fieldList(n)%vl_fillv
        rc = ESMF_SUCCESS
        return
      end if
    end do

    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="Field not found in fieldList "//trim(stateName), &
      line=__LINE__, &
      file=__FILE__)) &
      return

  end subroutine field_find_statename

  !-----------------------------------------------------------------------------
  function field_create_mesh(fld_name, mesh, did, memflg, rc) &
       result(field_create)
    ! return value
    type(ESMF_Field) :: field_create
    ! arguments
    character(*), intent(in)      :: fld_name
    type(ESMF_Mesh), intent(in)   :: mesh
    integer, intent(in)           :: did
    type(memory_flag), intent(in) :: memflg
    integer,          intent(out) :: rc
    ! local variables
    character(len=16)       :: cmemflg

    real(ESMF_KIND_R8), pointer :: p(:)   ! mesh element field is 1D local array


    rc = ESMF_SUCCESS
    ! print *, "=== WRFH: realize, try field create mesh: ", trim(fld_name)

    if (memflg .eq. MEMORY_POINTER) then
       ! field_create = ESMF_FieldCreate(name=fld_name, mesh=mesh, &
       !      typekind=ESMF_TYPEKIND_R8, &
       !      indexflag=ESMF_INDEX_DELOCAL, &
       !      meshloc=ESMF_MESHLOC_ELEMENT, &
       !      rc=rc)
       field_create = ESMF_FieldCreate(name=fld_name, mesh=mesh, &
            typekind=ESMF_TYPEKIND_R8, &
            indexflag=ESMF_INDEX_DELOCAL, &
            meshloc=ESMF_MESHLOC_ELEMENT, &
            rc=rc)
       call check(rc, __LINE__, file)
       call ESMF_FieldGet(field_create, farrayPtr=p, rc=rc)
       call check(rc, __LINE__, file)
       ! p = -9999.0d0
       p = 0.0

       ! print *, "=== WRFH: realize, success field created mesh: ", trim(fld_name)

       ! error stop "hi"

      ! select case (trim(fld_name))
        ! case ('smc')
        !   field_create = ESMF_FieldCreate(name=fld_name, mesh=mesh, &
        !        typekind=ESMF_TYPEKIND_R8, &
        !        indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !        ! farray=rt_domain(did)%smc(:,:,:), &
        !   ! field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !   !   farray=rt_domain(did)%smc(:,:,:), gridToFieldMap=(/1,2/), &
        !   !   ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
        !   !   indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('slc')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%sh2ox(:,:,:), gridToFieldMap=(/1,2/), &
        !     ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('stc')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%stc(:,:,:), gridToFieldMap=(/1,2/), &
        !     ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('sh2ox1')
        !   ! field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !   !   farray=rt_domain(did)%sh2ox(:,:,1), &
        !   !   indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   field_create = ESMF_FieldCreate(name=fld_name, mesh=mesh, &
        !       typekind=ESMF_TYPEKIND_R8, &
        !       indexflag=ESMF_INDEX_DELOCAL, &
        !       meshloc=ESMF_MESHLOC_ELEMENT, &
        !       rc=rc)
        !  call check(rc, __LINE__, file)

        !   call check(rc, __LINE__, file)
        ! case ('sh2ox2')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%sh2ox(:,:,2), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('sh2ox3')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%sh2ox(:,:,3), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('sh2ox4')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%sh2ox(:,:,4), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
      ! case ('smc1')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%smc(:,:,1), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('smc2')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%smc(:,:,2), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('smc3')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%smc(:,:,3), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('smc4')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%smc(:,:,4), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('smcmax1')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%smcmax1, &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
      ! case ('stc1')
          ! field_create = ESMF_FieldCreate(name=fld_name, mesh=mesh, &
          !      meshloc=ESMF_MESHLOC_ELEMENT, &
          !   ! farray=rt_domain(did)%stc(:,:,1), &
          !      indexflag=ESMF_INDEX_DELOCAL, rc=rc)
         ! field_create = ESMF_FieldCreate(name=fld_name, mesh=mesh, &
         !      typekind=ESMF_TYPEKIND_R8, &
         !      indexflag=ESMF_INDEX_DELOCAL, &
         !      meshloc=ESMF_MESHLOC_ELEMENT, &
         !      rc=rc)
         ! call check(rc, __LINE__, file)
      ! case ('stc2')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%stc(:,:,2), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('stc3')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%stc(:,:,3), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('stc4')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%stc(:,:,4), &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('vegtyp')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%vegtyp, &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('sfchead')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%overland%control%surface_water_head_lsm, &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('infxsrt')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%infxsrt, &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        ! case ('soldrain')
        !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !     farray=rt_domain(did)%soldrain, &
        !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !   call check(rc, __LINE__, file)
        !  case default
        !    stop "ADD BACK IN " // trim(fld_name)
        !    !     call ESMF_LogSetError(ESMF_FAILURE, &
        !       msg=method//": Field hookup missing: "//trim(fld_name), &
        !       file=filename, rcToReturn=rc)
        !   return
    !   end select
    ! elseif (memflg .eq. MEMORY_COPY) then
    !   select case (trim(fld_name))
    !     case ('smc','slc','stc')
    !       field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
    !         typekind=ESMF_TYPEKIND_FIELD, gridToFieldMap=(/1,2/), &
    !         ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
    !         indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    !       call check(rc, __LINE__, file)
    !     case default
    !       field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
    !         typekind=ESMF_TYPEKIND_FIELD, &
    !         indexflag=ESMF_INDEX_DELOCAL, rc=rc)
    !       call check(rc, __LINE__, file)
    !   end select
    !   call ESMF_FieldFill(field_create, dataFillScheme="const", &
    !     const1=ESMF_MISSING_VALUE, rc=rc)
    !   call check(rc, __LINE__, file)
    ! else
    !   cmemflg = memflg
    !   call ESMF_LogSetError(ESMF_FAILURE, &
    !     msg=method//": Field memory flag unknown: "//trim(cmemflg), &
    !     file=filename, rcToReturn=rc)
    !   return
    end if
    ! stop "CHECKING"
  end function


  function field_create_grid(fld_name,grid,did,memflg,rc) &
       result(field_create)
    use module_NoahMP_hrldas_driver, only: sfcrunoff, udrunoff
    ! return value
    type(ESMF_Field) :: field_create
    ! arguments
    character(*), intent(in)      :: fld_name
    type(ESMF_Grid), intent(in)   :: grid
    integer, intent(in)           :: did
    type(memory_flag), intent(in) :: memflg
    integer,          intent(out) :: rc
    ! local variables
    character(len=16)       :: cmemflg
    character(*), parameter :: method = "field_create->field_create_grid"

    rc = ESMF_SUCCESS

    ! print *, "filed_create_grid fld_name = ", trim(fld_name)
    if (memflg .eq. MEMORY_POINTER) then
      select case (trim(fld_name))
        !case ('smcratio')
        !  field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !    farray=rt_domain(did)%smc(:,:,:)*0.0+1.0, gridToFieldMap=(/1,2/), &
        !    ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
        !    indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !  call check(rc, __LINE__, file)
        !case ('smcret')
        !  field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !    farray=rt_domain(did)%smc(:,:,:), gridToFieldMap=(/1,2/), &
        !    ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
        !    indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !  call check(rc, __LINE__, file)
        !case ('slc')
        !  field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !    farray=rt_domain(did)%sh2ox(:,:,:), gridToFieldMap=(/1,2/), &
        !    ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
        !    indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !  call check(rc, __LINE__, file)
        !case ('stc')
        !  field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
        !    farray=rt_domain(did)%stc(:,:,:), gridToFieldMap=(/1,2/), &
        !    ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
        !    indexflag=ESMF_INDEX_DELOCAL, rc=rc)
        !  call check(rc, __LINE__, file)
        case ('sh2ox1')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%sh2ox(:,:,1), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('sh2ox2')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%sh2ox(:,:,2), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('sh2ox3')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%sh2ox(:,:,3), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('sh2ox4')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%sh2ox(:,:,4), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcratio1')
           ! print *, "rt_domain(did)%smc(1:2,1:2,1) = ", rt_domain(did)%smc(1:2,1:2,1)

           rt_domain(did)%smcratio1_buf = 1.0
           field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smcratio1_buf, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
           call check(rc, __LINE__, file)
        case ('smcratio2')
          rt_domain(did)%smcratio2_buf = 1.0
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smcratio2_buf, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcratio3')
          rt_domain(did)%smcratio3_buf = 1.0
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smcratio3_buf, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcratio4')
          rt_domain(did)%smcratio4_buf = 1.0
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smcratio4_buf, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcret1')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smc(:,:,1), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcret2')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smc(:,:,2), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcret3')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smc(:,:,3), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcret4')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smc(:,:,4), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('smcmax1')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%smcmax1, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('stc1')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%stc(:,:,1), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('stc2')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%stc(:,:,2), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('stc3')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%stc(:,:,3), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('stc4')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%stc(:,:,4), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('vegtyp')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%vegtyp, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
       case ('sfchead')
          ! print *, "sfchead=",&
          !      rt_domain(did)%overland%control%surface_water_head_lsm(0:5,0:5)
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%overland%control%surface_water_head_lsm, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
          ! stop "CHECKING"
        case ('infxsrt')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%infxsrt, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case ('soldrain')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            farray=rt_domain(did)%soldrain, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)

       ! WRFH TODO: these are writing toe the same vars as infxsrt/soldrain
          ! add for update the WRF state variable.
          ! THIS IS FROM module_wrf_HYDRO.F90:209
            ! do k = 1, nlst(did)%nsoil
            !     ! grid%TSLB(its:ite,k,jts:jte) = RT_DOMAIN(did)%STC(:,:,k)
            !     grid%smois(its:ite,k,jts:jte) = RT_DOMAIN(did)%smc(:,:,k)
            !     grid%sh2o(its:ite,k,jts:jte) = RT_DOMAIN(did)%sh2ox(:,:,k)
            ! end do
          ! update WRF variable after running routing model.
            ! grid%sfcheadrt(its:ite,jts:jte) = &
            !   rt_domain(did)%overland%control%surface_water_head_lsm

       !case ('sfcrunoff')
       !
       !   ! rt_domain(did)%overland%control%surface_water_head_lsm(:,:) = &
       !   !      -776
       !
       !   ! 102 x 60
       !   ! print *, shape(rt_domain(did)%overland%control%surface_water_head_lsm(:,:))
       !   stop "this should not be accessed"
       !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
       !     farray=rt_domain(did)%overland%control%surface_water_head_lsm(:,:), &
       !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
       !   call check(rc, __LINE__, file)
       !case ('udrunoff')
       !   field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
       !     ! farray=rt_domain(did)%udrunoff(:,:), &
       !     farray=rt_domain(did)%overland%control%surface_water_head_lsm(:,:), &
       !     indexflag=ESMF_INDEX_DELOCAL, rc=rc)
       !   call check(rc, __LINE__, file)
       case default
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg=method//": Field hookup missing: "//trim(fld_name), &
            file=filename, rcToReturn=rc)
          return
      end select
    elseif (memflg .eq. MEMORY_COPY) then
      select case (trim(fld_name))
        case ('smcratio','smcret','slc','stc')
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            typekind=ESMF_TYPEKIND_FIELD, gridToFieldMap=(/1,2/), &
            ungriddedLBound=(/1/), ungriddedUBound=(/nlst(did)%nsoil/), &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
        case default
          field_create = ESMF_FieldCreate(name=fld_name, grid=grid, &
            typekind=ESMF_TYPEKIND_FIELD, &
            indexflag=ESMF_INDEX_DELOCAL, rc=rc)
          call check(rc, __LINE__, file)
      end select
      call ESMF_FieldFill(field_create, dataFillScheme="const", &
        const1=ESMF_MISSING_VALUE, rc=rc)
      call check(rc, __LINE__, file)
    else
      cmemflg = memflg
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg=method//": Field memory flag unknown: "//trim(cmemflg), &
        file=filename, rcToReturn=rc)
      return
    end if
    ! print *, "end field_create_grid rc =", rc

  end function

  !-----------------------------------------------------------------------------

  subroutine state_fill_uniform(state, fillValue, rc)
    type(ESMF_State), intent(inout)        :: state
    real(ESMF_KIND_R8), intent(in)         :: fillValue
    integer, intent(out)                   :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: stat
    character(*), parameter :: method = "state_fill_uniform"

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item name memory failed.", &
      line=__LINE__, file=__FILE__)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item type memory failed.", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    call check(rc, __LINE__, file)

    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state, field=field, &
          itemName=itemNameList(n), rc=rc)
        call check(rc, __LINE__, file)
        call ESMF_FieldFill(field, dataFillScheme="const", &
          const1=fillValue, rc=rc)
        call check(rc, __LINE__, file)
        call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
        call check(rc, __LINE__, file)
      end if
    enddo

    deallocate(itemNameList)
    deallocate(itemTypeList)

  end subroutine state_fill_uniform

  !-----------------------------------------------------------------------------

  subroutine state_fill_prescribe(state, fieldList, rc)
    type(ESMF_State), intent(inout)        :: state
    type(cap_fld_type), intent(in)         :: fieldList(:)
    integer, intent(out)                   :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    real(ESMF_KIND_R8)                     :: filVal
    integer                                :: stat
    character(*), parameter :: method = "state_fill_prescribe"

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item name memory failed.", &
      line=__LINE__, file=__FILE__)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item type memory failed.", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    call check(rc, __LINE__, file)

    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        call field_find_statename(fieldList, &
          stateName=itemNameList(n), fillValue=filVal, rc=rc)
        call check(rc, __LINE__, file)
        call ESMF_StateGet(state, field=field, &
          itemName=itemNameList(n), rc=rc)
        call check(rc, __LINE__, file)
        call ESMF_FieldFill(field, dataFillScheme="const", &
          const1=filVal, rc=rc)
        call check(rc, __LINE__, file)
        call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
        call check(rc, __LINE__, file)
      end if
    enddo

    deallocate(itemNameList)
    deallocate(itemTypeList)

  end subroutine state_fill_prescribe

 !-----------------------------------------------------------------------------

  subroutine state_fill_file(state, filePrefix, rc)
    type(ESMF_State), intent(inout) :: state
    character(len=*), intent(in)    :: filePrefix
    integer, intent(out)            :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    character(len=64)                      :: fldName
    integer                                :: stat
    character(*), parameter :: method = "state_fill_file"

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item name memory failed.", &
      line=__LINE__, file=__FILE__)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item type memory failed.", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    call check(rc, __LINE__, file)

    do n=1, itemCount
      if ( itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state,field=field, &
          itemName=itemNameList(n), rc=rc)
        call check(rc, __LINE__, file)
        call NUOPC_GetAttribute(field, name="StandardName", &
          value=fldName, rc=rc)
        call check(rc, __LINE__, file)
        call ESMF_FieldRead(field, variableName=trim(fldName), &
          fileName=trim(filePrefix)//"_"//trim(itemNameList(n))//".nc", &
          iofmt=ESMF_IOFMT_NETCDF, rc=rc)
        call check(rc, __LINE__, file)
        call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
        call check(rc, __LINE__, file)
      end if
    enddo

    deallocate(itemNameList)
    deallocate(itemTypeList)

  end subroutine state_fill_file

  !-----------------------------------------------------------------------------

  subroutine state_copy_tohyd(state, did, rc)
    type(ESMF_State), intent(inout) :: state
    integer, intent(in)             :: did
    integer, intent(out)            :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: dimCount
    real(ESMF_KIND_FIELD), pointer         :: farrayPtr2d(:,:)
    real(ESMF_KIND_FIELD), pointer         :: farrayPtr3d(:,:,:)
    integer                                :: stat
    character(*), parameter :: method = "state_copy_tohyd"

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item name memory failed.", &
      line=__LINE__, file=__FILE__)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item type memory failed.", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    call check(rc, __LINE__, file)

    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state, field=field, &
          itemName=itemNameList(n), rc=rc)
        call check(rc, __LINE__, file)
        call ESMF_FieldGet(field, dimCount=dimCount, rc=rc)
        call check(rc, __LINE__, file)
        if (dimCount .eq. 2) then
          call ESMF_FieldGet(field, farrayPtr=farrayPtr2d, rc=rc)
          call check(rc, __LINE__, file)
        elseif (dimCount .eq. 3) then
          call ESMF_FieldGet(field, farrayPtr=farrayPtr3d, rc=rc)
          call check(rc, __LINE__, file)
        else
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg=method//": dimCount is not supported.", &
            file=filename, rcToReturn=rc)
          return
        end if
        select case (ItemNameList(n))
          !case ('smcratio')
          !  rt_domain(did)%smc = rt_domain(did)%smc * farrayPtr3d
          !case ('slc')
          !  rt_domain(did)%sh2ox = farrayPtr3d
          !case ('stc')
          !  rt_domain(did)%stc = farrayPtr3d
          case ('sh2ox1')
            rt_domain(did)%sh2ox(:,:,1) = farrayPtr2d
          case ('sh2ox2')
            rt_domain(did)%sh2ox(:,:,2) = farrayPtr2d
          case ('sh2ox3')
            rt_domain(did)%sh2ox(:,:,3) = farrayPtr2d
          case ('sh2ox4')
            rt_domain(did)%sh2ox(:,:,4) = farrayPtr2d
          case ('smcratio1')
            ! Incoming field is the MPAS current/previous smois ratio; apply it
            ! to the existing hydro smc value rather than overwriting it.
            print *, "ADCHECK: smc1_pre=", rt_domain(did)%smc(50,50,1), "ratio1=", farrayPtr2d(50,50)
            rt_domain(did)%smc(:,:,1) = rt_domain(did)%smc(:,:,1) * farrayPtr2d
            print *, "ADCHECK: smc1_post=", rt_domain(did)%smc(50,50,1)
          case ('smcratio2')
            rt_domain(did)%smc(:,:,2) = rt_domain(did)%smc(:,:,2) * farrayPtr2d
          case ('smcratio3')
            rt_domain(did)%smc(:,:,3) = rt_domain(did)%smc(:,:,3) * farrayPtr2d
          case ('smcratio4')
            rt_domain(did)%smc(:,:,4) = rt_domain(did)%smc(:,:,4) * farrayPtr2d
          case ('smcmax1')
            rt_domain(did)%smcmax1 = farrayPtr2d
          case ('stc1')
            rt_domain(did)%stc(:,:,1) = farrayPtr2d
          case ('stc2')
            rt_domain(did)%stc(:,:,2) = farrayPtr2d
          case ('stc3')
            rt_domain(did)%stc(:,:,3) = farrayPtr2d
          case ('stc4')
            rt_domain(did)%stc(:,:,4) = farrayPtr2d
          case ('vegtyp')
            rt_domain(did)%vegtyp = farrayPtr2d
          case ('sfchead')
            rt_domain(did)%overland%control%surface_water_head_lsm = &
              farrayPtr2d
          case ('infxsrt')
            rt_domain(did)%infxsrt = farrayPtr2d
          case ('soldrain')
            rt_domain(did)%soldrain = farrayPtr2d
          !case ('sfcrunoff')
          !   print *, "WRFH: check state_copy_tohyd for sfcrunoff is correct"
          !   rt_domain(did)%infxsrt = farrayPtr2d
          !case ('udrunoff')
          !   print *, "WRFH: check state_copy_tohyd for udrunoff is correct"
          !  rt_domain(did)%soldrain = farrayPtr2d
          case default
            call ESMF_LogSetError(ESMF_FAILURE, &
              msg=method//": Field hookup missing: "//trim(itemNameList(n)), &
              file=filename, rcToReturn=rc)
            return
        endselect
      end if
    enddo

    deallocate(itemNameList)
    deallocate(itemTypeList)

  end subroutine state_copy_tohyd

  !-----------------------------------------------------------------------------

  subroutine state_copy_frhyd(state, did, rc)
    type(ESMF_State), intent(inout)         :: state
    integer, intent(in)                     :: did
    integer, intent(out)                    :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: field
    integer                                :: dimCount
    real(ESMF_KIND_FIELD), pointer         :: farrayPtr2d(:,:)
    real(ESMF_KIND_FIELD), pointer         :: farrayPtr3d(:,:,:)
    integer                                :: stat
    character(len=16)                      :: cmissingv_flag
    character(*), parameter :: method = "state_copy_frhyd"

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item name memory failed.", &
      line=__LINE__, file=__FILE__)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item type memory failed.", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    call check(rc, __LINE__, file)

    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(state, field=field, &
          itemName=itemNameList(n), rc=rc)
        call check(rc, __LINE__, file)
        call ESMF_FieldGet(field, dimCount=dimCount, rc=rc)
        call check(rc, __LINE__, file)
        if (dimCount .eq. 2) then
          call ESMF_FieldGet(field, farrayPtr=farrayPtr2d, rc=rc)
          call check(rc, __LINE__, file)
        elseif (dimCount .eq. 3) then
          call ESMF_FieldGet(field, farrayPtr=farrayPtr3d, rc=rc)
          call check(rc, __LINE__, file)
        else
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg=method//": dimCount is not supported.", &
            file=filename, rcToReturn=rc)
          return
        end if
        select case (ItemNameList(n))
          !case ('smcret')
          !  farrayPtr3d = rt_domain(did)%smc
          !case ('slc')
          !  farrayPtr3d = rt_domain(did)%sh2ox
          !case ('stc')
          !  farrayPtr3d = rt_domain(did)%stc
          case ('sh2ox1')
            farrayPtr2d = rt_domain(did)%sh2ox(:,:,1)
          case ('sh2ox2')
            farrayPtr2d = rt_domain(did)%sh2ox(:,:,2)
          case ('sh2ox3')
            farrayPtr2d = rt_domain(did)%sh2ox(:,:,3)
          case ('sh2ox4')
            farrayPtr2d = rt_domain(did)%sh2ox(:,:,4)
          case ('smcret1')
            farrayPtr2d = rt_domain(did)%smc(:,:,1)
            print *, "ADCHECK: smc1_export=", rt_domain(did)%smc(50,50,1)
          case ('smcret2')
            farrayPtr2d = rt_domain(did)%smc(:,:,2)
          case ('smcret3')
            farrayPtr2d = rt_domain(did)%smc(:,:,3)
          case ('smcret4')
            farrayPtr2d = rt_domain(did)%smc(:,:,4)
          case ('smcmax1')
            farrayPtr2d = rt_domain(did)%smcmax1
          case ('stc1')
            farrayPtr2d = rt_domain(did)%stc(:,:,1)
          case ('stc2')
            farrayPtr2d = rt_domain(did)%stc(:,:,2)
          case ('stc3')
            farrayPtr2d = rt_domain(did)%stc(:,:,3)
          case ('stc4')
            farrayPtr2d = rt_domain(did)%stc(:,:,4)
          case ('vegtyp')
            farrayPtr2d = rt_domain(did)%vegtyp
          case ('sfchead')
            farrayPtr2d = rt_domain(did)%overland%control%surface_water_head_lsm
          case ('infxsrt')
            farrayPtr2d = rt_domain(did)%infxsrt
          case ('soldrain')
            farrayPtr2d = rt_domain(did)%soldrain
          !case ('sfcrunoff')
          !   print *, "WRFH: check state_copy_frhyd for sfcrunoff is correct"
          !   ! farrayPtr2d = rt_domain(did)%infxsrt ! original
          !   farrayPtr2d = rt_domain(did)%overland%control%surface_water_head_lsm(:,:)
          !case ('udrunoff')
          !   print *, "WRFH: check state_copy_frhyd for udrunoff is correct"
          !   farrayPtr2d = rt_domain(did)%soldrain
          case default
            call ESMF_LogSetError(ESMF_FAILURE, &
              msg=method//": Field hookup missing: "//trim(itemNameList(n)), &
              file=filename, rcToReturn=rc)
            return
        end select
      end if
    enddo

    deallocate(itemNameList)
    deallocate(itemTypeList)

  end subroutine state_copy_frhyd

  !-----------------------------------------------------------------------------

  subroutine state_check_missing(state, did, rc)
    type(ESMF_State), intent(inout)   :: state
    integer, intent(in)               :: did
    integer, intent(out)              :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    real(ESMF_KIND_R8), parameter          :: chkVal = real(ESMF_MISSING_VALUE)
    logical                                :: missng
    integer                                :: stat
    character(*), parameter :: method = "state_check_missing"

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item name memory failed.", &
      line=__LINE__, file=__FILE__)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item type memory failed.", &
      line=__LINE__, file=__FILE__)) return
    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    call check(rc, __LINE__, file)

    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        missng = .FALSE.
        select case (ItemNameList(n))
          !case ('smcratio')
          !  missng = any(rt_domain(did)%smc.eq.chkVal)
          !case ('smcret')
          !  missng = any(rt_domain(did)%smc.eq.chkVal)
          !case ('slc')
          !  missng = any(rt_domain(did)%sh2ox.eq.chkVal)
          !case ('stc')
          !  missng = any(rt_domain(did)%stc.eq.chkVal)
          case ('sh2ox1')
            missng = any(rt_domain(did)%sh2ox(:,:,1).eq.chkVal)
          case ('sh2ox2')
            missng = any(rt_domain(did)%sh2ox(:,:,2).eq.chkVal)
          case ('sh2ox3')
            missng = any(rt_domain(did)%sh2ox(:,:,3).eq.chkVal)
          case ('sh2ox4')
            missng = any(rt_domain(did)%sh2ox(:,:,4).eq.chkVal)
          case ('smcratio1')
            missng = any(rt_domain(did)%smc(:,:,1).eq.chkVal)
          case ('smcratio2')
            missng = any(rt_domain(did)%smc(:,:,2).eq.chkVal)
          case ('smcratio3')
            missng = any(rt_domain(did)%smc(:,:,3).eq.chkVal)
          case ('smcratio4')
            missng = any(rt_domain(did)%smc(:,:,4).eq.chkVal)
          case ('smcret1')
            missng = any(rt_domain(did)%smc(:,:,1).eq.chkVal)
          case ('smcret2')
            missng = any(rt_domain(did)%smc(:,:,2).eq.chkVal)
          case ('smcret3')
            missng = any(rt_domain(did)%smc(:,:,3).eq.chkVal)
          case ('smcret4')
            missng = any(rt_domain(did)%smc(:,:,4).eq.chkVal)
          case ('smcmax1')
            missng = any(rt_domain(did)%smcmax1.eq.chkVal)
          case ('stc1')
            missng = any(rt_domain(did)%stc(:,:,1).eq.chkVal)
          case ('stc2')
            missng = any(rt_domain(did)%stc(:,:,2).eq.chkVal)
          case ('stc3')
            missng = any(rt_domain(did)%stc(:,:,3).eq.chkVal)
          case ('stc4')
            missng = any(rt_domain(did)%stc(:,:,4).eq.chkVal)
          case ('vegtyp')
            missng = any(rt_domain(did)%vegtyp.eq.chkVal)
          case ('sfchead')
            missng = any(rt_domain(did)%overland%control%surface_water_head_lsm &
              .eq.chkVal)
          case ('infxsrt')
            missng = any(rt_domain(did)%infxsrt.eq.chkVal)
          case ('soldrain')
            missng = any(rt_domain(did)%soldrain.eq.chkVal)
          case default
            call ESMF_LogSetError(ESMF_FAILURE, &
              msg=method//": Field hookup missing: "//trim(itemNameList(n)), &
              file=filename, rcToReturn=rc)
            return
        endselect
        if (missng) then
          call ESMF_LogSetError(ESMF_FAILURE, &
            msg=method//": Missing value: "//trim(itemNameList(n)), &
            file=filename, rcToReturn=rc)
          return
        end if
      end if
    enddo
    deallocate(itemNameList)
    deallocate(itemTypeList)

  end subroutine state_check_missing

  !-----------------------------------------------------------------------------

  subroutine state_prescribe_missing(state, did, fieldList, rc)
    type(ESMF_State), intent(inout)   :: state
    integer, intent(in)               :: did
    type(cap_fld_type), intent(in)    :: fieldList(:)
    integer, intent(out)              :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    real(ESMF_KIND_R8), parameter          :: chkVal = real(ESMF_MISSING_VALUE)
    real(ESMF_KIND_R8)                     :: filVal
    integer                                :: stat
    character(*), parameter :: method = "state_prescribe_missing"

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
    call check(rc, __LINE__, file)

    allocate(itemNameList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item name memory failed.", &
      line=__LINE__, file=__FILE__)) return
    allocate(itemTypeList(itemCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of state item type memory failed.", &
      line=__LINE__, file=__FILE__)) return
    call ESMF_StateGet(state,itemNameList=itemNameList, &
      itemTypeList=itemTypeList, rc=rc)
    call check(rc, __LINE__, file)

    do n=1, itemCount
      if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
        call field_find_statename(fieldList, &
          stateName=itemNameList(n), fillValue=filVal, rc=rc)
        call check(rc, __LINE__, file)
        select case (itemNameList(n))
          !case ('smcratio')
          !  where (rt_domain(did)%smc.eq.chkVal) &
          !    rt_domain(did)%smc = filVal
          !case ('smcret')
          !  where (rt_domain(did)%smc.eq.chkVal) &
          !    rt_domain(did)%smc = filVal
          !case ('slc')
          !  where (rt_domain(did)%sh2ox.eq.chkVal) &
          !    rt_domain(did)%sh2ox = filVal
          !case ('stc')
          !  where (rt_domain(did)%stc.eq.chkVal) &
          !    rt_domain(did)%stc = filVal
          case ('sh2ox1')
            where (rt_domain(did)%sh2ox(:,:,1).eq.chkVal) &
              rt_domain(did)%sh2ox(:,:,1) = filVal
          case ('sh2ox2')
            where (rt_domain(did)%sh2ox(:,:,2).eq.chkVal) &
              rt_domain(did)%sh2ox(:,:,2) = filVal
          case ('sh2ox3')
            where (rt_domain(did)%sh2ox(:,:,3).eq.chkVal) &
              rt_domain(did)%sh2ox(:,:,3) = filVal
          case ('sh2ox4')
            where (rt_domain(did)%sh2ox(:,:,4).eq.chkVal) &
              rt_domain(did)%sh2ox(:,:,4) = filVal
          case ('smcratio1')
            where (rt_domain(did)%smc(:,:,1).eq.chkVal) &
              rt_domain(did)%smc(:,:,1) = filVal
          case ('smcratio2')
            where (rt_domain(did)%smc(:,:,2).eq.chkVal) &
              rt_domain(did)%smc(:,:,2) = filVal
          case ('smcratio3')
            where (rt_domain(did)%smc(:,:,3).eq.chkVal) &
              rt_domain(did)%smc(:,:,3) = filVal
          case ('smcratio4')
            where (rt_domain(did)%smc(:,:,4).eq.chkVal) &
              rt_domain(did)%smc(:,:,4) = filVal
          case ('smcret1')
            where (rt_domain(did)%smc(:,:,1).eq.chkVal) &
              rt_domain(did)%smc(:,:,1) = filVal
          case ('smcret2')
            where (rt_domain(did)%smc(:,:,2).eq.chkVal) &
              rt_domain(did)%smc(:,:,2) = filVal
          case ('smcret3')
            where (rt_domain(did)%smc(:,:,3).eq.chkVal) &
              rt_domain(did)%smc(:,:,3) = filVal
          case ('smcret4')
            where (rt_domain(did)%smc(:,:,4).eq.chkVal) &
              rt_domain(did)%smc(:,:,4) = filVal
          case ('smcmax1')
            where (rt_domain(did)%smcmax1.eq.chkVal) &
              rt_domain(did)%smcmax1 = filVal
          case ('stc1')
            where (rt_domain(did)%stc(:,:,1).eq.chkVal) &
              rt_domain(did)%stc(:,:,1) = filVal
          case ('stc2')
            where (rt_domain(did)%stc(:,:,2).eq.chkVal) &
              rt_domain(did)%stc(:,:,2) = filVal
          case ('stc3')
            where (rt_domain(did)%stc(:,:,3).eq.chkVal) &
              rt_domain(did)%stc(:,:,3) = filVal
          case ('stc4')
            where (rt_domain(did)%stc(:,:,4).eq.chkVal) &
              rt_domain(did)%stc(:,:,4) = filVal
          case ('vegtyp')
            where (rt_domain(did)%vegtyp.eq.chkVal) &
              rt_domain(did)%vegtyp = filVal
          case ('sfchead')
            where (rt_domain(did)%overland%control%surface_water_head_lsm &
              .eq.chkVal) &
              rt_domain(did)%overland%control%surface_water_head_lsm = filVal
          case ('infxsrt')
            where (rt_domain(did)%infxsrt.eq.chkVal) &
              rt_domain(did)%infxsrt = filVal
          case ('soldrain')
            where (rt_domain(did)%soldrain.eq.chkVal) &
              rt_domain(did)%soldrain = filVal
          case default
            call ESMF_LogSetError(ESMF_FAILURE, &
              msg=method//": Field hookup missing: "//trim(itemNameList(n)), &
              file=filename, rcToReturn=rc)
            return
        end select
      end if
    enddo

    deallocate(itemNameList)
    deallocate(itemTypeList)

  end subroutine state_prescribe_missing

  !-----------------------------------------------------------------------------

  subroutine model_debug(state, did, memflg, filePrefix, rc)
    type(ESMF_State), intent(inout)   :: state
    integer, intent(in)               :: did
    type(memory_flag)                 :: memflg
    character(len=*)                  :: filePrefix
    integer, intent(out)              :: rc
    ! local variables
    integer                                :: n
    integer                                :: itemCount
    character(len=64),allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: itemTypeList(:)
    type(ESMF_Field)                       :: cpyfield, outfield
    type(ESMF_Grid)                        :: grid
    integer                                :: stat
    character(len=16)                      :: cmemflg
    character(*), parameter :: method = "model_debug"

    rc = ESMF_SUCCESS

    if (memflg .eq. MEMORY_POINTER) then
      call NUOPC_Write(state, &
        fileNamePrefix=filePrefix, overwrite=.true., &
        status=ESMF_FILESTATUS_REPLACE, timeslice=1, rc=rc)
      call check(rc, __LINE__, file)
    elseif(memflg .eq. MEMORY_COPY) then
      call ESMF_StateGet(state,itemCount=itemCount, rc=rc)
      call check(rc, __LINE__, file)

      allocate(itemNameList(itemCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Allocation of state item name memory failed.", &
        line=__LINE__, file=__FILE__)) return
      allocate(itemTypeList(itemCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Allocation of state item type memory failed.", &
        line=__LINE__, file=__FILE__)) return
      call ESMF_StateGet(state,itemNameList=itemNameList, &
        itemTypeList=itemTypeList, rc=rc)
      call check(rc, __LINE__, file)

      do n=1, itemCount
        if (itemTypeList(n) == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(state, field=cpyfield, &
            itemName=itemNameList(n), rc=rc)
          call check(rc, __LINE__, file)
          call ESMF_FieldGet(cpyfield, grid=grid, rc=rc)
          call check(rc, __LINE__, file)
          outfield = field_create(itemNameList(n), grid=grid, did=did, &
            memflg=MEMORY_POINTER, rc=rc)
          call ESMF_FieldWrite(outfield, variableName=itemNameList(n), &
            fileName=trim(filePrefix)//"_"//trim(itemNameList(n))//".nc", &
            iofmt=ESMF_IOFMT_NETCDF, rc=rc)
          call ESMF_FieldDestroy(outfield, rc=rc)
          call check(rc, __LINE__, file)
        end if
      enddo
      deallocate(itemNameList)
      deallocate(itemTypeList)
    else
      cmemflg = memflg
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg=method//": Field memory flag unknown: "//trim(cmemflg), &
        file=filename, rcToReturn=rc)
      return
    end if

  end subroutine model_debug

  !-----------------------------------------------------------------------------

  subroutine printType(type)
    type(ESMF_StateItem_Flag), intent(in) :: type

    if (type == ESMF_STATEITEM_ARRAY) then
       print *, "type: ESMF_STATEITEM_ARRAY"
    end if

    if (type == ESMF_STATEITEM_ARRAYBUNDLE) then
       print *, "type: ESMF_STATEITEM_ARRAYBUNDLE"
    end if

    if (type == ESMF_STATEITEM_FIELD ) then
       print *, "type: ESMF_STATEITEM_FIELD "
    end if

    if (type == ESMF_STATEITEM_FIELDBUNDLE) then
       print *, "type: ESMF_STATEITEM_FIELDBUNDLE"
    end if

    if (type == ESMF_STATEITEM_ROUTEHANDLE) then
       print *, "type: ESMF_STATEITEM_ROUTEHANDLE"
    end if

    if (type == ESMF_STATEITEM_STATE) then
       print *, "type: ESMF_STATEITEM_STATE"
    end if


  end subroutine printType

subroutine dump_sig(st_name, f, importState)
character(len=*), intent(in) :: st_name
type(ESMF_Field), intent(inout) :: f
type(ESMF_State), intent(inout) :: importState
integer :: rc, tk, rank, stag, mloc
character(len=ESMF_MAXSTR) :: nm, stdn, units
type(ESMF_StateItem_Flag) :: itemType


call ESMF_FieldGet(f, name=nm, rank=rank,  rc=rc)
call check(rc, __LINE__, file)
! call ESMF_FieldGet(f, name=nm, typekind=tk, rank=rank, staggerloc=stag, meshloc=mloc, rc=rc)
call ESMF_AttributeGet(f, name='standard_name', value=stdn, rc=rc)
call check(rc, __LINE__, file)
call ESMF_AttributeGet(f, name='units', value=units, rc=rc)
call check(rc, __LINE__, file)
write(*,'(A)') '--- signature ['//trim(st_name)//'] ---'
write(*,'(A,1X,A)') 'name :', trim(nm)
write(*,'(A,1X,A)') 'standard_name:', trim(stdn)
write(*,'(A,1X,A)') 'units :', trim(units)
call ESMF_Stateget(importState, st_name, itemType, rc=rc)
call check(rc, __LINE__, file)
call printType(itemType)
write(*,'(A,1X,I0)') 'rank :', rank
! if (stag /= ESMF_STAGGERLOC_UNKNOWN) write(*,'(A,1X,I0)') 'staggerloc :', stag
! if (mloc /= ESMF_MESHLOC_UNKNOWN) write(*,'(A,1X,I0)') 'meshloc :', mloc
end subroutine dump_sig

end module wrfhydro_nuopc_fields
