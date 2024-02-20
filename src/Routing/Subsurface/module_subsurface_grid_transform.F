module module_subsurface_grid_transform
    implicit none

    type subsurface_grid_transform_interface

        ! temp array? field capacity disaggregated to routing grid, for each soil layer
        real, allocatable, dimension(:,:,:) :: smcrefrt

        ! Soil Moisture Content -- ???
        real, allocatable, dimension(:,:,:) :: smcrt

        ! Soil Moisture Content -- ??? Porosity
        real, allocatable, dimension(:,:,:) :: smcmaxrt

        ! Soil Moisture Content -- ??? Wilting Point
        real, allocatable, dimension(:,:,:) :: smcwltrt

        ! water tracer variables
        real, allocatable, dimension(:,:,:) :: smcrt_tr  ! need descrip & units
        real, allocatable, dimension(:,:,:) :: smcrt_sub  ! need descrip & units
        real, allocatable, dimension(:,:,:) :: smcmaxrt_sub  ! need descrip & units

    contains
        procedure :: init => subsurface_grid_transform_init
        procedure :: destroy => subsurface_grid_transfrom_destroy
    end type

    contains

    subroutine subsurface_grid_transform_init(this,ix,jx,nsoil,ntracer)
        implicit none
        class(subsurface_grid_transform_interface), intent(inout) :: this ! the type object being initialized
        integer, intent(in) :: ix                     ! x grid size
        integer, intent(in) :: jx                     ! y grid size
        integer, intent(in) :: nsoil                  ! number of soil layers
        ! for water tracers
        integer, optional, intent(in) :: ntracer      ! number of tracers

        logical :: allocation_error = .false.

        ! check allocation status of smcrefrt
        if ( .not. allocated(this%smcrefrt) ) then
            allocate( this%smcrefrt(ix,jx,nsoil) )
            !no initialization on this variable?
            !this%smcrefrt = 0.0
        else
            allocation_error = .true.
        end if

        ! check to see if smcrt is allocated
        if ( .not. allocated(this%smcrt) ) then
            allocate( this%smcrt(ix,jx,nsoil) )
            this%smcrt = 0.0
        else
            allocation_error = .true.
        end if

        ! check to see if smcmaxrt is allocated
        if ( .not. allocated(this%smcmaxrt) ) then
            allocate( this%smcmaxrt(ix,jx,nsoil) )
            this%smcmaxrt = 0.0
        end if

        ! check to see if smcwltrt is allocated
        if ( .not. allocated(this%smcwltrt) ) then
            allocate( this%smcwltrt(ix,jx,nsoil) )
            this%smcwltrt = 0.0
        else
            allocation_error = .true.
        end if

        !----- water tracer variables

        ! check to see if smcrt_tr is allocated
        if ( .not. allocated(this%smcrt_tr) ) then
            allocate( this%smcrt_tr(ix,jx,ntracer) )
            this%smcrt_tr = 0.0
        else
            allocation_error = .true.
        end if

        ! check to see if smcrt_sub is allocated
        if ( .not. allocated(this%smcrt_sub) ) then
            allocate( this%smcrt_sub(ix,jx,ntracer) )
            this%smcrt_sub = 0.0
        else
            allocation_error = .true.
        end if

        ! check to see if smcmaxrt_sub is allocated
        if ( .not. allocated(this%smcmaxrt_sub) ) then
            allocate( this%smcmaxrt_sub(ix,jx,ntracer) )
            this%smcmaxrt_sub = 0.0
        end if

        !----- end water tracer variables

        if ( allocation_error ) &
            write(0,*) "attempt to allocate data in members of subsurface grid transform structure&
            &that where already allocated. The allocated members where not changed"

    end subroutine subsurface_grid_transform_init

    subroutine subsurface_grid_transfrom_destroy(this)
        implicit none
        class(subsurface_grid_transform_interface), intent(inout) :: this ! the type object being initialized

        logical :: allocation_error = .false.

        ! check allocation status of smcrefrt
        if ( allocated(this%smcrefrt) ) then
            deallocate( this%smcrefrt )
        else
            allocation_error = .true.
        end if

        ! check to see if smcrt is allocated
        if ( allocated(this%smcrt) ) then
            deallocate( this%smcrt )
        else
            allocation_error = .true.
        end if

        ! check to see if smcmaxrt is allocated
        if ( allocated(this%smcmaxrt) ) then
            deallocate( this%smcmaxrt )
        else
            allocation_error = .true.
        end if

        ! check to see if smcwltrt is allocated
        if ( allocated(this%smcwltrt) ) then
            deallocate( this%smcwltrt )
        else
            allocation_error = .true.
        end if

!----- water tracer variables

        ! check to see if smcrt_tr is allocated
        if ( allocated(this%smcrt_tr) ) then
            deallocate( this%smcrt_tr )
        else
            allocation_error = .true.
        end if

        ! check to see if smcrt_sub is allocated
        if ( allocated(this%smcrt_sub) ) then
            deallocate( this%smcrt_sub )
        else
            allocation_error = .true.
        end if

        ! check to see if smcmaxrt_sub is allocated
        if ( allocated(this%smcmaxrt_sub) ) then
            deallocate( this%smcmaxrt_sub )
        else
            allocation_error = .true.
        end if

!----- end water tracer variables

        if ( allocation_error ) &
            write(0,*) "attempt to deallocate data in members of subsurface grid transform structure&
            &that where not allocated. The unallocated members where not changed"

    end subroutine

end module module_subsurface_grid_transform
