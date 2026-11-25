module module_subsurface_static_data
    implicit none

    type subsurface_static_interface

    integer :: ixrt
    integer :: jxrt
    integer :: nsoil
    real :: dt
    integer :: subrt_option
    integer :: subrt_gradient ! 1 = terrain slope (kinematic), 2 = water table

    contains
        procedure :: init => subsurface_static_data_init
        procedure :: destroy => subsurface_static_data_destroy
    end type subsurface_static_interface

    contains

    subroutine subsurface_static_data_init(this, ixrt, jxrt, nsoil, dt, subrt_option, subrt_gradient)
        implicit none
        class(subsurface_static_interface), intent(inout) :: this ! the type object being initialized
        integer, intent(in) :: ixrt                     ! x routing grid size
        integer, intent(in) :: jxrt                     ! y routing grid size
        integer, intent(in) :: nsoil                    ! nsoil
        real, intent(in) :: dt
        integer, intent(in) :: subrt_option
        integer, intent(in) :: subrt_gradient

        logical :: allocation_error = .false.

        this%ixrt = ixrt
        this%jxrt = jxrt
        this%nsoil = nsoil
        this%dt = dt
        this%subrt_option = subrt_option
        this%subrt_gradient = subrt_gradient

        if ( allocation_error ) &
            write(0,*) "attempt to allocate data in members of subsurface io structure&
            &that where already allocated. The allocated members where not changed"

    end subroutine subsurface_static_data_init

    subroutine subsurface_static_data_destroy(this)
        implicit none
        class(subsurface_static_interface), intent(inout) :: this ! the type object being initialized

        logical :: allocation_error = .false.

        if ( allocation_error ) &
            write(0,*) "attempt to deallocate data in members of subsurface io structure&
            &that where not allocated. The unallocated members where not changed"

    end subroutine subsurface_static_data_destroy

end module module_subsurface_static_data
