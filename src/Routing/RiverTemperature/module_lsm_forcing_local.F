!  Program Name: module_lsm_forcing
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

module module_lsm_forcing

   implicit none

   type :: weather_2D_inputs
      integer :: IX, JX
      real, dimension(:,:), allocatable   :: T2     ! Air temperature
      real, dimension(:,:), allocatable   :: q2x    ! Humidity
      real, dimension(:,:), allocatable   :: u      ! Near surface wind in the u-component
      real, dimension(:,:), allocatable   :: short  ! Shortwave radiation.

      contains
      procedure :: init => weather_2D_inputs_init
      procedure :: destroy => weather_2D_inputs_destroy
                
      ! For debugging               
      procedure :: test => test_function
      procedure :: inputs => weather_2D_inputs_variables
   end type weather_2D_inputs

   contains

   subroutine test_function(this)
      class(weather_2D_inputs), intent(inout) :: this
      print *, "in test_function"
      !print *, "this.IX is ", this.IX
      !print *, "this.JX is ", this.JX
      this.short = reshape ((/500.,500.,500., &
                              500.,500.,500., &
                              500.,500.,500./), shape(this.short))
   end subroutine

   subroutine weather_2D_inputs_init(this, IX, JX, T2, q2x, u, short)
      implicit none
      class(weather_2D_inputs), intent(inout) :: this  ! the type object being initialized
      integer, intent(in)                     :: IX,JX !
      real, intent(in), dimension(:,:)        :: T2    ! Air temperature
      real, intent(in), dimension(:,:)        :: q2x   ! Humidity
      real, intent(in), dimension(:,:)        :: u     ! Near surface wind in the u-component
      real, intent(in), dimension(:,:)        :: short ! Shortwave radiation
      ! Assign the values passed in to a particular
      this%IX        = IX
      this%JX        = JX
      this%T2        = T2
      this%q2x       = q2x
      this%u         = u
      this%short     = short
   end subroutine weather_2D_inputs_init

   subroutine weather_2D_inputs_destroy (this)
      implicit none
      ! the type object being destroyed
      class(weather_2D_inputs), intent(inout) :: this
   end subroutine weather_2D_inputs_destroy

   subroutine weather_2D_inputs_variables(this)
      implicit none
      class(weather_2D_inputs), intent(inout) :: this   ! the type object being initialized
      !real, intent(out), dimension(IX,JX)    :: T2     ! Air temperature
      !real, intent(out), dimension(IX,JX)    :: q2x    ! Humidity
      !real, intent(out), dimension(IX,JX)    :: u      ! Near surface wind in the u-component
      !real, intent(out), dimension(IX,JX)    :: short  ! Shortwave radiation.
      this.T2 = reshape ((/18.,18.,18., &
                           18.,18.,18., &
                           18.,18.,18./), shape(this.T2))
      this.q2x = reshape ((/85.,85.,85., &
                            85.,85.,85., &
                            85.,85.,85./), shape(this.q2x))
      this.u = reshape ((/1.,1.,1., &
                          1.,1.,1., &
                          1.,1.,1./), shape(this.u))
      this.short = reshape ((/500.,500.,500., &
                              500.,500.,500., &
                              500.,500.,500./), shape(this.short))
   end subroutine weather_2D_inputs_variables

end module module_lsm_forcing

