!  Program Name: module_sf_noahmpdrv
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

module module_sf_noahmpdrv

   ! -- 1st type, to be used by the weather disaggregate, for module_sf_noahmpdrv
   implicit none
   
   !integer, parameter    :: ims = 1
   !integer, parameter    :: ime = 3
   !integer, parameter    :: jms = 1
   !integer, parameter    :: jme = 3   
   !integer, parameter    :: nsoil = 2
   
   type :: NoahMP_vectorizing_tools
      integer :: ims,ime, jms,jme, nsoil
      real,dimension( :,:), allocatable      ::  TMN   ! deep soil temperature [K]
      real,dimension( :,:), allocatable      ::  TSK   ! surface radiative temperature [K]
      real,dimension( :,:,:), allocatable    ::  TSLB  ! soil temperature [K]   

      contains
      procedure :: init => NoahMP_vectorizing_tools_init
      procedure :: destroy => NoahMP_vectorizing_tools_destroy  
      ! for debugging/testing
      procedure :: noahMP_inputs => NoahMP_vectorizing_tools_variables
      
   end type NoahMP_vectorizing_tools
      
   contains 
      
   subroutine NoahMP_vectorizing_tools_init(this, ims,ime, jms,jme, TMN, TSK, TSLB, nsoil)
      implicit none
      class(NoahMP_vectorizing_tools), intent(inout) :: this  ! the type object being initialized
      integer, intent(in) :: ims,ime, jms,jme, nsoil
      real,dimension( ims:ime,jms:jme), intent(IN)         ::  TMN
      real,dimension( ims:ime,jms:jme), intent(IN)         ::  TSK
      real,dimension( ims:ime,1:nsoil,jms:jme), intent(IN) ::  TSLB
      this%ims = ims
      this%ime = ime
      this%jms = jms
      this%jme = jme
      this%TMN = TMN
      this%TSK = TSK
      this%TSLB = TSLB
      this%nsoil = nsoil
   end subroutine NoahMP_vectorizing_tools_init

   subroutine NoahMP_vectorizing_tools_destroy(this)
      implicit none
      ! the type object being destroyed
      class(NoahMP_vectorizing_tools), intent(inout) :: this 
   end subroutine NoahMP_vectorizing_tools_destroy

   subroutine NoahMP_vectorizing_tools_variables(this)
      implicit none
      class(NoahMP_vectorizing_tools), intent(inout) :: this       
      ! real, intent(out), dimension(ims:ime,jms:jme) ::  TMN, TSK, TSLB
      ! integer, intent(out) :: nsoil
      this.TMN =reshape ((/293.,293.,293., &
                           293.,293.,293., &
                           293.,293.,293./), shape(this.TMN))   !- 20C
      this.TSK =reshape ((/296.,296.,296., &
                           296.,296.,296., &
                           296.,296.,296./), shape(this.TSK))   !- 23C
      this.TSLB =reshape ((/294.,294.,294., &
                            294.,294.,294., &
                            294.,294.,294., &
                            294.,294.,294., &
                            294.,294.,294., &
                            294.,294.,294./), shape(this.TSLB))   !- 21C
      !nsoil = 2
      print *, "DONE; module_sf_noahmpdrv.f90 ; NoahMP "   
   end subroutine NoahMP_vectorizing_tools_variables

end module module_sf_noahmpdrv

