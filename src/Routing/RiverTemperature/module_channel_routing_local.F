!  Program Name: module_river_temperature main function
!  Author(s)/Contact(s): Reza Abdi (rabdi@ucar.edu)
!  Abstract:
!  History Log: First attept, Apr. 21
! --------------------------------------------------------------------------

module module_channel_routing

   implicit none

   !integer, parameter            :: IXRT = 6
   !integer, parameter            :: JXRT = 6
   !integer, parameter            :: NLINKS = 12

   type :: temperature_sim_H_inputs

      integer                                  :: NLINKS, IXRT,JXRT
      real, dimension(:), allocatable          :: ChanCondConst
      real, dimension(:,:),allocatable         :: QSTRMVOLRT        ! Overland flow
      real, dimension(:,:), allocatable        :: QINFLOWBASE       ! Groundwater flow
      real, dimension(:), allocatable          :: QSUM              ! Mass bal of node, flow from upstream
      real, dimension(:), allocatable          :: XS_Peremeter      ! wetted perimeter (m)
      real, dimension(:), allocatable          :: Bw                ! Bottom width (m)
      real, dimension(:), allocatable          :: Tw                ! top width before bankfull (m)
      !real, dimension(:),allocatable          :: velocity          ! Velocity (m/s)
      real, dimension(:), allocatable          :: velocity
      real, dimension(:), allocatable          :: HLINK             ! Water head (m)
      real, dimension(:), allocatable          :: ChSSlp            ! Channel slope (m/m)
      integer, dimension(:), allocatable       :: CHANXI, CHANYJ
      real, dimension(:), allocatable          :: ZELEV             ! elevation of nodes
      integer, dimension(:), allocatable       :: CHANLEN           ! channel lngth (m)
      integer, dimension(:,:), allocatable     :: CH_NETLNK         ! River channel IDs
      real, dimension(:,:), allocatable        :: QLINK
      integer                                  :: DT_STEPS          ! number of timestep in routing
      real                                     :: DTCT
      real, dimension(:),allocatable           :: node_area
      integer, dimension(:),allocatable        :: TO_NODE
      integer, dimension(:),allocatable        :: TYPEL
      
      
      contains
      procedure :: init => temperature_sim_H_inputs_init
      procedure :: destroy => temperature_sim_H_inputs_destroy
      
      ! For debugging and testing
      procedure :: inputs_hyd => vectors_and_variables_Hyd

   end type temperature_sim_H_inputs
        
   ! --------------------------------------------------------------------------
      
   type :: vectorizing_tools

      integer :: NLINKS, IXRT,JXRT
      integer, dimension(:,:), allocatable    :: CH_NETRT, LAKE_MSKRT
      integer, dimension(:), allocatable      :: CHANXI, CHANYJ
      
      contains
      procedure :: init => vectorizing_tools_init
      procedure :: destroy => vectorizing_tools_destroy 

      ! For debugging and testing
      procedure :: inputs_vec => vectors_and_variables_Vec 
      
   end type vectorizing_tools

   ! --------------------------------------------------------------------------

   contains

   subroutine temperature_sim_H_inputs_init(this, NLINKS, IXRT, JXRT, ChanCondConst, QSTRMVOLRT,  &
                  QINFLOWBASE, QSUM,  XS_Peremeter, Bw, Tw, velocity, HLINK, ChSSlp, &
                  CHANXI, CHANYJ, ZELEV, CHANLEN, CH_NETLNK, QLINK, DT_STEPS, DTCT, node_area, TO_NODE, TYPEL)

      implicit none
      class(temperature_sim_H_inputs), intent(inout) :: this             ! the type object being initialized
      integer, intent(in) :: NLINKS, IXRT,JXRT
      real, intent(in),  dimension(NLINKS)             :: ChanCondConst
      real,intent(in),   dimension(IXRT,JXRT)          :: QSTRMVOLRT     ! to be called like: QSTRMVOLRT(CHANXI(i),CHANYJ(i)); where i:1 to NLINKS
      real, intent(in),  dimension(IXRT,JXRT)          :: QINFLOWBASE   
      real, intent(in),  dimension(NLINKS)             :: QSUM              
      real, intent(in),  dimension(NLINKS)             :: XS_Peremeter   
      real, intent(in),  dimension(NLINKS)             :: Bw         
      real, intent(in),  dimension(NLINKS)             :: Tw         
      !real, intent(in),  dimension(:),allocatable     :: velocity      
      real, intent(in),  dimension(NLINKS)             :: velocity
      real, intent(in),  dimension(NLINKS)             :: HLINK         
      real, intent(in),  dimension(NLINKS)             :: ChSSlp   
      integer, intent(in),  dimension(NLINKS)          :: CHANXI, CHANYJ   
      real, intent(in),     dimension(NLINKS)          :: ZELEV
      integer, intent(in), dimension(NLINKS)           :: CHANLEN
      integer, intent(in), dimension(IXRT,JXRT)        :: CH_NETLNK
      real, intent(in),  dimension(IXRT,JXRT)          :: QLINK
      
      integer, intent(in)                              :: DT_STEPS       !-- number of timestep in routing
      real,    intent(in)                              :: DTCT
      real ,   intent(in), dimension(NLINKS)           :: node_area
      integer, intent(in), dimension(NLINKS)           :: TO_NODE
      integer, intent(IN), dimension(NLINKS)           :: TYPEL          !--link
      
      ! Assign the values passed in to a particular 
      this%NLINKS            = NLINKS
      this%IXRT              = IXRT
      this%JXRT              = JXRT
      this%ChanCondConst     = ChanCondConst
      this%QSTRMVOLRT        = QSTRMVOLRT
      this%QINFLOWBASE       = QINFLOWBASE
      this%QSUM              = QSUM 
      this%XS_Peremeter      = XS_Peremeter
      this%Bw                = Bw
      this%Tw                = Tw
      this%velocity          = velocity
      this%HLINK             = HLINK
      this%ChSSlp            = ChSSlp
      this%CHANXI            = CHANXI
      this%CHANYJ            = CHANYJ
      this%ZELEV             = ZELEV
      this%CHANLEN           = CHANLEN
      this%CH_NETLNK         = CH_NETLNK
      this%QLINK             = QLINK
      this%DT_STEPS          = DT_STEPS
      this%DTCT              = DTCT
      this%node_area         = node_area
      this%TO_NODE           = TO_NODE
      this%TYPEL             = TYPEL
      
   end subroutine temperature_sim_H_inputs_init

   ! --------------------------------------------------------------------------

   subroutine temperature_sim_H_inputs_destroy (this)

      implicit none
      ! the type object being destroyed

      class(temperature_sim_H_inputs), intent(inout) :: this 

   end subroutine temperature_sim_H_inputs_destroy 

   ! --------------------------------------------------------------------------
      
   subroutine vectorizing_tools_init(this, NLINKS, IXRT, JXRT, CH_NETRT, &
                                     CHANXI, CHANYJ, LAKE_MSKRT)

      implicit none
      class(vectorizing_tools), intent(inout)    :: this             ! the type object being initialized
      integer, intent(in)                        :: NLINKS, IXRT, JXRT
      integer ,intent(in), dimension(IXRT,JXRT)  :: CH_NETRT, LAKE_MSKRT
      integer, intent(in), dimension(NLINKS)     :: CHANXI, CHANYJ 
      
      this%NLINKS = NLINKS
      this%IXRT = IXRT
      this%JXRT = JXRT
      this%CH_NETRT = CH_NETRT
      this%CHANXI = CHANXI
      this%CHANYJ = CHANYJ
      this%LAKE_MSKRT = LAKE_MSKRT
      
   end subroutine vectorizing_tools_init

   ! --------------------------------------------------------------------------

   subroutine vectorizing_tools_destroy(this)

      implicit none
      ! the type object being destroyed
      class(vectorizing_tools), intent(inout) :: this 

   end subroutine vectorizing_tools_destroy 

   ! --------------------------------------------------------------------------

   subroutine vectors_and_variables_Hyd(This)

      !subroutine vectors_and_variables(ChanCondConst, QSTRMVOLRT, QINFLOWBASE, &
      !            QSUM, XS_Peremeter, Bw, Tw, velocity, HLINK, ChSSlp, &
      !            CHANXI, CHANYJ, ZELEV, CHANLEN, CH_NETLNK, QLINK, DT_STEPS, &
      !            DTCT, node_area, TO_NODE, TYPEL, CH_NETRT, LAKE_MSKRT)
      !- ChanCondConst: m/s; the i-Tree input is mm/second
      !!real, intent(in), dimension(NLINK)                :: ChanCondConst = (/0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18/)
      !real, intent(out), dimension(NLINKS)               :: ChanCondConst
      !!- QSTRMVOLRT:    m, only for one time step
      !real, intent(out), dimension(IXRT,JXRT)            :: QSTRMVOLRT
      !!- QINFLOWBASE:   m, only for one time step
      !real, intent(out), dimension(IXRT,JXRT)            :: QINFLOWBASE
      !real, intent(out), dimension(NLINKS)               :: QSUM          !- m
      !real, intent(out), dimension(NLINKS)               :: XS_Peremeter  !- m
      !real, intent(out), dimension(NLINKS)               :: Bw            !- m
      !real, intent(out), dimension(NLINKS)               :: Tw            !- m
      !real, intent(out), dimension(NLINKS)               :: velocity      !- m/s
      !real, intent(out), dimension(NLINKS)               :: HLINK         !- m
      !real, intent(out), dimension(NLINKS)               :: ChSSlp        !- m/m
      !real, intent(out), dimension(NLINKS)               :: ZELEV
      !integer, intent(out), dimension(NLINKS)            :: CHANLEN
      !!- CH_NETLNK: m, only for one time step
      !real, intent(out), dimension(IXRT,JXRT)            :: CH_NETLNK
      !real, intent(out), dimension(NLINKS,2)             :: QLINK
      !real, intent(out), dimension(NLINKS)               :: node_area
      !real, intent(out), dimension(NLINKS)               :: TO_NODE
      !integer, intent(out), dimension(IXRT,JXRT)         :: CH_NETRT                 
      !integer, intent(out), dimension(IXRT,JXRT)         :: LAKE_MSKRT              
      !integer, intent(out), dimension(NLINKS)            :: TYPEL
      !integer :: DT_STEPS
      !real    :: DTCT

      implicit none
      class(temperature_sim_H_inputs), intent(inout)  :: this  ! the type object being initialized

      !this.IXRT    = IXRT
      !this.JXRT    = JXRT
      !this.NLINKS  = NLINKS

      this.ChanCondConst = (/0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18,0.18/)
            
      this.QSTRMVOLRT  = reshape ((/0.,0.1,0.1,0.1,0.,0.,  &
                              0.,0.,0.,0.1,0.,0.,     &
                              0.,0.,0.,0.1,0.,0.,     &
                              0.1,0.1,0.1,0.1,0.1,0., &
                              0.,0.,0.,0.,0.1,0.,     &
                              0.,0.,0.,0.,0.1,0./), shape(this.QSTRMVOLRT))
 
      this.QINFLOWBASE = reshape  ((/0.,0.05,0.05,0.05,0.,0.,   &
                              0.,0.,0.,0.05,0.,0.,         &
                              0.,0.,0.,0.05,0.,0.,         &
                              0.05,0.05,0.05,0.05,0.05,0., &
                              0.,0.,0.,0.,0.05,0.,         &
                              0.,0.,0.,0.,0.05,0./), shape(this.QINFLOWBASE))
   
      this.QSUM = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0/)
                 
      this.XS_Peremeter = (/2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2,2.2/)

      this.Bw = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0/)

      this.Tw = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0/)

      this.velocity = (/5.,5.,5.,5.,5.,5.,5.,5.,5.,10.,10.,10./)

      this.HLINK = (/0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1/)

      this.ChSSlp = (/0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01/)

      this.CHANXI = (/2,3,4,4,4,1,2,3,4,5,5,5/)

      this.CHANYJ = (/1,1,1,2,3,4,4,4,4,4,5,6/)

      this.ZELEV = (/100.,97.5,95.,92.5,90.,87.5,85.,82.5,80.,77.5,75.,72.5/)

      this.CHANLEN = (/250,250,250,250,250,250,250,250,250,250,250,250/)

      this.CH_NETLNK = reshape ((/0,1,2,3,0,0,   &
                             0,0,0,4,0,0,   &
                             0,0,0,5,0,0,   &
                             6,7,8,9,10,0,  &
                             0,0,0,0,11,0,  &
                             0,0,0,0,12,0/), shape(this.CH_NETLNK))
              
      this.QLINK =reshape ((/0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, &
                   0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1/), shape(this.QLINK))
   
      this.DT_STEPS = 1

      this.DTCT = 0.001

      this.TYPEL =  (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)

      this.node_area = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)

      !- Not sure about the last variable!
      this.TO_NODE = (/2,3,4,5,9,7,8,9,10,11,12,13/)            

      ! For debugging
      print *, "DONE: module_channel_routing.f90 ; Hydraulic data "
           
   end subroutine vectors_and_variables_Hyd
       
   ! --------------------------------------------------------------------------

   subroutine vectors_and_variables_Vec(this)
               
      implicit none
      class(vectorizing_tools), intent(inout)         :: this             ! the type object being initialized
      !this.IXRT     = IXRT
      !this.JXRT     = JXRT
      !this.NLINKS   = NLINKS

      !integer, intent(out), dimension(NLINKS)            :: CHANXI
      !integer, intent(out), dimension(NLINKS)            :: CHANYJ
      !integer, intent(out), dimension(IXRT,JXRT)         :: CH_NETRT                 
      !integer, intent(out), dimension(IXRT,JXRT)         :: LAKE_MSKRT

      this.CHANXI = (/2,3,4,4,4,1,2,3,4,5,5,5/)
      this.CHANYJ = (/1,1,1,2,3,4,4,4,4,4,5,6/)

      this.CH_NETRT = reshape ((/-9999,0,0,0,-9999,-9999,        &
                                 -9999,-9999,-9999,0,-9999,-9999, &
                                 -9999,-9999,-9999,0,-9999,-9999, &
                                 0,0,0,0,0,-9999,                 &
                                 -9999,-9999,-9999,-9999,0,-9999, &
                                 -9999,-9999,-9999,-9999,0,-9999/), shape(this.CH_NETRT))

      !this.CH_NETRT = reshape ((/1,2,3,4,5,6, &
      !                      7,8,9,10,11,12, &
      !                      13,14,15,16,17,18, &
      !                      19,20,21,22,23,24, &
      !                      25,26,27,28,29,30,&
      !                      31,32,33,34,35,36/), shape(this.CH_NETRT))
   
      this.LAKE_MSKRT = reshape ((/1,-9999,-9999,-9999,1,1,         &
                              1,1,1,-9999,1,1,                 &
                              1,1,1,-9999,1,1,                 &
                              -9999,-9999,-9999,-9999,-9999,1, &
                              1,1,1,1,-9999,1,                 &
                              1,1,1,1,-9999,1/), shape(this.LAKE_MSKRT))

      ! For debugging
      print *, "DONE: module_channel_routing.f90; Vectorzing "

   end subroutine vectors_and_variables_Vec

end module module_channel_routing

