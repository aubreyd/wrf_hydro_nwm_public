!===================================================================================================
! Subroutine Name:
!   subroutine disaggregateWeather
! Author(s)/Contact(s):
!   R. Abdi <rabdi><ucar><edu>
! Abstract:
!   Disaggregates states and parameters from coarse grid to fine grid for the weather data
! History Log:
! Usage:
! Parameters:
! Input Files:
! Output Files:
! Condition codes:
! User controllable options: None.
! Notes:
!===================================================================================================

MODULE disaggregateWeatherModule

!use module_lsm_forcing, only: weather2d
!use config_base, only: noah_lsm, nlst
!use module_NoahMP_hrldas_driver
!use module_RT_data, only: rt_domain

implicit none

type  :: fineWeatherData_type
   !integer                                  :: NLINKS
   real, dimension(:), allocatable           :: TAIR_FINE_1D, RELH_FINE_1D, WINDSPD_FINE_1D, SOLARRAD_FINE_1D
   !real, dimension(:), allocatable          :: TMN_1D_FINE, TSLB_1D_FINE, TSK_1D_FINE, SFCP_1D
   real, dimension(:), allocatable           :: TMN_1D_FINE, TSLB_1D_FINE, TG_1D_FINE, SFCP_1D
   !type (weather_2D_inputs), pointer        :: weather_2D => null()
   !type (vectorizing_tools), pointer        :: vectorizer => null()
   !type (NoahMP_vectorizing_tools), pointer :: NMPvectorizer => null()
   logical :: pointer_allocation_guard = .false.
   contains
   procedure :: init => fineWeatherData_init
   procedure :: destroy => fineWeatherData_destroy
   !procedure :: run => fineWeatherData_run
   procedure :: runWeatherDisagg => disaggregateWeather
end type fineWeatherData_type

type(fineWeatherData_type) :: fineWeathData

contains

! --------------------------------------------------------------------------

subroutine fineWeatherData_init(this, Tair_fine_1D, RelH_fine_1D, windSpd_fine_1D, &
         !SolarRad_fine_1D, TMN_1D_fine, TSK_1D_fine,TSLB_1D_fine, SFCP_1D)
          SolarRad_fine_1D, TMN_1D_fine, TG_1D_fine,TSLB_1D_fine, SFCP_1D)
   implicit none
   class(fineWeatherData_type), intent(inout)  :: this                   ! the type object being initialized
   !character(len=20)                          :: Reach_number_string
   real,    intent(inout), dimension(:)        :: Tair_fine_1D, RelH_fine_1D, windSpd_fine_1D
   !real,    intent(inout), dimension(:)       :: SolarRad_fine_1D, TMN_1D_fine, TSK_1D_fine,TSLB_1D_fine
   real,    intent(inout), dimension(:)        :: SolarRad_fine_1D, TMN_1D_fine, TG_1D_fine,TSLB_1D_fine
   real,    intent(inout), dimension(:)        :: SFCP_1D
   !this%Tair_fine_1D       = Tair_fine_1D
   !this%RelH_fine_1D       = RelH_fine_1D
   !this%windSpd_fine_1D    = windSpd_fine_1D
   !this%SolarRad_fine_1D   = SolarRad_fine_1D
   !this%TMN_1D_fine        = TMN_1D_fine
   !this%TSK_1D_fine        = TSK_1D_fine
   !this%TSLB_1D_fine       = TSLB_1D_fine
end subroutine fineWeatherData_init

! --------------------------------------------------------------------------

subroutine fineWeatherData_destroy(this)
   implicit none
   class(fineWeatherData_type), intent(inout) :: this ! object being destroyed
end subroutine fineWeatherData_destroy

! --------------------------------------------------------------------------

subroutine disaggregateWeather(this, IX, JX, IXRT, JXRT, Tair, SolarRad, RelH, windSpd, sfcp, &
                               NLINKS, CH_NETRT, CHANXI, CHANYJ, LAKE_MSKRT, TMN, TSK, TSLB, nsoil, AGGFACTRT, &
                Tair_fine_1D, SolarRad_fine_1D, RelH_fine_1D, windSpd_fine_1D, TMN_1D_fine, &
                TSK_1D_fine, TSLB_1D_fine, SFCP_1D)
   implicit none
   class(fineWeatherData_type), intent(inout)                ::    this
   ! Input Variables ------------------------------------------------------------------------
   integer, intent(IN)                        :: IX, JX, IXRT, JXRT, NLINKS, nsoil, AGGFACTRT
   real, intent(IN), dimension(IX,JX)         :: Tair, SolarRad, RelH, windSpd, sfcp
   real, intent(IN), dimension(IX,JX)         :: TMN, TSK
   real, intent(IN), dimension(IX,1:nsoil,JX) :: TSLB
   integer,intent(in), dimension(IXRT,JXRT)   :: CH_NETRT, LAKE_MSKRT
   integer, intent(in), dimension(NLINKS)     :: CHANXI, CHANYJ
   real, intent(inout), dimension(NLINKS)     :: Tair_fine_1D, SolarRad_fine_1D, RelH_fine_1D, windSpd_fine_1D, TMN_1D_fine
   real, intent(inout), dimension(NLINKS)     :: TSK_1D_fine, TSLB_1D_fine, SFCP_1D
   ! Local Variables ------------------------------------------------------------------------
   integer        :: i, j                 ! coarse grid loop indices
   integer        :: IXXRT, JYYRT         ! fine grid i,j coordinates
   integer        :: AGGFACYRT, AGGFACXRT ! fine grid aggregation factors
   !integer       :: AGGFACTRT = 2        ! disaggregation factor, I assume as 2, 500m to 2x250m; should be checked later
   real,  dimension(IXRT,JXRT)   :: Tair_fine, SolarRad_fine, RelH_fine, windSpd_fine
   real,  dimension(IXRT,JXRT)   :: TMN_fine, TSK_fine, TSLB_fine, SFCP_fine

   do J = 1,JX ! Start coarse grid j loop
      do I = 1,IX ! Start coarse grid i loop
         ! For now, I hard-coded the AGGFACTRT, the disaggregation factor as 4
         ! => AGGFACTRT-1 = 4-1 = 3
         do AGGFACYRT = AGGFACTRT-1,0,-1 ! Start disagg fine grid j loop
            do AGGFACXRT = AGGFACTRT-1,0,-1 ! Start disagg fine grid i loop
               IXXRT = I * AGGFACTRT - AGGFACXRT ! Define fine grid i
               JYYRT = J * AGGFACTRT - AGGFACYRT ! Define fine grid j
               Tair_fine(IXXRT,JYYRT) =Tair(I,J)
               SolarRad_fine(IXXRT,JYYRT) = SolarRad(I,J)
               RelH_fine(IXXRT,JYYRT) = RelH(I,J)
               windSpd_fine(IXXRT,JYYRT) = windSpd(I,J)
               TMN_fine(IXXRT,JYYRT)  = TMN(I,J)
               TSK_fine(IXXRT,JYYRT)  = TSK(I,J)
               TSLB_fine(IXXRT,JYYRT) = TSLB(I,nsoil,J)
               SFCP_fine(IXXRT, JYYRT) = sfcp(I,J)
            end do ! end AGGFACXRT
         end do ! end AGGFACYRT
      end do ! end I
   end do ! end J

   do i = 1,NLINKS
      if ((CH_NETRT(CHANXI(i), CHANYJ(i) ) .ge. 0)) then !.and. &
            !(LAKE_MSKRT(CHANXI(i),CHANYJ(i)) .lt.0)) then
         Tair_fine_1D(i) = Tair_fine(CHANXI(i),CHANYJ(i))
         SolarRad_fine_1D(i) = SolarRad_fine(CHANXI(i),CHANYJ(i))
         RelH_fine_1D(i) = RelH_fine(CHANXI(i),CHANYJ(i))
         windSpd_fine_1D(i) = windSpd_fine(CHANXI(i),CHANYJ(i))
         TMN_1D_fine(i) = TMN_fine(CHANXI(i),CHANYJ(i))
         TSK_1D_fine(i) = TSK_fine(CHANXI(i),CHANYJ(i))
         TSLB_1D_fine(i) = TSLB_fine(CHANXI(i),CHANYJ(i))
         SFCP_1D(i) = (SFCP_fine(CHANXI(i),CHANYJ(i))/1000.0) !conversion from Pa to kPa
      end if
      !if ((CHANXI(i) .eq. 67 .and. CHANYJ(i) .eq. 48) .or. (CHANXI(i) .eq. 151 .and. CHANYJ(i) .eq. 49)) then
      !   print *, "ADCHECK_ORIENT_WEATH: i CHANXI CHANYJ", i, CHANXI(i), CHANYJ(i)
      !   print *, "  CH_NETRT, LAKE_MSKRT", CH_NETRT(CHANXI(i), CHANYJ(i) ), LAKE_MSKRT(CHANXI(i),CHANYJ(i))
      !   print *, "  T_air H wind", Tair_fine_1D(i), RelH_fine_1D(i), windSpd_fine_1D(i)
      !endif
   end do

end subroutine disaggregateWeather

end module disaggregateWeatherModule

