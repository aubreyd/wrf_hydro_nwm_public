!===================================================================================================
! Subroutine Name:
!   River Temperature Runs
! Author(s)/Contact(s):
!   R. Abdi <rabdi><ucar><edu>
! Abstract:
!   Calls the run subroutines for
!   1st: Weather Disaggregate and
!   2nd: River Temperature Physical Model
! History Log:
! Usage:
! Parameters:
! Input Files:
! Output Files:
! Condition codes:
! User controllable options: None.
! Notes:
!===================================================================================================

Module riverTemperatureExecutor

use module_lsm_forcing, only: weather2d
use config_base, only: noah_lsm, nlst
!use module_NoahMP_hrldas_driver
use module_RT_data, only: rt_domain
use disaggregateWeatherModule, only: fineWeathData
use module_river_temperature, only: riverTemp

implicit none


type  :: riverTemperatureRuns

	!type (fineWeathData), pointer   :: weatherDisaggregate => null()
	!type (riverTemp), pointer       :: riverTempSolver => null()

    contains
	procedure :: riverTemp_run => riverTemperature_run

end type riverTemperatureRuns

type(riverTemperatureRuns) :: riverTempRuns

contains

subroutine riverTemperature_run(this, did)

	implicit none
	INTEGER                                          ::  did  ! for the rt_domain
  	class(riverTemperatureRuns), intent(inout)       ::  this


	call fineWeathData%runWeatherDisagg (rt_domain(did)%IX,         &
                                                rt_domain(did)%JX,         &
						rt_domain(did)%IXRT,       &
						rt_domain(did)%JXRT,       &
						!weather2d,                &
						weather2d%T2,              &
						weather2d%short,           &
						weather2d%q2x,	           &
						weather2d%wind,	           &
            weather2d%sfcp,            &
						rt_domain(did)%NLINKS,	   &
						rt_domain(did)%overland%streams_and_lakes%CH_NETRT,    &
						rt_domain(did)%CHANXI,	                               &
						rt_domain(did)%CHANYJ,	                               &
						rt_domain(did)%overland%streams_and_lakes%lake_mask,   &
						!module_NoahMP_hrldas_driver%TMN,                      &
						!module_NoahMP_hrldas_driver%TSK,                      &
						!module_NoahMP_hrldas_driver%TSLB,                     &
						rt_domain(did)%TMN_temp,                               &
!            rt_domain(did)%TSK_temp,                               &
            rt_domain(did)%TG_temp,                                &
            rt_domain(did)%TSLB_temp,                              &
            noah_lsm%nsoil,                                        &
						nlst(did)%AGGFACTRT,                                   &
            rt_domain(did)%Tair_fine_1D,     &
						rt_domain(did)%SolarRad_fine_1D, &
						rt_domain(did)%RelH_fine_1D,     &  !To do: change this name to spec. humidity, since that's what it actually is  --TLE
						rt_domain(did)%windSpd_fine_1D,  &
						rt_domain(did)%TMN_1D_fine,      &
!            rt_domain(did)%TSK_1D_fine,      &
            rt_domain(did)%TG_1D_fine,       &
						rt_domain(did)%TSLB_1D_fine,     &
            rt_domain(did)%SFCP_1D           &
						)

  rt_domain(did)%RelH_fine_1D = 30.0
  rt_domain(did)%Tw = 3.0

  call riverTemp%runRiverTempModel(     nlst(did)%dxrt0,                &
					nlst(did)%DTRT_CH,               &
					nlst(did)%DTCT,	                 &
					rt_domain(did)%IXRT,             &
					rt_domain(did)%JXRT,             &
					rt_domain(did)%Tair_fine_1D,     &
					rt_domain(did)%RelH_fine_1D,     &
          rt_domain(did)%SFCP_1D,          &
					rt_domain(did)%windSpd_fine_1D,  &
					rt_domain(did)%C_Cloud,          &
					rt_domain(did)%TSLB_1D_fine,     &
					rt_domain(did)%SolarRad_fine_1D, &
					rt_domain(did)%XS_Peremeter,     &
					rt_domain(did)%HLINK,	           &
!					rt_domain(did)%ChSSlp,	         &
          rt_domain(did)%So,               &
					rt_domain(did)%ZELEV,	           &
					rt_domain(did)%velocity,         &
					rt_domain(did)%Tw,        	     &
!					rt_domain(did)%QSUM,	           &
!					rt_domain(did)%QINFLOWBASE,      &
!          rt_domain(did)%QSTRMVOLRT_TS,    &
          rt_domain(did)%InflowUpstream,   &
          rt_domain(did)%InflowGw,         &
          rt_domain(did)%InflowLand,       &
					rt_domain(did)%TMN_1D_fine,      &
					rt_domain(did)%TSLB_1D_fine,     &
!					rt_domain(did)%TSK_1D_fine,      &
          rt_domain(did)%TG_1D_fine,       &
          rt_domain(did)%InflowUpstreamTemp,&
					rt_domain(did)%S_Albedo,         &
					rt_domain(did)%F_Shade,          &
					rt_domain(did)%G_Chann_K,        &
					rt_domain(did)%C_Emb,	           &
					rt_domain(did)%P_Size,           &
					rt_domain(did)%CHANXI,           &
					rt_domain(did)%CHANYJ,	         &
					rt_domain(did)%CH_NETLNK,        &
					rt_domain(did)%overland%streams_and_lakes%CH_NETRT,        &
					rt_domain(did)%QLINK,	           &
					rt_domain(did)%node_area,        &
					!rt_domain(did)%FROM_NODE,       &
					rt_domain(did)%TO_NODE,	         &
					rt_domain(did)%TYPEL,	           &
					rt_domain(did)%T_QSUM,	         &
          rt_domain(did)%T_QSUM_prev,      &
          rt_domain(did)%NetRad,           &
          rt_domain(did)%NetRad_prev,      &
          rt_domain(did)%node_order,       &
          rt_domain(did)%N_UPSTR,          &
          rt_domain(did)%water_temp_counter, &
					rt_domain(did)%NLINKS            &
					)


end subroutine riverTemperature_run

End Module riverTemperatureExecutor
