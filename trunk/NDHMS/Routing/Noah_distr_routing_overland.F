!DJG ------------------------------------------------
!DJG   SUBROUTINE OverlandRouting
!DJG ------------------------------------------------

!subroutine OverlandRouting (DT, DTRT_TER, rt_option, ixrt, jxrt,LAKE_MSKRT, &
    !           INFXSUBRT, RETDEPRT,OVROUGHRT,SOXRT, SOYRT, SFCHEADSUBRT,DHRT, &
    !           CH_NETRT, QSTRMVOLRT,LAKE_INFLORT,QBDRYRT, &
    !           QSTRMVOLTRT,QBDRYTRT, LAKE_INFLOTRT, q_sfcflx_x,q_sfcflx_y, &
    !           dist, SO8RT, SO8RT_D, &
    !           SMCTOT2,suminfxs1,suminfxsrt,smctot1,dsmctot )


subroutine OverlandRouting( &
    ovrt_data, & ! overland data structure
    DT, &
    DTRT_TER, &
    rt_option, &
    ixrt, &       ! routing grid x size
    jxrt, &      ! routing grid y size
    q_sfcflx_x, &! accumulated x flux
    q_sfcflx_y  &! accumulated y flux
    )
#ifdef MPP_LAND
    use module_mpp_land, only:  mpp_land_max_int1,  sum_real1, my_id, io_id, numprocs
    use overland_data
#endif
    implicit none

    type (overland_struct), intent(inout) :: ovrt_data
    REAL, INTENT(IN) :: DT, DTRT_TER
    integer, INTENT(IN) :: ixrt, jxrt, rt_option
    !REAL, INTENT(out), DIMENSION(IXRT,JXRT)   :: SFCHEADSUBRT  ! moved into overland_data%control
    real, intent(inout) :: q_sfcflx_x,q_sfcflx_y

    ! REMOVED VARIABLES
    !INTEGER, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: LAKE_MSKRT

    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)   :: INFXSUBRT,  &
        !          RETDEPRT,OVROUGHRT,SOXRT, SOYRT
    !REAL, INTENT(OUT), DIMENSION(IXRT,JXRT) :: SFCHEADSUBRT,DHRT
    !INTEGER, INTENT(IN), DIMENSION(IXRT,JXRT) :: CH_NETRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: QSTRMVOLRT,LAKE_INFLORT,QBDRYRT, &
        !         QSTRMVOLTRT,QBDRYTRT, LAKE_INFLOTRT, q_sfcflx_x,q_sfcflx_y

    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT,9):: dist
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT,8)  :: SO8RT
    !INTEGER SO8RT_D(IXRT,JXRT,3)

    integer  :: i,j


    !real            :: smctot2,smctot1,dsmctot
    !real            :: suminfxsrt,suminfxs1
    ! local variable
    real            :: chan_in1,chan_in2
    real            :: lake_in1,lake_in2
    real            :: qbdry1,qbdry2
    integer :: sfcrt_flag



    !DJG Third, Call Overland Flow Routing Routine...
#ifdef HYDRO_D
    print *, "Beginning OV_routing..."
    print *, "Routing method is ",rt_option, " direction."
#endif

    !DJG debug...OV Routing...
    ovrt_data%mass_balance%pre_infiltration_excess = 0.
    chan_in1=0.
    lake_in1=0.
    qbdry1=0.
    do i=1,IXRT
        do j=1,JXRT
            ovrt_data%mass_balance%pre_infiltration_excess = &
                ovrt_data%mass_balance%pre_infiltration_excess + ( ovrt_data%control%infiltration_excess(I,J)/float(IXRT*JXRT) )
            chan_in1 = chan_in1 + ( ovrt_data%streams_and_lakes%surface_water_to_channel(I,J)/float(IXRT*JXRT) )
            lake_in1 = lake_in1 + ( ovrt_data%streams_and_lakes%surface_water_to_lake(I,J)/float(IXRT*JXRT) )
            qbdry1 = qbdry1 + ( ovrt_data%control%boundary_flux(I,J)/float(IXRT*JXRT) )
        end do
    end do

#ifdef MPP_LAND
    ! not tested
    CALL sum_real1(ovrt_data%mass_balance%pre_infiltration_excess)
    CALL sum_real1(chan_in1)
    CALL sum_real1(lake_in1)
    CALL sum_real1(qbdry1)
    ovrt_data%mass_balance%pre_infiltration_excess = ovrt_data%mass_balance%pre_infiltration_excess/float(numprocs)
    chan_in1 = chan_in1/float(numprocs)
    lake_in1 = lake_in1/float(numprocs)
    qbdry1 = qbdry1/float(numprocs)
#endif


    !DJG.7.20.2007 - Global check for infxs>retdep & skip if not...(set sfcrt_flag)
    !DJG.7.20.2007 - this check will skip ov rtng when no flow is present...

    sfcrt_flag = 0

    do j=1,jxrt
        do i=1,ixrt
            if( ovrt_data%control%infiltration_excess(i,j) .gt. ovrt_data%properties%retention_depth(i,j) ) then
                sfcrt_flag = 1
                exit
            end if
        end do
        if(sfcrt_flag .eq. 1) exit
    end do

#ifdef MPP_LAND
    call mpp_land_max_int1(sfcrt_flag)
#endif
    !DJG.7.20.2007 - Global check for infxs>retdep & skip if not...(IF)

    if (sfcrt_flag.eq.1) then  !If/then for sfc_rt check...
#ifdef HYDRO_D
        write(6,*) "calling OV_RTNG "
#endif
        CALL OV_RTNG(           &
            ovrt_data, &
            DT,                   &
            DTRT_TER,             &
            IXRT,                 &
            JXRT,                 &
            rt_option,            &
            q_sfcflx_x,           &
            q_sfcflx_y)
    else
        ovrt_data%control%surface_water_head_routing = ovrt_data%control%infiltration_excess
#ifdef HYDRO_D
        print *, "No water to route overland..."
#endif
    end if  !Endif for sfc_rt check...

    !DJG.7.20.2007 - Global check for infxs>retdep & skip if not...(ENDIF)

#ifdef HYDRO_D
    print *, "OV routing called and returned..."
#endif

    !DJG Debug...OV Routing...
    ovrt_data%mass_balance%post_infiltration_excess = 0.
    chan_in2=0.
    lake_in2=0.
    qbdry2=0.
    do i=1,IXRT
        do j=1,JXRT
            ovrt_data%mass_balance%post_infiltration_excess = &
                ovrt_data%mass_balance%post_infiltration_excess + (ovrt_data%control%surface_water_head_routing(I,J)/float(IXRT*JXRT))
            chan_in2=chan_in2 + (ovrt_data%streams_and_lakes%surface_water_to_channel(I,J)/float(IXRT*JXRT))
            lake_in2=lake_in2 + (ovrt_data%streams_and_lakes%surface_water_to_lake(I,J)/float(IXRT*JXRT))
            qbdry2=qbdry2 + (ovrt_data%control%boundary_flux(I,J)/float(IXRT*JXRT))
        end do
    end do
#ifdef MPP_LAND
    ! not tested
    CALL sum_real1(ovrt_data%mass_balance%post_infiltration_excess)
    CALL sum_real1(chan_in2)
    CALL sum_real1(lake_in2)
    CALL sum_real1(qbdry2)
    ovrt_data%mass_balance%post_infiltration_excess = ovrt_data%mass_balance%post_infiltration_excess/float(numprocs)
    chan_in2 = chan_in2/float(numprocs)
    lake_in2 = lake_in2/float(numprocs)
    qbdry2 = qbdry2/float(numprocs)
#endif

#ifdef HYDRO_D
#ifdef MPP_LAND
    if(my_id .eq. IO_id) then
#endif
        print *, "OV Routing Mass Bal: "
        print *, "WB_OV!InfxsDiff", ovrt_data%mass_balance%post_infiltration_excess - ovrt_data%mass_balance%pre_infiltration_excess
        print *, "WB_OV!Infxs1", ovrt_data%mass_balance%pre_infiltration_excess
        print *, "WB_OV!Infxs2", ovrt_data%mass_balance%post_infiltration_excess
        print *, "WB_OV!ChaninDiff", chan_in2-chan_in1
        print *, "WB_OV!Chanin1", chan_in1
        print *, "WB_OV!Chanin2", chan_in2
        print *, "WB_OV!LakeinDiff", lake_in2-lake_in1
        print *, "WB_OV!Lakein1", lake_in1
        print *, "WB_OV!Lakein2", lake_in2
        print *, "WB_OV!QbdryDiff", qbdry2-qbdry1
        print *, "WB_OV!Qbdry1", qbdry1
        print *, "WB_OV!Qbdry2", qbdry2
        print *, "WB_OV!Residual", (ovrt_data%mass_balance%post_infiltration_excess - ovrt_data%mass_balance%pre_infiltration_excess)&
            -(chan_in2-chan_in1) -(lake_in2 - lake_in1) - (qbdry2 - qbdry1)
#ifdef MPP_LAND
    endif
#endif
#endif


end subroutine OverlandRouting

!DJG ------------------------------------------------
!DJG   SUBROUTINE OV_RTNG
!DJG ------------------------------------------------

!SUBROUTINE OV_RTNG(DT,DTRT_TER,IXRT,JXRT,INFXSUBRT,      &
    !  SFCHEADSUBRT,DHRT,CH_NETRT,RETDEPRT,OVROUGHRT,      &
    !  QSTRMVOLRT,QBDRYRT,QSTRMVOLTRT,QBDRYTRT,SOXRT,     &
    !  SOYRT,dist,LAKE_MSKRT,LAKE_INFLORT,LAKE_INFLOTRT,  &
    !  SO8RT,SO8RT_D,rt_option,q_sfcflx_x,q_sfcflx_y)

subroutine ov_rtng( &
        ovrt_data, &
        DT, &
        DTRT_TER, &
        IXRT, &
        JXRT, &
        rt_option, &
        q_sfcflx_x, &
        q_sfcflx_y &
        )

    !yyww
#ifdef MPP_LAND
    use module_mpp_land, only: left_id,down_id,right_id, &
        up_id,mpp_land_com_real, my_id, &
        mpp_land_sync
#endif
    use overland_data
    IMPLICIT NONE

    !DJG --------DECLARATIONS----------------------------

    type (overland_struct), intent(inout) :: ovrt_data
    INTEGER, INTENT(IN)			:: IXRT,JXRT
    REAL, INTENT(IN)			:: DT,DTRT_TER
    integer, intent(in) :: rt_option

    !INTEGER, INTENT(IN), DIMENSION(IXRT,JXRT) :: CH_NETRT
    !INTEGER, INTENT(IN), DIMENSION(IXRT,JXRT) :: LAKE_MSKRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)	:: INFXSUBRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)	:: SOXRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)	:: SOYRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT,9):: dist
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT)	:: RETDEPRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)	:: OVROUGHRT

    !REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)	:: SFCHEADSUBRT
    !REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)	:: DHRT

    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: QSTRMVOLRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: LAKE_INFLORT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: QBDRYRT
    REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: q_sfcflx_x,q_sfcflx_y
    !REAL, INTENT(INOUT)     :: QSTRMVOLTRT,QBDRYTRT,LAKE_INFLOTRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT,8)  :: SO8RT

    !DJG Local Variables

    INTEGER :: KRT,I,J,ct

    REAL, DIMENSION(IXRT,JXRT)	:: INFXS_FRAC
    REAL	:: DT_FRAC,SUM_INFXS,sum_head
    !INTEGER SO8RT_D(IXRT,JXRT,3), rt_option




    !DJG ----------------------------------------------------------------------
    ! DJG BEGIN 1-D or 2-D OVERLAND FLOW ROUTING LOOP
    !DJG ---------------------------------------------------------------------
    !DJG  Loop over 'routing time step'
    !DJG  Compute the number of time steps based on NOAH DT and routing DTRT_TER

    DT_FRAC=INT(DT/DTRT_TER)

#ifdef HYDRO_D
    write(6,*) "OV_RTNG  DT_FRAC, DT, DTRT_TER",DT_FRAC, DT, DTRT_TER
    write(6,*) "IXRT, JXRT = ",ixrt,jxrt
#endif

    !DJG NOTE: Applying all infiltration excess water at once then routing
    !DJG       Pre-existing SFHEAD gets combined with Precip. in the
    !DJG       calculation of INFXS1 during subroutine SRT.f.
    !DJG debug


    !DJG Assign all infiltration excess to surface head...
    ovrt_data%control%surface_water_head_routing=ovrt_data%control%infiltration_excess

    !DJG Divide infiltration excess over all routing time-steps
    !	     INFXS_FRAC=INFXSUBRT/(DT/DTRT_TER)

    !DJG Set flux accumulation fields to 0. before each loop...
    q_sfcflx_x = 0.
    q_sfcflx_y = 0.
    ct =0


    !DJG Execute routing time-step loop...


    DO KRT=1,DT_FRAC

        DO J=1,JXRT
            DO I=1,IXRT

                !DJG Removed 4_29_05, sfhead now updated in route_overland subroutine...
                !           SFCHEADSUBRT(I,J)=SFCHEADSUBRT(I,J)+DHRT(I,J)
                !!           SFCHEADSUBRT(I,J)=SFCHEADSUBRT(I,J)+DHRT(I,J)+INFXS_FRAC(I,J)
                !           DHRT(I,J)=0.

                !DJG ERROR Check...

                IF (ovrt_data%control%surface_water_head_routing(I,J).lt.0.) THEN
#ifdef HYDRO_D
                    print *, "ywcheck 2 ERROR!!!: Neg. Surface Head Value at (i,j):",    &
                        i,j,ovrt_data%control%surface_water_head_routing(I,J)
                    print *, "RETDEPRT(I,J) = ",ovrt_data%properties%retention_depth(I,J), "KRT=",KRT
                    print *, "INFXSUBRT(i,j)=",ovrt_data%control%infiltration_excess(i,j)
                    print *, "jxrt=",jxrt," ixrt=",ixrt
#endif
                END IF

                !DJG Remove surface water from channel cells
                !DJG Channel inflo cells specified as nonzeros from CH_NET
                !DJG 9/16/04  Channel Extractions Removed until stream model implemented...



                !yw            IF (CH_NETRT(I,J).ne.-9999) THEN
                IF (ovrt_data%streams_and_lakes%CH_NETRT(I,J).ge.0) THEN
                    ct = ct +1

                    !DJG Temporary test to up the retention depth of channel grid cells to 'soak'
                    !more water into valleys....set retdep = retdep*100 (=5 mm)

                    !	     RETDEPRT(I,J) = RETDEPRT(I,J) * 100.0    !DJG TEMP HARDWIRE!!!!
                    !	     RETDEPRT(I,J) = 10.0    !DJG TEMP HARDWIRE!!!!

                    ! AD hardwire to force channel retention depth to be 5mm.
                    ovrt_data%properties%retention_depth(I,J) = 5.0

                    IF (ovrt_data%control%surface_water_head_routing(I,J) .GT. ovrt_data%properties%retention_depth(I,J)) THEN
                        !!               QINFLO(CH_NET(I,J)=QINFLO(CH_NET(I,J)+SFCHEAD(I,J) - RETDEPRT(I,J)
                        ovrt_data%streams_and_lakes%accumulated_surface_water_to_channel = ovrt_data%streams_and_lakes%accumulated_surface_water_to_channel + &
                            (ovrt_data%control%surface_water_head_routing(I,J) - ovrt_data%properties%retention_depth(I,J))
                        ovrt_data%streams_and_lakes%surface_water_to_channel(I,J) = ovrt_data%streams_and_lakes%surface_water_to_channel(I,J) + &
                            (ovrt_data%control%surface_water_head_routing(I,J) - ovrt_data%properties%retention_depth(I,J))

                        ! if(QSTRMVOLRT(I,J) .gt. 0) then
                        !     print *, "QSTRVOL GT 0", QSTRMVOLRT(I,J),I,J
                        !  endif

                        ovrt_data%control%surface_water_head_routing(I,J) = ovrt_data%properties%retention_depth(I,J)
                    END IF
                END IF

                !DJG Lake inflow withdrawl from surface head...(4/29/05)


                IF (ovrt_data%streams_and_lakes%lake_mask(I,J) .gt. 0) THEN
                    IF (ovrt_data%control%surface_water_head_routing(I,J) .GT. ovrt_data%properties%retention_depth(I,J)) THEN
                        ovrt_data%streams_and_lakes%accumulated_surface_water_to_lake = ovrt_data%streams_and_lakes%accumulated_surface_water_to_lake + &
                            (ovrt_data%control%surface_water_head_routing(I,J) - ovrt_data%properties%retention_depth(I,J))
                        ovrt_data%streams_and_lakes%surface_water_to_lake(I,J) = ovrt_data%streams_and_lakes%surface_water_to_lake(I,J) + &
                            (ovrt_data%control%surface_water_head_routing(I,J)- ovrt_data%properties%retention_depth(I,J))
                        ovrt_data%control%surface_water_head_routing(I,J) = ovrt_data%properties%retention_depth(I,J)
                    END IF
                END IF



            END DO
        END DO

        !yw check         call MPP_LAND_COM_REAL(QSTRMVOLRT,IXRT,JXRT,99)
        !DJG----------------------------------------------------------------------
        !DJG CALL OVERLAND FLOW ROUTING SUBROUTINE
        !DJG----------------------------------------------------------------------

        !DJG Debug...


        if(rt_option .eq. 1) then
            CALL ROUTE_OVERLAND1(DTRT_TER, &
                ovrt_data%properties%distance_to_neighbor, &
                ovrt_data%control%surface_water_head_routing, &
                ovrt_data%properties%retention_depth, &
                ovrt_data%properties%roughness, &
                IXRT, JXRT, &
                ovrt_data%control%boundary_flux, &
                ovrt_data%control%boundary_flux_total,    &
                ovrt_data%properties%surface_slope, &
                ovrt_data%properties%max_surface_slope_index, &
                q_sfcflx_x, q_sfcflx_y)
        else
            CALL ROUTE_OVERLAND2(DTRT_TER, &
                ovrt_data%properties%distance_to_neighbor, &
                ovrt_data%control%surface_water_head_routing, &
                ovrt_data%control%DHRT, &
                ovrt_data%properties%surface_slope_x, &
                ovrt_data%properties%surface_slope_y, &
                ovrt_data%properties%retention_depth, &
                ovrt_data%properties%roughness, &
                IXRT, JXRT, &
                ovrt_data%control%boundary_flux, &
                ovrt_data%control%boundary_flux_total,  &
                q_sfcflx_x,q_sfcflx_y)
        end if

    END DO          ! END routing time steps

#ifdef HYDRO_D
    print *, "End of OV_routing call..."
#endif

    !----------------------------------------------------------------------
    ! END OVERLAND FLOW ROUTING LOOP
    !     CHANNEL ROUTING TO FOLLOW
    !----------------------------------------------------------------------

    !DJG ----------------------------------------------------------------
END SUBROUTINE OV_RTNG

!===================================================================================================
! Subroutine Name:
!   subroutine ROUTE_OVERLAND1
! Author(s)/Contact(s):
!   D. Gochis <gochis><ucar><edu>
! Abstract:
!   Route excess rainfall over the watershed using a 1d (steepest slope) diffusion routing scheme.
! History Log:
!   ??/??/?? -Created, DG.
!   01/14/20 -Modified, AD.
! Usage:
! Parameters:
! Input Files:
! Output Files:
! Condition codes:
! User controllable options: None.
! Notes:
!  - Adapted from CASC2D source code
!  - Returns: qsfc=DQOV   which in turn becomes DH in head calc.
!  - dh from original code has been replaced by qsfc
!  - dhh replaced by qqsfc
!===================================================================================================
SUBROUTINE ROUTE_OVERLAND1(dt, &
              dist, &
              h, &
              retent_dep, &
              roughness, &
              XX, YY, &
              QBDRY, &
              QBDRYT, &
              SO8RT, &
              SO8RT_D, &
              q_sfcflx_x, q_sfcflx_y)

#ifdef MPP_LAND
   use module_mpp_land, only: left_id, down_id, right_id, &
                               up_id, mpp_land_com_real, &
                               my_id, mpp_land_com_real8, &
                               mpp_land_sync
#endif

   IMPLICIT NONE

   ! Passed variables
   integer, intent(in) :: XX,YY ! x and y dims on the routing grid
   real, intent(in)    :: dt ! overland routing timestep (sec)

   ! Passed arrays
   real, intent(in), dimension(XX,YY,9)    :: dist ! distance to neighbor cells (m)
   real, intent(in), dimension(XX,YY,8)    :: SO8RT ! terrain slope in all directions (m/m)
   integer, intent(in), dimension(XX,YY,3) :: SO8RT_D ! steepest terrain slope cell (i, j, index)
   real, intent(in), dimension(XX,YY)      :: retent_dep ! retention depth to hold before water is routed (mm)
   real, intent(in), dimension(XX,YY)      :: roughness ! manning's roughness coefficient (sec/m^[1/3]?)
   real, intent(inout), dimension(XX,YY)   :: h ! depth of overland water (mm)
   real, intent(inout), dimension(XX,YY)   :: QBDRY ! boundary cell flux (mm)
   real, intent(inout)                     :: QBDRYT ! total boundary flux (m)
   real, intent(inout), dimension(XX,YY)   :: q_sfcflx_x, q_sfcflx_y ! directional fluxes (m)

   ! Local variables
   integer :: i, j ! indices for local loops
   integer :: i_max_neigh, j_max_neigh, indx_max_neigh ! i, j, and direction index for steepest slope neighbor
   real    :: hmax ! maximum depth for mass flux limit check (m)
   real    :: alfax ! alpha parameter (m^[1/3]/sec)
   real    :: sfx ! friction slope (m/m)
   real    :: hh ! depth of water available for overland routing (m)
   real    :: hh53 ! h to the beta (5/3) exponent (m^[5/3])
   real    :: qqsfc ! depth of water leaving cell (m) 
   real    :: h_out ! tracking variable for water depth movement (mm) 
   real    :: one_over_dist ! 1/neighbor_distance to ease slope calculations (1/m)

   ! Local arrays
   real*8, dimension(XX,YY)  :: DH ! head change tracking for self (mm)
   real*8, dimension(XX,YY)  :: DH_from_neigh ! head change tracking for neighbor (mm)
   real*8, dimension(XX,YY)  :: QBDRY_from_neigh ! temporary boundary flux accumulator (mm)

   ! Initialize variables
   DH = 0.
   DH_from_neigh = 0.
   QBDRY_from_neigh = 0.

   ! Begin Routing of Excess Rainfall over the Watershed

   ! Loop to route water
   do j=2,YY-1
      do i=2,XX-1

         if (h(I,J) .GT. retent_dep(I,J)) then

            i_max_neigh = SO8RT_D(i,j,1)
            j_max_neigh = SO8RT_D(i,j,2)
            indx_max_neigh = SO8RT_D(i,j,3)
            one_over_dist = 1.0/dist(i,j,indx_max_neigh)
            sfx = so8RT(i,j,indx_max_neigh) - (h(i_max_neigh,j_max_neigh)-h(i,j)) * 0.001 * one_over_dist

            hmax = h(i,j)*0.001  !Specify max head for mass flux limit

            if (sfx .lt. 1E-20) then
               ! ADNOTE: Updates i_max_neigh, j_max_neigh, indx_max_neigh, and sfx for steepest slope neighbor.
               !call GETMAX8DIR(i_max_neigh,j_max_neigh,I,J,H,RETENT_DEP,so8rt,dist(i,j,:),sfx,XX,YY)
               call GetSfcMax8dir(i, j, XX, YY, H, so8rt, dist(i,j,:), &
                                  i_max_neigh, j_max_neigh, indx_max_neigh, sfx)
               ! ADNOTE: Shouldn't this be updated after the neighbor check?
               ! one_over_dist = 1.0/dist(i,j,index)
            end if

            !ADNOTE: GetSfcMax8dir will return i_max_neigh=-1 if steepest slope cannot be found
            if (i_max_neigh > 0) then  ! do the rest if the lowest grid can be found.

               if (sfx .lt. 1E-20) then
#ifdef HYDRO_D
                  print*, "Message: sfx reset to 1E-20. sfx =",sfx
                  print*, "i,j,indx_max_neigh,i_max_neigh,j_max_neigh",i,j,indx_max_neigh,i_max_neigh,j_max_neigh
                  print*, "so8RT(i,j,indx_max_neigh), h(i_max_neigh,j_max_neigh), h(i,j), dist(i,j,indx_max_neigh) ", &
                           so8RT(i,j,indx_max_neigh), h(i_max_neigh,j_max_neigh), h(i,j), dist(i,j,indx_max_neigh)
#endif
                  sfx = 1E-20
               end if

               ! See CASC2D manual (Julien & Saghafian 1991)
               ! http://hdl.handle.net/10217/173165
               ! Depth-discharge relationship:
               !   q_x = alpha_x * ( h ^ beta)
               ! where alpha & beta parameters for turbulent flow over a rough boundary are:
               !   alpha_x = sqrt(slopef_x) / n
               !   beta = 5/3
               ! Variables:
               !   q_x = unit flow in x direction
               !   h = surface flow depth
               !   slopef_x = friction slope in x direction
               !   n = manning's n
               alfax = sqrt(sfx) / roughness(i,j)
               hh = (h(i,j)-retent_dep(i,j)) * 0.001
               hh53 = hh**(5./3.)

               ! Calculate q-flux
               ! Units check: (m^(1/3))/sec * (m^(5/3)) * sec * 1/m = m
               ! ADNOTE: Is this the right "dist" to use? Should this be flowpath width vs. length?
               ! OK since identical for square cells, but need to watch for other shapes.
               qqsfc = alfax * hh53 * dt * one_over_dist ! meters

               ! Courant check (simple mass limit on overland flow)
               ! ADNOTE: This check needs updating... not sure where this eqn came from, units don't match.
               ! Another option: qqsfc = min(qqsfc, 0.6*hmax)
               if (qqsfc .ge. (hmax*dt*one_over_dist)) qqsfc = hmax*dt*one_over_dist

               ! Accumulate directional fluxes on routing subgrid
               ! ADNOTE: Why do these calcs for 1-d flow?
               if (i_max_neigh .gt. i) then
                  q_sfcflx_x(I,J) = q_sfcflx_x(I,J) + qqsfc * &
                                      (1.0 - 0.5 * (ABS(j-j_max_neigh)))
               else if (i_max_neigh .lt. i) then
                  q_sfcflx_x(I,J) = q_sfcflx_x(I,J) - 1.0 * &
                                      qqsfc * (1.0 - 0.5 * (ABS(j-j_max_neigh)))
               else
                  q_sfcflx_x(I,J) = q_sfcflx_x(I,J) + 0.
               end if
               if (j_max_neigh .gt. j) then
                  q_sfcflx_y(I,J) = q_sfcflx_y(I,J) + qqsfc * &
                                      (1.0 - 0.5 * (ABS(i-i_max_neigh)))
               else if (j_max_neigh .lt. j) then
                  q_sfcflx_y(I,J) = q_sfcflx_y(I,J) - 1.0 * &
                                      qqsfc * (1.0 - 0.5 * (ABS(i-i_max_neigh)))
               else
                  q_sfcflx_y(I,J) = q_sfcflx_y(I,J) + 0.
               end if

               !DJG put adjustment in for (h) due to qqsfc
               !yw changed as following:
               h_out = qqsfc * 1000. ! mm
               if ((h(i,j) - h_out) <0 )  then
#ifdef HYDRO_D
                  print*, "Error Warning: surface head is negative:  ",i,j,i_max_neigh,j_max_neigh, &
                            h(i,j) - h_out
#endif
                  h_out = h(i,j)
               end if
               ! Remove water from self and add to neighbor
               DH(i,j) = DH(i,j) - h_out ! mm
               DH_from_neigh(i_max_neigh,j_max_neigh) = DH_from_neigh(i_max_neigh,j_max_neigh) + h_out ! mm
               !yw end change

               !DG Boundary adjustments here
               !DG Constant Flux Condition
#ifdef MPP_LAND
               if ( ((i_max_neigh.eq.XX) .and. (right_id .lt. 0)) .or. &
                    ((i_max_neigh.eq.1)  .and. (left_id  .lt. 0)) .or. &
                    ((j_max_neigh.eq.1)  .and. (down_id  .lt. 0)) .or. &
                    ((j_max_neigh.eq.YY) .and. (up_id    .lt. 0)) ) then
#else
               if ( (i_max_neigh.eq.XX) .or. (i_max_neigh.eq.1) .or. &
                    (j_max_neigh.eq.1)  .or. (j_max_neigh.eq.YY ) ) then
#endif
                  ! For boundary cells, move water from DH to QBDRY.
                  ! ADNOTE: Note that QBDRY is negative so by subtracting we are really adding to the negative flux.
                  ! ADNOTE: Not sure why we are using h_out for one and qqsfc for the other...?
                  DH_from_neigh(i_max_neigh,j_max_neigh) = DH_from_neigh(i_max_neigh,j_max_neigh) - h_out ! mm
                  QBDRY_from_neigh(i_max_neigh,j_max_neigh) = QBDRY_from_neigh(i_max_neigh,j_max_neigh) - qqsfc*1000. ! mm
                  QBDRYT = QBDRYT - qqsfc ! m
               end if

            end if ! End if lowest grid cell can be found
         end if ! End if overland flow depth > retention depth

      end do
   end do

#ifdef MPP_LAND
   ! use double precision to solve the underflow problem.
   call MPP_LAND_COM_REAL8(DH_from_neigh, XX, YY, 1)
   call MPP_LAND_COM_REAL8(QBDRY_from_neigh, XX, YY, 1)
#endif
   QBDRY = QBDRY + QBDRY_from_neigh ! mm
   DH = DH + DH_from_neigh ! mm

#ifdef MPP_LAND
   call MPP_LAND_COM_REAL8(DH, XX, YY, 99)
   call MPP_LAND_COM_REAL(QBDRY, XX, YY, 99)
#endif

   ! Make final full matrix net adjustment to overland depth
   H = H + DH

   return

end subroutine ROUTE_OVERLAND1


!DJG ----------------------------------------------------------------------

!DJG     SUBROUTINE ROUTE_OVERLAND2
!DJG ----------------------------------------------------------------

SUBROUTINE ROUTE_OVERLAND2 (dt,                               &
        &          dist,h,qsfc,sox,soy,                                   &
        &          retent_dep,dist_rough,XX,YY,QBDRY,QBDRYT,               &
        &          q_sfcflx_x,q_sfcflx_y)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  Subroutine to route excess rainfall over the watershed
    !     using a 2d diffusion routing scheme.
    !
    !  Called from: main.f
    !
    !      Will try to formulate this to be called from NOAH
    !
    !  Returns: qsfc=DQOV   which in turn becomes DH in head calc.
    !
    !  Created:  Adaptded from CASC2D source code
    !  NOTE: dh from original code has been replaced by qsfc
    !        dhh replaced by qqsfc
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPP_LAND
    use module_mpp_land, only: left_id,down_id,right_id,&
        up_id,mpp_land_com_real,MPP_LAND_UB_COM, &
        MPP_LAND_LR_COM,mpp_land_com_integer
#endif

    IMPLICIT NONE


    !! Declare Passed variables

    real :: gsize
    INTEGER, INTENT(IN) :: XX,YY
    REAL, INTENT(IN) :: dt , dist(XX,YY,9)

    !! Declare passed arrays

    REAL, INTENT(INOUT), DIMENSION(XX,YY) :: h
    REAL, INTENT(INOUT), DIMENSION(XX,YY) :: qsfc
    REAL, INTENT(IN), DIMENSION(XX,YY) :: sox
    REAL, INTENT(IN), DIMENSION(XX,YY) :: soy
    REAL, INTENT(IN), DIMENSION(XX,YY) :: retent_dep
    REAL, INTENT(IN), DIMENSION(XX,YY) :: dist_rough
    REAL, INTENT(INOUT), DIMENSION(XX,YY) :: QBDRY
    REAL, INTENT(INOUT), DIMENSION(XX,YY) :: q_sfcflx_x,q_sfcflx_y
    REAL, INTENT(INOUT) :: QBDRYT
    REAL  :: DH(XX,YY)

    !!! Declare Local Variables

    REAL :: dhdx,dhdy,alfax,alfay
    REAL :: hh53,qqsfc,hh,dt_new
    REAL :: sfx,sfy
    REAL :: tmp_adjust

    INTEGER :: i,j

    !!! Initialize variables




    !!! Begin Routing of Excess Rainfall over the Watershed


    DH = 0
    !!! Loop to route water in x-direction
    do j=1,YY
        do i=1,XX
            ! check for boundary gridpoint?
            if (i.eq.XX) GOTO 998
            gsize = dist(i,j,3)

            ! check for detention storage?
            if (h(i,j).lt.retent_dep(i,j).AND.     &
                h(i+1,j).lt.retent_dep(i+1,j)) GOTO 998

            dhdx = (h(i+1,j)/1000. - h(i,j)/1000.) / gsize  ! gisze-(m),h-(mm)

            sfx = (sox(i,j)-dhdx+1E-30)
            if (abs(sfx).lt.1E-20) sfx=1E-20
            alfax = ((abs(sfx))**0.5)/dist_rough(i,j)
            if (sfx.lt.0.) then
                hh=(h(i+1,j)-retent_dep(i+1,j))/1000.
            else
                hh=(h(i,j)-retent_dep(i,j))/1000.
            end if

            if ((retent_dep(i,j).gt.0.).AND.(hh.le.0.)) GOTO 998
            if (hh.lt.0.) then
                GOTO 998
            end if

            hh53=hh**(5./3.)

            ! Calculate q-flux... (units (m))
            qqsfc = (sfx/abs(sfx))*alfax*hh53*dt/gsize
            q_sfcflx_x(I,J) = q_sfcflx_x(I,J) + qqsfc

            !DJG put adjustment in for (h) due to qqsfc

            !yw changed as following:
            tmp_adjust=qqsfc*1000
            if(tmp_adjust .le. 0 ) GOTO 998
            if((h(i,j) - tmp_adjust) <0 )  then
#ifdef HYDRO_D
                print*, "WARNING: surface head is negative:  ",i,j
#endif
                tmp_adjust = h(i,j)
            end if
            if((h(i+1,j) + tmp_adjust) <0) then
#ifdef HYDRO_D
                print*, "WARNING: surface head is negative: ",i+1,j
#endif
                tmp_adjust = -1*h(i+1,j)
            end if
            Dh(i,j) = Dh(i,j)-tmp_adjust
            Dh(i+1,j) = Dh(i+1,j) + tmp_adjust
            !yw end change



            !DG Boundary adjustments here
            !DG Constant Flux Condition
#ifdef MPP_LAND
            if ((i.eq.1).AND.(sfx.lt.0).and. &
                    (left_id .lt. 0) ) then
#else
                if ((i.eq.1).AND.(sfx.lt.0)) then
#endif
                    Dh(i,j) = Dh(i,j) + qqsfc*1000.
                    QBDRY(I,J)=QBDRY(I,J) + qqsfc*1000.
                    QBDRYT=QBDRYT + qqsfc*1000.
#ifdef MPP_LAND
                else if ( (i.eq.(XX-1)).AND.(sfx.gt.0) &
                        .and. (right_id .lt. 0) ) then
#else
                else if ((i.eq.(XX-1)).AND.(sfx.gt.0)) then
#endif
                    tmp_adjust = qqsfc*1000.
                    if(h(i+1,j).lt.tmp_adjust) tmp_adjust = h(i+1,j)
                    Dh(i+1,j) = Dh(i+1,j) - tmp_adjust
                    !DJG Re-assign h(i+1) = 0.0 when <0.0 (from rounding/truncation error)
                    QBDRY(I+1,J)=QBDRY(I+1,J) - tmp_adjust
                    QBDRYT=QBDRYT - tmp_adjust
                end if


                998     continue

                !! End loop to route sfc water in x-direction
        end do
    end do

    H = H + DH
#ifdef MPP_LAND
    call MPP_LAND_LR_COM(H,XX,YY,99)
    call MPP_LAND_LR_COM(QBDRY,XX,YY,99)
#endif


    DH = 0
    !!!! Loop to route water in y-direction
    do j=1,YY
        do i=1,XX

            !! check for boundary grid point?
            if (j.eq.YY) GOTO 999
            gsize = dist(i,j,1)


            !! check for detention storage?
            if (h(i,j).lt.retent_dep(i,j).AND.     &
                h(i,j+1).lt.retent_dep(i,j+1)) GOTO 999

            dhdy = (h(i,j+1)/1000. - h(i,j)/1000.) / gsize

            sfy = (soy(i,j)-dhdy+1E-30)
            if (abs(sfy).lt.1E-20) sfy=1E-20
            alfay = ((abs(sfy))**0.5)/dist_rough(i,j)
            if (sfy.lt.0.) then
                hh=(h(i,j+1)-retent_dep(i,j+1))/1000.
            else
                hh=(h(i,j)-retent_dep(i,j))/1000.
            end if

            if ((retent_dep(i,j).gt.0.).AND.(hh.le.0.)) GOTO 999
            if (hh.lt.0.) then
                GOTO 999
            end if

            hh53=hh**(5./3.)

            ! Calculate q-flux...
            qqsfc = (sfy/abs(sfy))*alfay*hh53*dt / gsize
            q_sfcflx_y(I,J) = q_sfcflx_y(I,J) + qqsfc


            !DJG put adjustment in for (h) due to qqsfc
            !yw	  h(i,j) = h(i,j)-qqsfc*1000.
            !yw          h(i,j+1) = h(i,j+1) + qqsfc*1000.
            !yw changed as following:
            tmp_adjust=qqsfc*1000
            if(tmp_adjust .le. 0 ) GOTO 999

            if((h(i,j) - tmp_adjust) <0 )  then
#ifdef HYDRO_D
                print *, "WARNING: surface head is negative:  ",i,j
#endif
                tmp_adjust = h(i,j)
            end if
            if((h(i,j+1) + tmp_adjust) <0) then
#ifdef HYDRO_D
                print *, "WARNING: surface head is negative: ",i,j+1
#endif
                tmp_adjust = -1*h(i,j+1)
            end if
            Dh(i,j) = Dh(i,j)-tmp_adjust
            Dh(i,j+1) = Dh(i,j+1) + tmp_adjust
            !yw end change

            !qsfc(i,j) = qsfc(i,j)-qqsfc
            !qsfc(i,j+1) = qsfc(i,j+1) + qqsfc
            !!DG Boundary adjustments here
            !!DG Constant Flux Condition
#ifdef MPP_LAND
            if ((j.eq.1).AND.(sfy.lt.0)   &
                    .and. (down_id .lt. 0) ) then
#else
                if ((j.eq.1).AND.(sfy.lt.0)) then
#endif
                    Dh(i,j) = Dh(i,j) + qqsfc*1000.
                    QBDRY(I,J)=QBDRY(I,J) + qqsfc*1000.
                    QBDRYT=QBDRYT + qqsfc*1000.
#ifdef MPP_LAND
                else if ((j.eq.(YY-1)).AND.(sfy.gt.0) &
                        .and. (up_id .lt. 0) ) then
#else
                else if ((j.eq.(YY-1)).AND.(sfy.gt.0)) then
#endif
                    tmp_adjust = qqsfc*1000.
                    if(h(i,j+1).lt.tmp_adjust) tmp_adjust = h(i,j+1)
                    Dh(i,j+1) = Dh(i,j+1) - tmp_adjust
                    !DJG Re-assign h(j+1) = 0.0 when <0.0 (from rounding/truncation error)
                    QBDRY(I,J+1)=QBDRY(I,J+1) - tmp_adjust
                    QBDRYT=QBDRYT - tmp_adjust
                end if

                999     continue

                !!!! End loop to route sfc water in y-direction
        end do
    end do

    H = H +DH
#ifdef MPP_LAND
    call MPP_LAND_UB_COM(H,XX,YY,99)
    call MPP_LAND_UB_COM(QBDRY,XX,YY,99)
#endif
    return

    !DJG ----------------------------------------------------------------------
end subroutine ROUTE_OVERLAND2

