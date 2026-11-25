subroutine subsurfaceRouting ( subrt_data, subrt_static, subrt_input, subrt_output)
#ifdef MPP_LAND
    use module_mpp_land, only:  mpp_land_com_real, mpp_land_com_integer
    use module_subsurface_data
    use module_subsurface_static_data
    use module_subsurface_input
    use module_subsurface_output
#endif
    implicit none
    type (subsurface_struct), intent(inout) :: subrt_data
    type (subsurface_static_interface), intent(inout) :: subrt_static
    type (subsurface_input_interface), intent(inout) :: subrt_input
    type (subsurface_output_interface), intent(inout) :: subrt_output
    !integer, INTENT(IN) :: ixrt, jxrt , nsoil, rt_option
    !REAL, INTENT(IN)                          :: DT
    !real,INTENT(IN), DIMENSION(NSOIL)      :: SLDPTH
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT) :: sub_resid
    !real,INTENT(INOUT), DIMENSION(subrt_static%IXRT,subrt_static%JXRT)::INFXSUBRT
    !real,INTENT(INOUT) :: QSUBBDRYTRT
    !REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)   :: QSUBBDRYRT, QSUBRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT,NSOIL) :: SMCRT, SMCWLTRT, SMCMAXRT, SMCREFRT


    !INTEGER :: SO8RT_D(IXRT,JXRT,3)
    !REAL :: SO8RT(IXRT,JXRT,8)
    !REAL, INTENT(IN)                          :: dist(ixrt,jxrt,9)
    !  -----local array ----------
    !REAL, DIMENSION(IXRT,JXRT)   :: ZWATTABLRT
    REAL, DIMENSION(subrt_static%ixrt,subrt_static%jxrt)   :: CWATAVAIL
    INTEGER, DIMENSION(subrt_static%ixrt,subrt_static%jxrt) :: SATLYRCHK

    CWATAVAIL = 0.
    CALL FINDZWAT(subrt_static%ixrt,subrt_static%jxrt,subrt_static%NSOIL, &
        subrt_data%grid_transform%smcrt, &
        subrt_data%grid_transform%smcmaxrt, &
        subrt_data%grid_transform%smcrefrt, &
        subrt_data%grid_transform%smcwltrt, subrt_data%properties%zsoil,SATLYRCHK,subrt_data%properties%zwattablrt, &
        CWATAVAIL,subrt_data%properties%sldpth, subrt_static%subrt_wtable)
#ifdef MPP_LAND
    call MPP_LAND_COM_REAL(subrt_data%properties%zwattablrt,subrt_static%ixrt,subrt_static%jxrt,99)
    call MPP_LAND_COM_REAL(CWATAVAIL,subrt_static%ixrt,subrt_static%jxrt,99)
    call MPP_LAND_COM_INTEGER(SATLYRCHK,subrt_static%ixrt,subrt_static%jxrt,99)
#endif


    !DJG Second, Call subsurface routing routine...
#ifdef HYDRO_D
    print *, "Beginning SUB_routing..."
    print *, "Routing method is ",subrt_static%subrt_option, " direction."
#endif

    !!!! Find saturated layer depth...
    ! Loop through domain to determine sat. layers and assign wat tbl depth...
    !    and water available for subsfc routing (CWATAVAIL)...
    ! This subroutine returns: ZWATTABLRT, CWATAVAIL and SATLYRCHK


    CALL SUBSFC_RTNG( subrt_data, subrt_static, subrt_input, subrt_output, CWATAVAIL, SATLYRCHK)

#ifdef HYDRO_D
    print *, "SUBROUTE routing called and returned..."
#endif

end subroutine subsurfaceRouting

!DJG ------------------------------------------------
!DJG   SUBROUTINE SUBSFC_RTNG
!DJG ------------------------------------------------

SUBROUTINE SUBSFC_RTNG(subrt_data, subrt_static, subrt_input, subrt_output, CWATAVAIL, SATLYRCHK)

    !       use module_mpp_land, only: write_restart_rt_3, write_restart_rt_2, &
        !            my_id
#ifdef MPP_LAND
    use module_mpp_land, only: MPP_LAND_COM_REAL, sum_real1, &
        my_id, io_id, numprocs
    use module_subsurface_data
    use module_subsurface_static_data
    use module_subsurface_input
    use module_subsurface_output
    use module_hydro_stop, only:HYDRO_stop
#endif
    IMPLICIT NONE

    !DJG -------- DECLARATIONS ------------------------

    type (subsurface_struct), intent(inout) :: subrt_data
    type (subsurface_static_interface), intent(inout) :: subrt_static
    type (subsurface_input_interface), intent(inout) :: subrt_input
    type (subsurface_output_interface), intent(inout) :: subrt_output

    !INTEGER, INTENT(IN) :: IXRT,JXRT,NSOIL
    !INTEGER, INTENT(IN) :: IXRT,JXRT

    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)    :: SOXRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)    :: SOYRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)    :: LATKSATRT
    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)    :: SOLDEPRT

    !REAL, INTENT(IN), DIMENSION(IXRT,JXRT)   :: ZWATTABLRT
    REAL, INTENT(IN), DIMENSION(subrt_static%ixrt,subrt_static%jxrt)   :: CWATAVAIL
    INTEGER, INTENT(IN), DIMENSION(subrt_static%ixrt,subrt_static%jxrt) :: SATLYRCHK


    !REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)   :: QSUBRT
    !REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)   :: QSUBBDRYRT

    !REAL, INTENT(IN)                          :: dist(ixrt,jxrt,9)
    !REAL, INTENT(IN)                          :: DT
    !REAL, INTENT(IN), DIMENSION(subrt_static%nsoil)        :: ZSOIL
    !REAL, INTENT(IN), DIMENSION(subrt_static%nsoil) 	  :: SLDPTH
    !REAL, INTENT(INOUT)                       :: QSUBBDRYTRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: INFXSUBRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT,subrt_static%nsoil) :: SMCMAXRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT,subrt_static%nsoil) :: SMCREFRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT,subrt_static%nsoil) :: SMCRT
    !REAL, INTENT(INOUT), DIMENSION(IXRT,JXRT,subrt_static%nsoil) :: SMCWLTRT

    REAL, DIMENSION(subrt_static%IXRT,subrt_static%JXRT)	:: ywtmp
    !DJG Local Variables

    INTEGER	:: I,J,KK
    !djg        INTEGER, DIMENSION(IXRT,JXRT) :: SATLYRCHK

    REAL 	:: GRDAREA
    REAL	:: SUBFLO
    REAL	:: WATAVAIL

    !INTEGER :: SO8RT_D(IXRT,JXRT,3)
    !REAL :: SO8RT(IXRT,JXRT,8)
    !INTEGER, INTENT(IN)  ::  rt_option
    integer  ::  index

    INTEGER :: DT_STEPS             !-- number of timestep in routing
    REAL :: SUBDT                !-- subsurface routing timestep
    INTEGER :: KRT                  !-- routing counter
    REAL, DIMENSION(subrt_static%IXRT,subrt_static%JXRT,subrt_static%nsoil) :: SMCTMP  !--temp store of SMC
    REAL, DIMENSION(subrt_static%IXRT,subrt_static%JXRT) :: ZWATTABLRTTMP ! temp store of ZWAT
    REAL, DIMENSION(subrt_static%IXRT,subrt_static%JXRT) :: INFXSUBRTTMP ! temp store of infilx
    !djg        REAL, DIMENSION(IXRT,JXRT) :: CWATAVAIL ! temp specif. of wat avial



    !DJG Debug Variables...
    REAL :: qsubchk,qsubbdrytmp
    REAL :: junk1,junk2,junk3,junk5,junk6,junk7
    INTEGER, PARAMETER :: double=8
    REAL (KIND=double) :: smctot1a,smctot2a
    INTEGER :: kx,count

#ifdef HYDRO_D
    ! ADCHANGE: Water balance variables
    real   :: smctot1,smctot2
    real   :: suminfxsrt1,suminfxsrt2
    real   :: qbdry1,qbdry2
    real   :: sumqsubrt1, sumqsubrt2
#endif

    !DJG -----------------------------------------------------------------
    !DJG  SUBSURFACE ROUTING LOOP
    !DJG    - SUBSURFACE ROUTING RUN ON NOAH TIMESTEP
    !DJG    - SUBSURFACE ROUITNG ONLY PERFORMED ON SATURATED LAYERS
    !DJG -----------------------------------------------------------------

    !Z. Cui -----------------------------------------------------------------
    !Move SUBDT initialization up to avoid crash when HYDRO_D is defined
    SUBDT = subrt_static%dt                !-- initialize the routing timestep to DT
    !Z. Cui -----------------------------------------------------------------
#ifdef HYDRO_D
    ! ADCHANGE: START Initial water balance variables
    ! ALL VARS in MM
    suminfxsrt1 = 0.
    qbdry1 = 0.
    smctot1 = 0.
    sumqsubrt1 = 0.
    do i=1,subrt_static%IXRT
        do j=1,subrt_static%JXRT
            suminfxsrt1 = suminfxsrt1 + subrt_input%infiltration_excess(I,J) / float(subrt_static%IXRT*subrt_static%JXRT)
            qbdry1 = qbdry1 + subrt_data%state%qsubbdryrt(i,j) / subrt_data%properties%distance_to_neighbor(i,j,9)*SUBDT*1000. / float(subrt_static%IXRT*subrt_static%JXRT)
            sumqsubrt1 = sumqsubrt1 + subrt_data%state%qsubrt(i,j) / subrt_data%properties%distance_to_neighbor(i,j,9)*SUBDT*1000. / float(subrt_static%IXRT*subrt_static%JXRT)
            do kk=1,subrt_static%nsoil
                smctot1 = smctot1 + subrt_data%grid_transform%smcrt(I,J,KK)*subrt_data%properties%sldpth(KK)*1000. / float(subrt_static%IXRT*subrt_static%JXRT)
            end do
        end do
    end do

#ifdef MPP_LAND
    ! not tested
    CALL sum_real1(suminfxsrt1)
    CALL sum_real1(qbdry1)
    CALL sum_real1(sumqsubrt1)
    CALL sum_real1(smctot1)
    suminfxsrt1 = suminfxsrt1/float(numprocs)
    qbdry1 = qbdry1/float(numprocs)
    sumqsubrt1 = sumqsubrt1/float(numprocs)
    smctot1 = smctot1/float(numprocs)
#endif
    ! END Initial water balance variables
#endif


    !yw GRDAREA=DXRT*DXRT
    ! GRDAREA=subrt_data%properties%distance_to_neighbor(i,j,9)


    !DJG debug subsfc...
    subflo = 0.0

    !DJG Set up mass balance checks...
    !         CWATAVAIL = 0.            !-- initialize subsurface watavail


    !!!! Find saturated layer depth...
    ! Loop through domain to determine sat. layers and assign wat tbl depth...
    !    and water available for subsfc routing (CWATAVAIL)...
    !
    !   This function call appears to have been moved to SubsurfaceRouting()
    !
    !
    !         CALL FINDZWAT(IXRT,JXRT,subrt_static%nsoil,SMCRT,SMCMAXRT,SMCREFRT, &
        !                             SMCWLTRT,ZSOIL,SATLYRCHK,subrt_data%properties%zwattablrt, &
        !                             CWATAVAIL,SLDPTH)




    !DJG debug variable...

    !DJG Courant check temp variable setup...
    ZWATTABLRTTMP = subrt_data%properties%zwattablrt !-- temporary storage of water table level




    !!!! Call subsurface routing subroutine...
#ifdef HYDRO_D
    print *, "calling subsurface routing subroutine...Opt. ",subrt_static%subrt_option
#endif

    ! Subsurface routing method is selected by subrt_option in hydro.namelist
    ! (independent of ovrt_option, which controls overland routing):
    !   1 = steepest descent (D8) single-direction routing
    !   2 = 2-D routing across the 4 cell faces (x and y directions)
    !   3 = multi-direction D8 (D8 total outflow split among downslope neighbors)
    ! The hydraulic gradient is selected by subrt_gradient:
    !   1 = terrain slope (kinematic; flow directions fixed by terrain)
    !   2 = water table (terrain slope adjusted by water-table-depth gradient)
    if (subrt_static%subrt_option .eq. 3) then
       CALL ROUTE_SUBSURFACE3(subrt_data%properties%distance_to_neighbor, &
                              subrt_data%properties%zwattablrt, &
                              subrt_data%properties%lksatrt, subrt_data%properties%nexprt, &
                              subrt_data%properties%soldeprt, &
                              subrt_data%properties%mfdexprt, &
                              subrt_static%IXRT, subrt_static%JXRT, &
                              subrt_data%properties%surface_slope, &
                              CWATAVAIL, SUBDT, subrt_static%subrt_gradient, &
                              subrt_data%state%qsubbdryrt, subrt_data%state%qsubbdrytrt, &
                              subrt_data%state%qsubrt)
    else if (subrt_static%subrt_option .eq. 2) then
       CALL ROUTE_SUBSURFACE2(subrt_data%properties%distance_to_neighbor, &
                              subrt_data%properties%zwattablrt, &
                              subrt_data%properties%lksatrt, subrt_data%properties%nexprt, &
                              subrt_data%properties%soldeprt, &
                              subrt_static%IXRT, subrt_static%JXRT, &
                              subrt_data%properties%surface_slope, &
                              CWATAVAIL, SUBDT, subrt_static%subrt_gradient, &
                              subrt_data%state%qsubbdryrt, subrt_data%state%qsubbdrytrt, &
                              subrt_data%state%qsubrt)
    else
       CALL ROUTE_SUBSURFACE1(subrt_data%properties%distance_to_neighbor, &
                              subrt_data%properties%zwattablrt, &
                              subrt_data%properties%lksatrt, subrt_data%properties%nexprt, &
                              subrt_data%properties%soldeprt, &
                              subrt_static%IXRT, subrt_static%JXRT, &
                              subrt_data%properties%surface_slope, &
                              subrt_data%properties%max_surface_slope_index, &
                              CWATAVAIL, SUBDT, subrt_static%subrt_gradient, &
                              subrt_data%state%qsubbdryrt, subrt_data%state%qsubbdrytrt, &
                              subrt_data%state%qsubrt)
    end if

#ifdef HYDRO_D
    write(6,*) "finish calling ROUTE_SUBSURFACE ", subrt_static%subrt_option
#endif


    !!!! Update soil moisture fields with subsurface flow...

    !!!! Loop through subsurface routing domain...
    DO I=1,subrt_static%IXRT
        DO J=1,subrt_static%JXRT

            !!DJG Check for courant condition violation...put limit on qsub
            !!DJG QSUB HAS units of m^3/s SUBFLO has units of m

            ! ADCHANGE: Moved this constraint to the ROUTE_SUBSURFACE routines
            !IF (CWATAVAIL(i,j).le.ABS(qsubrt(i,j))/subrt_data%properties%distance_to_neighbor(i,j,9)*SUBDT) THEN
            !  subrt_data%state%qsubrt(i,j) = -1.0*CWATAVAIL(i,j)
            !  SUBFLO = subrt_data%state%qsubrt(i,j)  !Units of qsubrt converted via CWATAVAIL
            !ELSE
            SUBFLO = subrt_data%state%qsubrt(i,j) / subrt_data%properties%distance_to_neighbor(i,j,9) * SUBDT !Convert qsubrt from m^3/s to m
            !END IF

            WATAVAIL = 0.  !Initialize to 0. for every cell...


            !!DJG Begin loop through soil profile to adjust soil water content
            !!DJG based on subsfc flow (SUBFLO)...

            IF (SUBFLO.GT.0) THEN ! Increase soil moist for +SUBFLO (Inflow)

                ! Loop through soil layers from bottom to top
                DO KK=subrt_static%nsoil,1,-1


                    ! Check for saturated layers
                    IF (subrt_data%grid_transform%smcrt(I,J,KK) .GE. subrt_data%grid_transform%smcmaxrt(I,J,KK)) THEN
                        IF (subrt_data%grid_transform%smcrt(I,J,KK) .GT. subrt_data%grid_transform%smcmaxrt(I,J,KK)) THEN
                            print *, "FATAL ERROR: Subsfc acct. SMCMAX exceeded...", &
                                subrt_data%grid_transform%smcrt(I,J,KK), subrt_data%grid_transform%smcmaxrt(I,J,KK),KK,i,j
                            call hydro_stop("In SUBSFC_RTNG() - SMCMAX exceeded")
                        ELSE
                        END IF
                    ELSE
                        WATAVAIL = (subrt_data%grid_transform%smcmaxrt(I,J,KK)-subrt_data%grid_transform%smcrt(I,J,KK))*subrt_data%properties%sldpth(KK)
                        IF (WATAVAIL.GE.SUBFLO) THEN
                            subrt_data%grid_transform%smcrt(I,J,KK) = subrt_data%grid_transform%smcrt(I,J,KK) + SUBFLO/subrt_data%properties%sldpth(KK)
                            SUBFLO = 0.
                        ELSE
                            SUBFLO = SUBFLO - WATAVAIL
                            subrt_data%grid_transform%smcrt(I,J,KK) = subrt_data%grid_transform%smcmaxrt(I,J,KK)
                        END IF
                    END IF

                    IF (SUBFLO.EQ.0.) EXIT
                    !                IF (SUBFLO.EQ.0.) goto 669

                END DO      ! END DO FOR SOIL LAYERS

                669           continue

                ! If all layers sat. add remaining subflo to infilt. excess...
                IF (KK.eq.0.AND.SUBFLO.gt.0.) then
                    subrt_input%infiltration_excess(I,J) = subrt_input%infiltration_excess(I,J) + SUBFLO*1000.    !Units = mm
                    SUBFLO=0.
                END IF

                !DJG Error trap...
                if (subflo.ne.0.) then
#ifdef HYDRO_D
                    print *, "Subflo (+) not expired...:",subflo, i, j, kk, subrt_data%grid_transform%smcrt(i,j,1), &
                        subrt_data%grid_transform%smcrt(i,j,2), subrt_data%grid_transform%smcrt(i,j,3), &
                        subrt_data%grid_transform%smcrt(i,j,4), subrt_data%grid_transform%smcrt(i,j,5), &
                        subrt_data%grid_transform%smcrt(i,j,6), subrt_data%grid_transform%smcrt(i,j,7), &
                        subrt_data%grid_transform%smcrt(i,j,8), "SMCMAX", subrt_data%grid_transform%smcmaxrt(i,j,1)
#endif
                end if


            ELSE IF (SUBFLO.LT.0) THEN    ! Decrease soil moist for -SUBFLO (Drainage)


                !DJG loop from satlyr back down and subtract out subflo as necess...
                !    now set to SMCREF, 8/24/07
                !DJG and then using unsat cond as opposed to Ksat...

                DO KK=SATLYRCHK(I,J),subrt_static%nsoil
                    WATAVAIL = (subrt_data%grid_transform%smcrt(I,J,KK) - subrt_data%grid_transform%smcrefrt(I,J,KK)) * subrt_data%properties%sldpth(KK)
                    IF (WATAVAIL.GE.ABS(SUBFLO)) THEN
                        !?yw mod                 IF (WATAVAIL.GE.(ABS(SUBFLO)+0.000001) ) THEN
                        subrt_data%grid_transform%smcrt(I,J,KK) = subrt_data%grid_transform%smcrt(I,J,KK) + SUBFLO/subrt_data%properties%sldpth(KK)
                        SUBFLO=0.
                    ELSE     ! Since subflo is small on a time-step following is unlikely...
                        subrt_data%grid_transform%smcrt(I,J,KK) = subrt_data%grid_transform%smcrefrt(I,J,KK)
                        SUBFLO=SUBFLO+WATAVAIL
                    END IF
                    IF (SUBFLO.EQ.0.) EXIT
                    !                IF (SUBFLO.EQ.0.) goto 668

                END DO  ! END DO FOR SOIL LAYERS
                668        continue


                !DJG Error trap...
                if(abs(subflo) .le. 1.E-7 )  subflo = 0.0  !truncate residual to 1E-7 prec.

                if (subflo.ne.0.) then
#ifdef HYDRO_D
                    print *, "Subflo (-) not expired:",i,j,subflo,CWATAVAIL(i,j)
                    print *, "zwatabl = ", subrt_data%properties%zwattablrt(i,j)
                    print *, "QSUBRT(I,J)=",subrt_data%state%qsubrt(i,j)
                    print *, "WATAVAIL = ",WATAVAIL, "kk=",kk
                    print *
#endif
                end if



            END IF  ! end if for +/- SUBFLO soil moisture accounting...




        END DO        ! END DO X dim
    END DO          ! END DO Y dim
    !!!! End loop through subsurface routing domain...

#ifdef MPP_LAND
    do i = 1, subrt_static%nsoil
        call MPP_LAND_COM_REAL(subrt_data%grid_transform%smcrt(:,:,i),subrt_static%IXRT,subrt_static%JXRT,99)
    end DO
#endif

#ifdef HYDRO_D
    ! ADCHANGE: START Final water balance variables
    ! ALL VARS in MM
    suminfxsrt2 = 0.
    qbdry2 = 0.
    smctot2 = 0.
    sumqsubrt2 = 0.
    do i=1,subrt_static%IXRT
        do j=1,subrt_static%JXRT
            suminfxsrt2 = suminfxsrt2 + subrt_input%infiltration_excess(I,J) / float(subrt_static%IXRT*subrt_static%JXRT)  !
            qbdry2 = qbdry2 + subrt_data%state%qsubbdryrt(i,j)/subrt_data%properties%distance_to_neighbor(i,j,9)*SUBDT*1000. / float(subrt_static%IXRT*subrt_static%JXRT)
            sumqsubrt2 = sumqsubrt2 + subrt_data%state%qsubrt(i,j)/subrt_data%properties%distance_to_neighbor(i,j,9)*SUBDT*1000. / float(subrt_static%IXRT*subrt_static%JXRT)
            do kk=1,subrt_static%nsoil
                smctot2 = smctot2 + subrt_data%grid_transform%smcrt(I,J,KK)*subrt_data%properties%sldpth(KK)*1000. / float(subrt_static%IXRT*subrt_static%JXRT)
            end do
        end do
    end do

#ifdef MPP_LAND
    ! not tested
    CALL sum_real1(suminfxsrt2)
    CALL sum_real1(qbdry2)
    CALL sum_real1(sumqsubrt2)
    CALL sum_real1(smctot2)
    suminfxsrt2 = suminfxsrt2/float(numprocs)
    qbdry2 = qbdry2/float(numprocs)
    sumqsubrt2 = sumqsubrt2/float(numprocs)
    smctot2 = smctot2/float(numprocs)
#endif

#ifdef MPP_LAND
    if (my_id .eq. IO_id) then
#endif
        print *, "SUBSFC Routing Mass Bal: "
        print *, "WB_SUB!QsubDiff", sumqsubrt2-sumqsubrt1
        print *, "WB_SUB!Qsub1", sumqsubrt1
        print *, "WB_SUB!Qsub2", sumqsubrt2
        print *, "WB_SUB!InfxsDiff", suminfxsrt2-suminfxsrt1
        print *, "WB_SUB!Infxs1", suminfxsrt1
        print *, "WB_SUB!Infxs2", suminfxsrt2
        print *, "WB_SUB!QbdryDiff", qbdry2-qbdry1
        print *, "WB_SUB!Qbdry1", qbdry1
        print *, "WB_SUB!Qbdry2", qbdry2
        print *, "WB_SUB!SMCDiff", smctot2-smctot1
        print *, "WB_SUB!SMC1", smctot1
        print *, "WB_SUB!SMC2", smctot2
        print *, "WB_SUB!Residual", sumqsubrt1 - ( (suminfxsrt2-suminfxsrt1) &
            + (smctot2-smctot1) )
#ifdef MPP_LAND
    endif
#endif
    ! END Final water balance variables
#endif


    !DJG ----------------------------------------------------------------
END SUBROUTINE SUBSFC_RTNG
!DJG ----------------------------------------------------------------


!DJG ------------------------------------------------------------------------
!DJG  SUBSURFACE FINDZWAT
!DJG ------------------------------------------------------------------------
! Finds the contiguous stack of layers (from the bottom up) with soil moisture at or
! above field capacity (SMCREF), the water available for lateral routing in that stack
! (CWATAVAIL, water above SMCREF), and the water table depth (ZWATTABLRT):
!   WTOPT = 1: depth to the top of the highest layer in the stack (layer-based; the
!              water table jumps by a full layer thickness as layers cross SMCREF)
!   WTOPT = 2: continuous; the mobile water in the stack is expressed as an equivalent
!              saturated thickness, hsat = sum((SMC - SMCREF) / (SMCMAX - SMCREF) * dz),
!              and the water table depth is soil depth - hsat. The water table then
!              changes smoothly with storage (specific yield = SMCMAX - SMCREF).
! SATLYRCHK and CWATAVAIL are the same for both options.
SUBROUTINE FINDZWAT(IXRT,JXRT,NSOIL,SMCRT,SMCMAXRT,SMCREFRT, &
        SMCWLTRT,ZSOIL,SATLYRCHK,ZWATTABLRT,CWATAVAIL,&
        SLDPTH,WTOPT)

    IMPLICIT NONE

    !DJG -------- DECLARATIONS ------------------------

    INTEGER, INTENT(IN) :: IXRT,JXRT,NSOIL
    REAL, INTENT(IN), DIMENSION(IXRT,JXRT,NSOIL) :: SMCMAXRT
    REAL, INTENT(IN), DIMENSION(IXRT,JXRT,NSOIL) :: SMCREFRT
    REAL, INTENT(IN), DIMENSION(IXRT,JXRT,NSOIL) :: SMCRT
    REAL, INTENT(IN), DIMENSION(IXRT,JXRT,NSOIL) :: SMCWLTRT
    REAL, INTENT(IN), DIMENSION(NSOIL)        :: ZSOIL
    REAL, INTENT(IN), DIMENSION(NSOIL)        :: SLDPTH
    REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)   :: ZWATTABLRT
    REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)   :: CWATAVAIL
    INTEGER, INTENT(OUT), DIMENSION(IXRT,JXRT) :: SATLYRCHK
    INTEGER, INTENT(IN) :: WTOPT ! 1 = top of saturated layer, 2 = continuous from storage

    !DJG Local Variables
    INTEGER :: KK,i,j
    REAL    :: hsat ! equivalent saturated thickness of mobile water (m)
    REAL    :: drainPor ! drainable porosity, SMCMAX - SMCREF (m3/m3)

    !!!! Find saturated layer depth...
    ! Loop through domain to determine sat. layers and assign wat tbl depth...


    SATLYRCHK = 0  !set flag for sat. layers
    CWATAVAIL = 0.  !set wat avail for subsfc rtng = 0.

    DO J=1,JXRT
        DO I=1,IXRT

            ! Loop through soil layers from bottom to top
            DO KK=NSOIL,1,-1

                ! Check for saturated layers
                ! Add additional logical check to ensure water is 'available' for routing,
                !  (i.e. not 'frozen' or otherwise immobile)
                !                IF (SMCRT(I,J,KK).GE.SMCMAXRT(I,J,KK).AND.SMCMAXRT(I,J,KK) &
                    !                  .GT.SMCWLTRT(I,J,KK)) THEN
                IF ( (SMCRT(I,J,KK).GE.SMCREFRT(I,J,KK)).AND.(SMCREFRT(I,J,KK) &
                        .GT.SMCWLTRT(I,J,KK)) ) THEN
                    ! Add additional check to ensure saturation from bottom up only...8/8/05
                    IF((SATLYRCHK(I,J).EQ.KK+1) .OR. (KK.EQ.NSOIL) ) SATLYRCHK(I,J) = KK
                END IF

            END DO


            ! Designate ZWATTABLRT based on highest sat. layer and
            ! Define amount of water avail for subsfc routing on each gridcell (CWATAVAIL)
            !  note: using a 'field capacity' value of SMCREF as lower limit...

            IF (SATLYRCHK(I,J).ne.0) then
                IF (SATLYRCHK(I,J).ne.1) then  ! soil column is partially sat.
                    ZWATTABLRT(I,J) = -ZSOIL(SATLYRCHK(I,J)-1)
                    !DJG 2/16/2016 fix                  DO KK=SATLYRCHK(I,J),NSOIL
                    !old                   CWATAVAIL(I,J) = (SMCRT(I,J,SATLYRCHK(I,J))-&
                        !old                                    SMCREFRT(I,J,SATLYRCHK(I,J))) * &
                        !old                                    (ZSOIL(SATLYRCHK(I,J)-1)-ZSOIL(NSOIL))
                    !DJG 2/16/2016 fix                    CWATAVAIL(I,J) = CWATAVAIL(I,J)+(SMCRT(I,J,KK)- &
                        !DJG 2/16/2016 fix                                     SMCREFRT(I,J,KK))*SLDPTH(KK)
                    !DJG 2/16/2016 fix                  END DO


                ELSE  ! soil column is fully saturated to sfc.
                    ZWATTABLRT(I,J) = 0.
                    !DJG 2/16/2016 fix                  DO KK=SATLYRCHK(I,J),NSOIL
                    !DJG 2/16/2016 fix                    CWATAVAIL(I,J) = (SMCRT(I,J,KK)-SMCREFRT(I,J,KK))*SLDPTH(KK)
                    !DJG 2/16/2016 fix                  END DO
                END IF
                !DJG 2/16/2016 fix accumulation of CWATAVAIL...
                DO KK=SATLYRCHK(I,J),NSOIL
                    CWATAVAIL(I,J) = CWATAVAIL(I,J)+(SMCRT(I,J,KK)- &
                        SMCREFRT(I,J,KK))*SLDPTH(KK)
                END DO

                ! Continuous water table from mobile storage in the saturated stack
                IF (WTOPT .eq. 2) THEN
                    hsat = 0.
                    DO KK=SATLYRCHK(I,J),NSOIL
                        drainPor = SMCMAXRT(I,J,KK) - SMCREFRT(I,J,KK)
                        IF (drainPor .gt. 1.e-6) THEN
                            hsat = hsat + SLDPTH(KK) * &
                                min(1., max(0., (SMCRT(I,J,KK) - SMCREFRT(I,J,KK)) / drainPor))
                        ELSE
                            ! no drainable pore space: layer at/above SMCREF counts as saturated
                            hsat = hsat + SLDPTH(KK)
                        END IF
                    END DO
                    ZWATTABLRT(I,J) = min(-ZSOIL(NSOIL), max(0., -ZSOIL(NSOIL) - hsat))
                END IF
            ELSE  ! no saturated layers...
                ZWATTABLRT(I,J) = -ZSOIL(NSOIL)
                SATLYRCHK(I,J) = NSOIL + 1
            END IF

        END DO
    END DO


    !DJG ----------------------------------------------------------------
END SUBROUTINE FINDZWAT
!DJG ----------------------------------------------------------------

!===================================================================================================
! Subroutine Name:
!   subroutine ROUTE_SUBSURFACE1
! Author(s)/Contact(s):
!   D. Gochis <gochis><ucar><edu>
! Abstract:
!   Calculates single direction (steepest slope) subsurface flow
! History Log:
!   3/27/03 -Created, DG.
!   1/05/04 -Modified, DG.
!   1/07/19 -Modified, AD.
!   2026    -Added gradOpt: 1 = kinematic (terrain slope, fixed direction),
!            2 = water-table gradient (direction re-evaluated each step, as since 2019).
! Usage:
! Parameters:
! Input Files:
! Output Files:
! Condition codes:
! User controllable options: None.
! Notes:
!  - Adapted from Wigmosta, 1994
!  - Returns qsub=DQSUB which in turn becomes SUBFLO in head calc.
!===================================================================================================
SUBROUTINE ROUTE_SUBSURFACE1(dist, z, latksat, nexp, soldep, &
                             XX, YY, SO8RT, SO8RT_D, &
                             CWATAVAIL, SUBDT, gradOpt, QSUBDRY, QSUBDRYT, qsub)
#ifdef MPP_LAND
   use module_mpp_land, only: left_id, down_id, right_id, up_id, &
                              mpp_land_com_real8, my_id, mpp_land_com_real
#endif
   use module_hydro_stop, only:HYDRO_stop

   implicit none

   ! Passed variables
   integer, intent(in) :: XX,YY
   ! Passed arrays
   real, intent(in), dimension(XX,YY,9)    :: dist ! distance to neighbor cells (m)
   real, intent(in), dimension(XX,YY)      :: z ! depth to water table (m)
   real, intent(in), dimension(XX,YY)      :: latksat
                                              ! lateral saturated hydraulic conductivity (m/s)
   real, intent(in), dimension(XX,YY)      :: nexp ! latksat decay coefficient
   real, intent(in), dimension(XX,YY)      :: soldep ! soil depth (m)
   real, intent(in), dimension(XX,YY,8)    :: SO8RT ! terrain slope in all directions (m/m)
   integer, intent(in), dimension(XX,YY,3) :: SO8RT_D ! steepest terrain slope cell (i, j, index)
   real, intent(in), dimension(XX,YY)      :: CWATAVAIL ! water available for routing (m)
   real, intent(in)                        :: SUBDT ! subsurface routing timestep (s)
   integer, intent(in)                     :: gradOpt
                                              ! hydraulic gradient: 1 = terrain slope (kinematic),
                                              ! 2 = water table (terrain + water-table-depth gradient)
   real, intent(inout), dimension(XX,YY)   :: QSUBDRY ! subsurface flow at domain boundary (m3/s)
   real, intent(inout)                     :: QSUBDRYT
                                              ! total subsurface flow outside of domain (m3/s)
   real, intent(out), dimension(XX,YY)     :: qsub ! subsurface flow outside of cell (m3/s)

   ! Local variables
   integer :: i, j ! indices for local loops
   integer :: IXX0, JYY0, index ! i, j, and direction index for steepest slope neighbor
   real    :: beta ! total head slope (m/m)
   real    :: hh ! interim variable for subsurface solution per DHSVM (dimensionless)
   real    :: gamma ! interim variable for subsurface solution per DHSVM (m3/s)
   real    :: ksat ! saturated hydraulic conductivity (m/s)
   real    :: qqsub ! net subsurface flow out of cell (m3/s?)
   real    :: waterToRoute ! total water to route from cell (m)
   real    :: width ! flow width perpendicular to the flow direction (m)
   ! Local arrays
   real*8, DIMENSION(XX,YY) :: qsub_tmp, QSUBDRY_tmp ! temp trackers for fluxes (m3/s)
   ! Local parameters
   real :: tmp_dist(9)

   ! Initialize temp variables
   qsub_tmp = 0.
   QSUBDRY_tmp = 0.

#ifdef HYDRO_D
   write(6,*) "call subsurface routing xx= , yy =", yy, xx
#endif

   do j=2,YY-1 ! start j loop

      do i=2,XX-1 ! start i loop

         if (i.ge.2.AND.i.le.XX-1.AND.j.ge.2.AND.j.le.YY-1) then !if grdcl chk
            ! check for boundary grid point?

            ! Steepest terrain slope neighbor
            IXX0 = SO8RT_D(i,j,1)
            JYY0 = SO8RT_D(i,j,2)
            index = SO8RT_D(i,j,3)

            if (gradOpt .eq. 2) then
               ! Water-table gradient: find the neighbor with the steepest terrain+head
               ! slope. Updates IXX0, JYY0, index, and beta for that neighbor.
               ! Note that beta will be -1 if steepest slope cannot be found.
               tmp_dist = dist(i,j,:)
               call GETSUB8(i, j, XX, YY, Z, so8rt, tmp_dist, &
                            IXX0, JYY0, index, beta)
            else
               ! Kinematic: hydraulic gradient = terrain slope, direction fixed by terrain
               ! (no flow from cells with no downslope neighbor, beta <= 0)
               beta = SO8RT(i,j,index)
            endif

            if (dist(i,j,index) .le. 0) then
                write(6,*) "FATAL ERROR: dist(i,j,index) is <= zero "
                call hydro_stop("In ROUTE_SUBSURFACE1() - dist(i,j,index) is <= zero ")
            endif
            if (soldep(i,j) .eq. 0) then
                call hydro_stop("In ROUTE_SUBSURFACE1() - soldep is = zero")
            endif

            if (beta .gt. 0) then            !if-then for flux calc
               if (beta .lt. 1E-20 ) then
#ifdef HYDRO_D
                  print*, "Message: beta need to be reset to 1E-20. beta = ",beta
#endif
                  beta = 1E-20
               endif

               ! Do the rest if the lowest grid can be found.
               hh = ( 1 - ( z(i,j) / soldep(i,j) ) ) ** nexp(i,j)
               ksat = latksat(i,j)

               if (hh .lt. 0.) then
                  print *, "hsub<0 at gridcell...", i,j,hh,z(i+1,j),z(i,j), &
                       soldep(i,j)
                  call hydro_stop("In ROUTE_SUBSURFACE1() - hsub<0 at gridcell ")
               endif

               ! Calculate flux from cell
               ! AD_NOTE: gamma and qqsub are negative when flow is out of cell
               ! Flow width is cell area / flow length (dx for cardinal flow, dx/sqrt(2)
               ! for diagonal flow on a square grid), i.e. the spacing between adjacent
               ! D8 flow paths. Previously the flow length dist(i,j,index) was used as
               ! the width, which doubled diagonal fluxes.
               width = dist(i,j,9) / dist(i,j,index)
               gamma = -1.0 * ( (width * ksat * soldep(i,j)) / nexp(i,j) ) * beta
               qqsub = gamma * hh

               ! AD_NOTE: Moved this water available constraint from outside qsub calc loop
               ! to inside to better account for adjustments to adjacent cells.
               ! Calculate total water to route (where dist(i,j,9) is cell area):
               waterToRoute = ABS(qqsub) / dist(i,j,9) * SUBDT
               if ( (qqsub .le. 0.0) .and. (CWATAVAIL(i,j) .lt. waterToRoute) ) THEN
                  qqsub = -1.0 * CWATAVAIL(i,j) * dist(i,j,9) / SUBDT
               endif

               ! Remove from cell qsub to track net fluxes over full i, j loop
               ! (remember: qqsub is negative when flow is out of cell!)
               qsub(i,j) = qsub(i,j) + qqsub
               ! Add to neighbor cell qsub to track net fluxes over full i, j loop
               qsub_tmp(ixx0,jyy0) = qsub_tmp(ixx0,jyy0) - qqsub

               !DJG Error Checks...
               if (qqsub .gt. 0) then
                  print*, "FATAL ERROR: qqsub should be negative, qqsub =",qqsub,&
                     "gamma=",gamma,"hh=",hh,"beta=",beta,&
                     "so8RT=",so8RT(i,j,index),"latksat=",ksat, &
                     "tan(beta)=",tan(beta),i,j,z(i,j),z(IXX0,JYY0)
                  print*, "ixx0=",ixx0, "jyy0=",jyy0
                  print*, "soldep =", soldep(i,j), "nexp=",nexp(i,j)
                  call hydro_stop("In ROUTE_SUBSURFACE1() - qqsub should be negative")
               endif

               ! Make boundary adjustments if cells are on edge of local domain
#ifdef MPP_LAND
               if ( ((ixx0.eq.XX).and.(right_id .lt. 0)) .or. &
                  ((ixx0.eq.1) .and.(left_id  .lt. 0)) .or. &
                  ((jyy0.eq.1) .and.(down_id  .lt. 0)) .or. &
                  ((JYY0.eq.YY).and.(up_id    .lt. 0)) ) then
#else
               if ((ixx0.eq.1).or.(ixx0.eq.xx).or.(jyy0.eq.1).or.(jyy0.eq.yy)) then
#endif
                  ! If on edge, move flux tracking BACK from neighbor cell qsub and into
                  ! boundary qsub (note: qsubdry is negative and qqsub is negative, so we
                  ! are making it a bigger sink here)
                  qsub_tmp(ixx0,jyy0) = qsub_tmp(ixx0,jyy0) + qqsub
                  QSUBDRY_tmp(ixx0,jyy0) = QSUBDRY_tmp(ixx0,jyy0) + qqsub
                  ! Add to total BOUNDARY qsub
                  QSUBDRYT = QSUBDRYT + qqsub
               endif

            endif  ! end if for flux calc

         endif   ! end if for gridcell check

      end do  ! end i loop

   end do   ! end j loop

#ifdef MPP_LAND
   call MPP_LAND_COM_REAL8(qsub_tmp,XX,YY,1)
   call MPP_LAND_COM_REAL8(QSUBDRY_tmp,XX,YY,1)
#endif

   ! Sum grids to get net of self and neighbor fluxes after looping through all cells
   ! (i.e., net of out-flux from one cell and in-flux from its neighbor cell)
   qsub = qsub + qsub_tmp
   QSUBDRY = QSUBDRY + QSUBDRY_tmp


#ifdef MPP_LAND
   call MPP_LAND_COM_REAL(qsub,XX,YY,99)
   call MPP_LAND_COM_REAL(QSUBDRY,XX,YY,99)
#endif

   return

end subroutine ROUTE_SUBSURFACE1

!===================================================================================================



!===================================================================================================
! Subroutine Name:
!   subroutine ROUTE_SUBSURFACE2
! Author(s)/Contact(s):
!   D. Gochis <gochis><ucar><edu> (original), rewritten 2026
! Abstract:
!   Calculates 2-D subsurface flow across the four faces (N, E, S, W) of each cell.
!   Flux across each face follows the Dupuit-Forchheimer / DHSVM (Wigmosta, 1994)
!   transmissivity form used by ROUTE_SUBSURFACE1:
!      q = T * beta * w,   T = (ksat * soldep / n) * (1 - z/soldep)**n
!   where beta is the hydraulic gradient across the face (the terrain slope for the
!   kinematic option, or the terrain slope minus the water-table-depth gradient for the
!   water-table option), w is the face width, and T is taken from the upwind cell. Each face flux is computed once, by its upwind cell, so mass is
!   exchanged pairwise and conserved.
! History Log:
!   3/27/03 -Created, DG.
!   1/05/04 -Modified, DG.
!   2026    -Rewritten: fixed head-slope sign (was terrain slope PLUS water-table-depth
!            gradient), upwind transmissivity uses upwind cell soil depth and n exponent
!            (was hard-coded n=1), face width now perpendicular cell dimension, added
!            CWATAVAIL limit on total cell outflow, boundary handling consistent with
!            ROUTE_SUBSURFACE1. Selected with subrt_option = 2 in hydro.namelist.
! Usage:
! Parameters:
! Input Files:
! Output Files:
! Condition codes:
! User controllable options: None.
! Notes:
!  - gradOpt selects the hydraulic gradient: 1 = terrain slope across each face
!    (kinematic), 2 = terrain slope minus the water-table-depth gradient.
!  - Only interior cells (2..XX-1, 2..YY-1) act as sources, matching ROUTE_SUBSURFACE1.
!    Flux into a cell on the outer domain boundary is removed as boundary flux (QSUBDRY).
!  - If a cell's total outflow over SUBDT would exceed its available water (CWATAVAIL),
!    all of its face outflows are scaled down proportionally.
!  - Returns qsub (net flux into cell, m3/s), which becomes SUBFLO in SUBSFC_RTNG.
!===================================================================================================
SUBROUTINE ROUTE_SUBSURFACE2(dist, z, latksat, nexp, soldep, &
                             XX, YY, SO8RT, &
                             CWATAVAIL, SUBDT, gradOpt, QSUBDRY, QSUBDRYT, qsub)
#ifdef MPP_LAND
   use module_mpp_land, only: left_id, down_id, right_id, up_id, &
                              mpp_land_com_real8, mpp_land_com_real
#endif
   use module_hydro_stop, only:HYDRO_stop

   implicit none

   ! Passed variables
   integer, intent(in) :: XX,YY
   ! Passed arrays
   real, intent(in), dimension(XX,YY,9)    :: dist ! distance to neighbor cells (m), 9 = cell area (m2)
   real, intent(in), dimension(XX,YY)      :: z ! depth to water table (m)
   real, intent(in), dimension(XX,YY)      :: latksat
                                              ! lateral saturated hydraulic conductivity (m/s)
   real, intent(in), dimension(XX,YY)      :: nexp ! latksat decay coefficient
   real, intent(in), dimension(XX,YY)      :: soldep ! soil depth (m)
   real, intent(in), dimension(XX,YY,8)    :: SO8RT ! terrain slope in all directions (m/m)
   real, intent(in), dimension(XX,YY)      :: CWATAVAIL ! water available for routing (m)
   real, intent(in)                        :: SUBDT ! subsurface routing timestep (s)
   integer, intent(in)                     :: gradOpt
                                              ! hydraulic gradient: 1 = terrain slope (kinematic),
                                              ! 2 = water table (terrain + water-table-depth gradient)
   real, intent(inout), dimension(XX,YY)   :: QSUBDRY ! subsurface flow at domain boundary (m3/s)
   real, intent(inout)                     :: QSUBDRYT
                                              ! total subsurface flow outside of domain (m3/s)
   real, intent(out), dimension(XX,YY)     :: qsub ! net subsurface flow into cell (m3/s)

   ! Local parameters: face neighbors N, E, S, W and their index in dist/SO8RT
   integer, parameter :: nface = 4
   integer, parameter :: faceDI(nface) = (/ 0, 1, 0, -1 /)
   integer, parameter :: faceDJ(nface) = (/ 1, 0, -1, 0 /)
   integer, parameter :: faceK(nface)  = (/ 1, 3, 5, 7 /)

   ! Local variables
   integer :: i, j, f, ni, nj ! indices
   real    :: area ! cell area (m2)
   real    :: hh ! relative saturated thickness term (1 - z/D)**n (dimensionless)
   real    :: trans ! transmissivity of this (upwind) cell (m2/s)
   real    :: beta ! total head slope across a face (m/m), > 0 means flow out of cell
   real    :: width ! face width (m)
   real    :: qout(nface) ! outflow across each face (m3/s, positive out)
   real    :: qouttot ! total outflow from cell (m3/s)
   real    :: scl ! scale factor to limit outflow to available water
   real    :: qq ! signed face flux (m3/s, negative = out of this cell)
   logical :: onBoundary
   ! Local arrays
   real*8, dimension(XX,YY) :: qsub_tmp, QSUBDRY_tmp ! temp trackers for fluxes (m3/s)

   qsub_tmp = 0.
   QSUBDRY_tmp = 0.

   do j = 2, YY-1
      do i = 2, XX-1

         if (soldep(i,j) .le. 0.) then
            call hydro_stop("In ROUTE_SUBSURFACE2() - soldep is <= zero")
         endif

         hh = 1.0 - z(i,j) / soldep(i,j)
         if (hh .lt. 0.) then
            print *, "hsub<0 at gridcell...", i, j, hh, z(i,j), soldep(i,j)
            call hydro_stop("In ROUTE_SUBSURFACE2() - hsub<0 at gridcell")
         endif
         hh = hh ** nexp(i,j)
         trans = latksat(i,j) * soldep(i,j) / nexp(i,j) * hh

         area = dist(i,j,9)
         qout = 0.
         if (trans .gt. 0. .and. CWATAVAIL(i,j) .gt. 0.) then
            do f = 1, nface
               ni = i + faceDI(f)
               nj = j + faceDJ(f)
               if (dist(i,j,faceK(f)) .le. 0.) then
                  call hydro_stop("In ROUTE_SUBSURFACE2() - dist(i,j,face) is <= zero")
               endif
               if (gradOpt .eq. 2) then
                  ! Total head slope = terrain slope - water-table-depth gradient
                  ! (head = elev - z, so a deeper water table here lowers outflow)
                  beta = SO8RT(i,j,faceK(f)) - (z(i,j) - z(ni,nj)) / dist(i,j,faceK(f))
               else
                  ! Kinematic: hydraulic gradient = terrain slope across the face
                  beta = SO8RT(i,j,faceK(f))
               endif
               if (beta .gt. 0.) then
                  ! face width = cell dimension perpendicular to the flow direction
                  width = area / dist(i,j,faceK(f))
                  qout(f) = trans * beta * width
               endif
            end do
         endif

         qouttot = sum(qout)
         if (qouttot .le. 0.) cycle

         ! Limit total outflow over the timestep to the water available in the cell
         scl = 1.0
         if (qouttot * SUBDT / area .gt. CWATAVAIL(i,j)) then
            scl = CWATAVAIL(i,j) * area / SUBDT / qouttot
         endif

         do f = 1, nface
            if (qout(f) .le. 0.) cycle
            ni = i + faceDI(f)
            nj = j + faceDJ(f)
            qq = -1.0 * qout(f) * scl   ! negative = out of cell (i,j)

            qsub_tmp(i,j)   = qsub_tmp(i,j)   + qq
            qsub_tmp(ni,nj) = qsub_tmp(ni,nj) - qq

            ! Flux into a cell on the outer domain boundary leaves the domain
#ifdef MPP_LAND
            onBoundary = ((ni .eq. XX) .and. (right_id .lt. 0)) .or. &
                         ((ni .eq. 1)  .and. (left_id  .lt. 0)) .or. &
                         ((nj .eq. 1)  .and. (down_id  .lt. 0)) .or. &
                         ((nj .eq. YY) .and. (up_id    .lt. 0))
#else
            onBoundary = (ni .eq. 1) .or. (ni .eq. XX) .or. (nj .eq. 1) .or. (nj .eq. YY)
#endif
            if (onBoundary) then
               qsub_tmp(ni,nj)    = qsub_tmp(ni,nj) + qq
               QSUBDRY_tmp(ni,nj) = QSUBDRY_tmp(ni,nj) + qq
               QSUBDRYT           = QSUBDRYT + qq
            endif
         end do

      end do
   end do

#ifdef MPP_LAND
   ! Sum flux contributions written into halo cells back onto their owning tiles
   call MPP_LAND_COM_REAL8(qsub_tmp,XX,YY,1)
   call MPP_LAND_COM_REAL8(QSUBDRY_tmp,XX,YY,1)
#endif

   qsub = qsub_tmp
   QSUBDRY = QSUBDRY + QSUBDRY_tmp

#ifdef MPP_LAND
   call MPP_LAND_COM_REAL(qsub,XX,YY,99)
   call MPP_LAND_COM_REAL(QSUBDRY,XX,YY,99)
#endif

end subroutine ROUTE_SUBSURFACE2
!===================================================================================================

!===================================================================================================
! Subroutine Name:
!   subroutine ROUTE_SUBSURFACE3
! Author(s)/Contact(s):
!   WRF-Hydro team, 2026
! Abstract:
!   Multi-direction D8 subsurface flow. The total outflow from each cell is computed
!   exactly as in ROUTE_SUBSURFACE1 (steepest slope, D8 flow width, same water-available
!   limit), and is then partitioned among all downslope neighbors (of 8) with weights
!      f_k = beta_k**p / sum(beta_j**p),   beta_k > 0
!   where beta_k is the hydraulic gradient to neighbor k (terrain slope, or terrain slope
!   minus water-table-depth gradient, per gradOpt) and p is the slope exponent MFDEXP
!   (read from HYDROTBL_F; default 4). As p -> infinity this reduces to single-direction
!   D8; smaller p spreads flow more. Because the weights change smoothly with the
!   gradients, flow does not switch abruptly between nearly equal neighbors.
! Usage:
!   Selected with subrt_option = 3 in hydro.namelist.
! Notes:
!  - Only interior cells (2..XX-1, 2..YY-1) act as sources; flux into a cell on the outer
!    domain boundary is removed as boundary flux (QSUBDRY), as in ROUTE_SUBSURFACE1.
!  - Weights are computed as (beta_k/beta_max)**p to avoid overflow/underflow at large p.
!  - The last receiving neighbor gets the remainder of the total, so mass is conserved
!    exactly.
!===================================================================================================
SUBROUTINE ROUTE_SUBSURFACE3(dist, z, latksat, nexp, soldep, mfdexp, &
                             XX, YY, SO8RT, &
                             CWATAVAIL, SUBDT, gradOpt, QSUBDRY, QSUBDRYT, qsub)
#ifdef MPP_LAND
   use module_mpp_land, only: left_id, down_id, right_id, up_id, &
                              mpp_land_com_real8, mpp_land_com_real
#endif
   use module_hydro_stop, only:HYDRO_stop

   implicit none

   ! Passed variables
   integer, intent(in) :: XX,YY
   ! Passed arrays
   real, intent(in), dimension(XX,YY,9)    :: dist ! distance to neighbor cells (m), 9 = cell area (m2)
   real, intent(in), dimension(XX,YY)      :: z ! depth to water table (m)
   real, intent(in), dimension(XX,YY)      :: latksat
                                              ! lateral saturated hydraulic conductivity (m/s)
   real, intent(in), dimension(XX,YY)      :: nexp ! latksat decay coefficient
   real, intent(in), dimension(XX,YY)      :: soldep ! soil depth (m)
   real, intent(in), dimension(XX,YY)      :: mfdexp ! slope exponent for flow partitioning
   real, intent(in), dimension(XX,YY,8)    :: SO8RT ! terrain slope in all directions (m/m)
   real, intent(in), dimension(XX,YY)      :: CWATAVAIL ! water available for routing (m)
   real, intent(in)                        :: SUBDT ! subsurface routing timestep (s)
   integer, intent(in)                     :: gradOpt
                                              ! hydraulic gradient: 1 = terrain slope (kinematic),
                                              ! 2 = water table (terrain + water-table-depth gradient)
   real, intent(inout), dimension(XX,YY)   :: QSUBDRY ! subsurface flow at domain boundary (m3/s)
   real, intent(inout)                     :: QSUBDRYT
                                              ! total subsurface flow outside of domain (m3/s)
   real, intent(out), dimension(XX,YY)     :: qsub ! net subsurface flow into cell (m3/s)

   ! Neighbor offsets, N first and clockwise (same order as SO8RT and dist)
   integer, parameter :: nbDI(8) = (/ 0, 1, 1, 1, 0, -1, -1, -1 /)
   integer, parameter :: nbDJ(8) = (/ 1, 1, 0, -1, -1, -1, 0, 1 /)

   ! Local variables
   integer :: i, j, k, ni, nj, kmax, klast
   real    :: betaK(8) ! hydraulic gradient to each neighbor (m/m)
   real    :: beta ! steepest gradient (m/m)
   real    :: hh, ksat, gamma, width, waterToRoute
   real    :: qqsub ! total outflow from cell (m3/s, negative = out of cell)
   real*8  :: wgt(8), wsum, qk, qrem
   logical :: onBoundary
   real*8, dimension(XX,YY) :: qsub_tmp, QSUBDRY_tmp ! temp trackers for fluxes (m3/s)

   qsub_tmp = 0.
   QSUBDRY_tmp = 0.
   qsub = 0.

   do j = 2, YY-1
      do i = 2, XX-1

         ! Hydraulic gradient to each neighbor
         do k = 1, 8
            if (gradOpt .eq. 2) then
               betaK(k) = SO8RT(i,j,k) - (z(i,j) - z(i+nbDI(k),j+nbDJ(k))) / dist(i,j,k)
            else
               betaK(k) = SO8RT(i,j,k)
            endif
         end do

         ! Steepest neighbor (first one wins ties, as in ROUTE_SUBSURFACE1)
         kmax = 1
         do k = 2, 8
            if (betaK(k) .gt. betaK(kmax)) kmax = k
         end do
         beta = betaK(kmax)

         if (dist(i,j,kmax) .le. 0) then
            call hydro_stop("In ROUTE_SUBSURFACE3() - dist(i,j,index) is <= zero ")
         endif
         if (soldep(i,j) .eq. 0) then
            call hydro_stop("In ROUTE_SUBSURFACE3() - soldep is = zero")
         endif

         if (beta .le. 0) cycle  ! no downslope neighbor
         if (beta .lt. 1E-20) beta = 1E-20

         ! Total outflow: identical to ROUTE_SUBSURFACE1
         hh = ( 1 - ( z(i,j) / soldep(i,j) ) ) ** nexp(i,j)
         ksat = latksat(i,j)
         if (hh .lt. 0.) then
            print *, "hsub<0 at gridcell...", i, j, hh, z(i,j), soldep(i,j)
            call hydro_stop("In ROUTE_SUBSURFACE3() - hsub<0 at gridcell ")
         endif
         width = dist(i,j,9) / dist(i,j,kmax)
         gamma = -1.0 * ( (width * ksat * soldep(i,j)) / nexp(i,j) ) * beta
         qqsub = gamma * hh
         waterToRoute = ABS(qqsub) / dist(i,j,9) * SUBDT
         if ( (qqsub .le. 0.0) .and. (CWATAVAIL(i,j) .lt. waterToRoute) ) THEN
            qqsub = -1.0 * CWATAVAIL(i,j) * dist(i,j,9) / SUBDT
         endif
         if (qqsub .ge. 0.) cycle

         ! Partition weights over downslope neighbors
         wgt = 0.d0
         klast = 0
         do k = 1, 8
            if (betaK(k) .gt. 0.) then
               wgt(k) = (dble(betaK(k)) / dble(beta)) ** dble(mfdexp(i,j))
               if (wgt(k) .gt. 0.d0) klast = k
            endif
         end do
         wsum = sum(wgt)

         ! Remove total from this cell
         qsub(i,j) = qsub(i,j) + qqsub

         ! Distribute to neighbors (qqsub negative = out of cell (i,j));
         ! the last receiver gets the remainder so the shares sum exactly to qqsub
         qrem = dble(qqsub)
         do k = 1, 8
            if (wgt(k) .le. 0.d0) cycle
            if (k .eq. klast) then
               qk = qrem
            else
               qk = dble(qqsub) * (wgt(k) / wsum)
               qrem = qrem - qk
            endif
            ni = i + nbDI(k)
            nj = j + nbDJ(k)
            qsub_tmp(ni,nj) = qsub_tmp(ni,nj) - qk

            ! Flux into a cell on the outer domain boundary leaves the domain
#ifdef MPP_LAND
            onBoundary = ((ni .eq. XX) .and. (right_id .lt. 0)) .or. &
                         ((ni .eq. 1)  .and. (left_id  .lt. 0)) .or. &
                         ((nj .eq. 1)  .and. (down_id  .lt. 0)) .or. &
                         ((nj .eq. YY) .and. (up_id    .lt. 0))
#else
            onBoundary = (ni .eq. 1) .or. (ni .eq. XX) .or. (nj .eq. 1) .or. (nj .eq. YY)
#endif
            if (onBoundary) then
               qsub_tmp(ni,nj)    = qsub_tmp(ni,nj) + qk
               QSUBDRY_tmp(ni,nj) = QSUBDRY_tmp(ni,nj) + qk
               QSUBDRYT           = QSUBDRYT + qk
            endif
         end do

      end do
   end do

#ifdef MPP_LAND
   ! Sum flux contributions written into halo cells back onto their owning tiles
   call MPP_LAND_COM_REAL8(qsub_tmp,XX,YY,1)
   call MPP_LAND_COM_REAL8(QSUBDRY_tmp,XX,YY,1)
#endif

   qsub = qsub + qsub_tmp
   QSUBDRY = QSUBDRY + QSUBDRY_tmp

#ifdef MPP_LAND
   call MPP_LAND_COM_REAL(qsub,XX,YY,99)
   call MPP_LAND_COM_REAL(QSUBDRY,XX,YY,99)
#endif

end subroutine ROUTE_SUBSURFACE3
!===================================================================================================
