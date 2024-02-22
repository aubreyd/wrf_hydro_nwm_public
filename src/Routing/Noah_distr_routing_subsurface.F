module module_subsurfaceRouting
contains
subroutine subsurfaceRouting ( subrt_data, subrt_static, subrt_input, subrt_output, &
                               NTRACER, WVTFLAG, PMOPT, NSUB) ! wt variables
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

!-----------------------------------------------------
! Needed for water tracer capability
!-----------------------------------------------------
    INTEGER,           INTENT(IN)                   :: NTRACER
    INTEGER, OPTIONAL, INTENT(IN)                   :: WVTFLAG
    INTEGER, OPTIONAL, INTENT(IN)                   :: PMOPT
    INTEGER, OPTIONAL, INTENT(IN), DIMENSION(subrt_static%NSOIL) :: NSUB

    !INTEGER :: SO8RT_D(IXRT,JXRT,3)
    !REAL :: SO8RT(IXRT,JXRT,8)
    !REAL, INTENT(IN)                          :: dist(ixrt,jxrt,9)
    !  -----local array ----------
    !REAL, DIMENSION(IXRT,JXRT)   :: ZWATTABLRT
    REAL, DIMENSION(subrt_static%ixrt,subrt_static%jxrt)   :: CWATAVAIL
    INTEGER, DIMENSION(subrt_static%ixrt,subrt_static%jxrt) :: SATLYRCHK

    print *, "ADCHECK: NTRACER", NTRACER
    if (present(WVTFLAG)) print *, "ADCHECK: WVTFLAG", WVTFLAG
    if (present(PMOPT)) print *, "ADCHECK: PMOPT", PMOPT
    if (present(NSUB)) print *, "ADCHECK: NSUB", NSUB

    CWATAVAIL = 0.
    CALL FINDZWAT(subrt_static%ixrt,subrt_static%jxrt,subrt_static%NSOIL, &
        subrt_data%grid_transform%smcrt, &
        subrt_data%grid_transform%smcmaxrt, &
        subrt_data%grid_transform%smcrefrt, &
        subrt_data%grid_transform%smcwltrt, subrt_data%properties%zsoil,SATLYRCHK,subrt_data%properties%zwattablrt, &
        CWATAVAIL,subrt_data%properties%sldpth)
#ifdef MPP_LAND
    call MPP_LAND_COM_REAL(subrt_data%properties%zwattablrt,subrt_static%ixrt,subrt_static%jxrt,99)
    call MPP_LAND_COM_REAL(CWATAVAIL,subrt_static%ixrt,subrt_static%jxrt,99)
    call MPP_LAND_COM_INTEGER(SATLYRCHK,subrt_static%ixrt,subrt_static%jxrt,99)
#endif


    !DJG Second, Call subsurface routing routine...
#ifdef HYDRO_D
    print *, "Beginning SUB_routing..."
    print *, "Routing method is ",subrt_static%rt_option, " direction."
#endif

    !!!! Find saturated layer depth...
    ! Loop through domain to determine sat. layers and assign wat tbl depth...
    !    and water available for subsfc routing (CWATAVAIL)...
    ! This subroutine returns: ZWATTABLRT, CWATAVAIL and SATLYRCHK

    if (present(WVTFLAG)) then
       print *, "ADCHECK: Calling SUBSFC_RTNG with tracers"
       CALL SUBSFC_RTNG( subrt_data, subrt_static, subrt_input, subrt_output, &
                         CWATAVAIL, SATLYRCHK &
                         ,NTRACER &
                         ,WVTFLAG             & !wt vars
                         ,PMOPT               & !wt vars
                         ,NSUB                & !wt vars
                         !,QSUBRT_TR_GAIN     = QSUBRT_TR_GAIN      & !wt vars
                         !,QSUBBDRYRT_TR_PASS = QSUBBDRYRT_TR_PASS  & !wt vars
                         !,QSUBBDRYTRT_TR_PASS= QSUBBDRYTRT_TR_PASS & !wt vars
                         !,SMCRT_TR_PASS      = SMCRT_TR_PASS       & !wt vars
                         !,SMCRT_SUB_PASS     = SMCRT_SUB_PASS      & !wt vars
                         !,INFXSUBRT_TR_PASS  = INFXSUBRT_TR_PASS   & !wt vars
                         !,SMCMAXRT_SUB_PASS  = SMCMAXRT_SUB_PASS   & !wt vars
                       )
    else
       print *, "ADCHECK: Calling SUBSFC_RTNG without tracers"
       CALL SUBSFC_RTNG( subrt_data, subrt_static, subrt_input, subrt_output, &
                         CWATAVAIL, SATLYRCHK &
                         ,NTRACER)
    end if

#ifdef HYDRO_D
    print *, "SUBROUTE routing called and returned..."
#endif

end subroutine subsurfaceRouting

!DJG ------------------------------------------------
!DJG   SUBROUTINE SUBSFC_RTNG
!DJG ------------------------------------------------

SUBROUTINE SUBSFC_RTNG(subrt_data, subrt_static, subrt_input, subrt_output, CWATAVAIL, SATLYRCHK &
                       ,NTRACER, WVTFLAG, PMOPT, NSUB       & !wt vars
                       !,QSUBRT_TR_GAIN, QSUBBDRYRT_TR_PASS  & !wt vars
                       !,QSUBBDRYTRT_TR_PASS, SMCRT_TR_PASS  & !wt vars
                       !,SMCRT_SUB_PASS, INFXSUBRT_TR_PASS   & !wt vars
                       !,SMCMAXRT_SUB_PASS                   & !wt vars
                       )

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

!-------------------------------------------------
! Needed for water tracer capability
!-------------------------------------------------
    INTEGER,          INTENT(IN)                             :: NTRACER
    INTEGER,OPTIONAL, INTENT(IN)                             :: WVTFLAG
    INTEGER,OPTIONAL, INTENT(IN)                             :: PMOPT
    INTEGER,OPTIONAL, INTENT(IN),   DIMENSION(subrt_static%nsoil)         :: NSUB
    !REAL,   OPTIONAL, INTENT(OUT),  DIMENSION(subrt_static%ixrt,subrt_static%jxrt,NTRACER) :: QSUBRT_TR_GAIN
    !REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(subrt_static%ixrt,subrt_static%jxrt,NTRACER) :: QSUBBDRYRT_TR_PASS
    !REAL,   OPTIONAL, INTENT(INOUT)                          :: QSUBBDRYTRT_TR_PASS
    !REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(subrt_static%ixrt,subrt_static%jxrt,NTRACER) :: SMCRT_TR_PASS
    !REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(subrt_static%ixrt,subrt_static%jxrt,NTRACER) :: SMCRT_SUB_PASS
    !REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(subrt_static%ixrt,subrt_static%jxrt,NTRACER) :: INFXSUBRT_TR_PASS
    !REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(subrt_static%ixrt,subrt_static%jxrt,NTRACER) :: SMCMAXRT_SUB_PASS
    !-----local vars------
    REAL,   DIMENSION(subrt_static%ixrt,subrt_static%jxrt,NTRACER)  :: QSUBRT_SUB_GAIN
    REAL,   DIMENSION(NTRACER)            :: SUBFLO_TR
    REAL,   DIMENSION(NTRACER)            :: SUBFLO_SUB
    REAL                                  :: PROPOR_TR
    INTEGER                               :: KSUB, KSOIL
    INTEGER                               :: ISUB
    REAL                                  :: SUBFLO_INIT, SUBFLO_EXCESS
    INTEGER                               :: INFXSUBRT_TR_FILL
!-------------------------------------------------

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

    if (present(WVTFLAG))  subflo_tr = 0.0

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
    print *, "calling subsurface routing subroutine...Opt. ",subrt_static%rt_option
#endif

    ! ADCHANGE: IMPORTANT!
    ! 2D subsurface option currently has bug so forcing to option 1 in this routine to
    ! allow users to still have option to use 2d overland (both are controlled by same
    ! rt_option flag). Remove this hard-coded option when rt_option=2 is fixed for subsurface.
    !     if(subrt_static%rt_option .eq. 1) then
    if (present(WVTFLAG)) then
       print *, "ADCHECK: Before ROUTE_SUBSURFACE1 tracers"
       CALL ROUTE_SUBSURFACE1(subrt_data%properties%distance_to_neighbor, &
                           subrt_data%properties%zwattablrt, &
                           subrt_data%properties%lksatrt, subrt_data%properties%nexprt, &
                           subrt_data%properties%soldeprt, &
                           subrt_static%IXRT, subrt_static%JXRT, &
                           subrt_data%properties%surface_slope, &
                           subrt_data%properties%max_surface_slope_index, &
                           CWATAVAIL, SUBDT, &
                           subrt_data%state%qsubbdryrt, subrt_data%state%qsubbdrytrt, &
                           subrt_data%state%qsubrt, &
                           NTRACER, subrt_static%nsoil &
                          ,WVTFLAG          &
                          ,PMOPT            &
                          ,NSUB             &
                          ,subrt_data%properties%sldpth           &
                          ,SATLYRCHK        &
                          ,subrt_data%state%QSUBRT_SUB  &
                          ,subrt_data%state%QSUBRT_TR   &
                          ,subrt_data%state%QSUBBDRYRT_TR  &
                          ,subrt_data%state%QSUBBDRYTRT_TR &
                          ,subrt_data%grid_transform%SMCRT_TR    &
                          ,subrt_data%grid_transform%SMCRT_SUB   &
                          )
    else
       print *, "ADCHECK: Before ROUTE_SUBSURFACE1 no tracers"
       CALL ROUTE_SUBSURFACE1(subrt_data%properties%distance_to_neighbor, &
                           subrt_data%properties%zwattablrt, &
                           subrt_data%properties%lksatrt, subrt_data%properties%nexprt, &
                           subrt_data%properties%soldeprt, &
                           subrt_static%IXRT, subrt_static%JXRT, &
                           subrt_data%properties%surface_slope, &
                           subrt_data%properties%max_surface_slope_index, &
                           CWATAVAIL, SUBDT, &
                           subrt_data%state%qsubbdryrt, subrt_data%state%qsubbdrytrt, &
                           subrt_data%state%qsubrt, &
                           NTRACER, subrt_static%nsoil)
    end if

    !     else
    !        CALL ROUTE_SUBSURFACE2(subrt_data%properties%distance_to_neighbor,subrt_data%properties%zwattablrt, &
        !               subrt_data%state%qsubrt, subrt_data%properties%surface_slope_x, subrt_data%properties%surface_slope_y, &
        !               subrt_data%properties%lksatrt, subrt_data%properties%soldeprt, subrt_static%IXRT, subrt_static%JXRT, subrt_data%state%qsubbdryrt, subrt_data%state%qsubbdrytrt, &
        !               CWATAVAIL,SUBDT)
    !     end if

#ifdef HYDRO_D
    write(6,*) "finish calling ROUTE_SUBSURFACE ", subrt_static%rt_option
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

            print *, "ADCHECK: Before SUBFLO_SUB"
            !-----wt cacl-----
            if (present(WVTFLAG)) then
               INFXSUBRT_TR_FILL = 0         !initialize INFXSUBRT_TR_FILL
               SUBFLO_INIT       = SUBFLO    !save this value for tracer purpose
               SUBFLO_EXCESS     = 0.0
               SUBFLO_SUB        = 0.0
               SUBFLO_TR         = 0.0
               do KSOIL  = 1, subrt_static%nsoil 
                  do KSUB = 1, NSUB(KSOIL)
                     ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
                     SUBFLO_SUB(ISUB) = subrt_data%state%QSUBRT_SUB(I,J,ISUB) / subrt_data%properties%distance_to_neighbor(i,j,9) * SUBDT
                        ! -> can only be >0 no matter what SUBLO is
                     SUBFLO_TR(ISUB)  = subrt_data%state%QSUBRT_TR(I,J,ISUB) / subrt_data%properties%distance_to_neighbor(i,j,9) * SUBDT
                  end do
               end do
            end if
            !-----------------

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
                   !-----wt calc-----
                   if (present(WVTFLAG)) then
                      INFXSUBRT_TR_FILL = 1       !to mark that all layers saturated, and to add remaining subflo to infilt. excess
                      SUBFLO_EXCESS = SUBFLO
                   end if
                   !-----------------
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

!-------------different pathways for tracer transport in subsurface routing----------
            if (present(WVTFLAG)) then
               if (present(PMOPT)) then
                  if (PMOPT .eq. 15 .or. PMOPT .eq. 25) then
                     ! add bottom-to-top flow
                     call implement_tracerflow_bottom2top(I,J, subrt_static%nsoil, &
                        NTRACER, subrt_data%properties%sldpth, NSUB, &
                        subrt_data%grid_transform%SMCRT(I,J,:), &
                        subrt_data%grid_transform%SMCRT_SUB(I,J,:), subrt_data%grid_transform%SMCRT_TR(I,J,:),        &
                        SUBFLO_INIT, SUBFLO_SUB, SUBFLO_TR, &
                        INFXSUBRT_TR_FILL,subrt_input%infiltration_excess_tr(I,J))
                  else
                     if (PMOPT .eq. 16 .or. PMOPT .eq. 26) then
                        ! add laminar flow
                        call implement_tracerflow_laminar(I,J, subrt_static%nsoil, &
                           NTRACER, subrt_data%properties%sldpth, NSUB, &
                           subrt_data%grid_transform%SMCRT(I,J,:), &
                           subrt_data%grid_transform%SMCRT_SUB(I,J,:), subrt_data%grid_transform%SMCRT_TR(I,J,:),        &
                           SUBFLO_INIT, SUBFLO_SUB, SUBFLO_TR, &
                           INFXSUBRT_TR_FILL,SUBFLO_EXCESS,subrt_input%infiltration_excess_tr(I,J))
                     else
                        if (PMOPT .eq. 17 .or. PMOPT .eq. 27) then
                           ! add top-to-bottom flow
                           call implement_tracerflow_top2bottom(I,J, subrt_static%nsoil, &
                              NTRACER, subrt_data%properties%sldpth, NSUB, &
                              subrt_data%grid_transform%SMCRT(I,J,:), &
                              subrt_data%grid_transform%SMCRT_SUB(I,J,:), subrt_data%grid_transform%SMCRT_TR(I,J,:),        &
                              SUBFLO_INIT, SUBFLO_SUB, SUBFLO_TR, &
                              INFXSUBRT_TR_FILL,SUBFLO_EXCESS,subrt_input%infiltration_excess_tr(I,J))
                        else
                        ! add default flow (fill water gaps)
                           call implement_tracerflow_fill(I,J, subrt_static%nsoil, &
                              NTRACER, subrt_data%properties%sldpth, NSUB, &
                              subrt_data%grid_transform%SMCRT(I,J,:), &
                              subrt_data%grid_transform%SMCRT_SUB(I,J,:), subrt_data%grid_transform%SMCRT_TR(I,J,:),        &
                              SUBFLO_INIT, SUBFLO_SUB, SUBFLO_TR, &
                              INFXSUBRT_TR_FILL,subrt_input%infiltration_excess_tr(I,J))
                        end if
                     end if
                  end if
               end if
            end if
!------------------------------------------------------------------------------------

        END DO        ! END DO X dim
    END DO          ! END DO Y dim
    !!!! End loop through subsurface routing domain...

#ifdef MPP_LAND
    do i = 1, subrt_static%nsoil
        call MPP_LAND_COM_REAL(subrt_data%grid_transform%smcrt(:,:,i),subrt_static%IXRT,subrt_static%JXRT,99)
    end DO
    !-----wt calc-----
    if (present(WVTFLAG)) then
       do i = 1, NTRACER
          call MPP_LAND_COM_REAL(subrt_data%grid_transform%SMCRT_TR(:,:,i), subrt_static%IXRT, subrt_static%JXRT, 99)
          call MPP_LAND_COM_REAL(subrt_data%grid_transform%SMCRT_SUB(:,:,i), subrt_static%IXRT, subrt_static%JXRT, 99)
       end do
    end if
    !----------------
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
SUBROUTINE FINDZWAT(IXRT,JXRT,NSOIL,SMCRT,SMCMAXRT,SMCREFRT, &
        SMCWLTRT,ZSOIL,SATLYRCHK,ZWATTABLRT,CWATAVAIL,&
        SLDPTH)

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

    !DJG Local Variables
    INTEGER :: KK,i,j

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
                             CWATAVAIL, SUBDT, QSUBDRY, QSUBDRYT, qsub &
                            ,NTRACER, NSOIL        & !wt vars
                            ,WVTFLAG, PMOPT, NSUB, SLDPTH, SATLYRCHK   & !wt vars
                            ,QSUBRT_SUB_GAIN,QSUBRT_TR_GAIN         & !wt vars
                            ,QSUBBDRYRT_TR_PASS, QSUBBDRYTRT_TR_PASS & !wt vars
                            ,SMCRT_TR_PASS, SMCRT_SUB_PASS           & !wt vars
                            )

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
   real, intent(inout), dimension(XX,YY)   :: QSUBDRY ! subsurface flow at domain boundary (m3/s)
   real, intent(inout)                     :: QSUBDRYT
                                              ! total subsurface flow outside of domain (m3/s)
   real, intent(out), dimension(XX,YY)     :: qsub ! subsurface flow outside of cell (m3/s)
!-------------------------------------------------
! Needed for water tracer capability
!-------------------------------------------------
   INTEGER,          INTENT(IN)                             :: NTRACER
   INTEGER,          INTENT(IN)                             :: NSOIL
   INTEGER,OPTIONAL, INTENT(IN)                             :: WVTFLAG
   INTEGER,OPTIONAL, INTENT(IN)                             :: PMOPT
   INTEGER,OPTIONAL, INTENT(IN),   DIMENSION(NSOIL)         :: NSUB
   REAL,   OPTIONAL, INTENT(IN),   DIMENSION(NSOIL)         :: SLDPTH
   INTEGER,OPTIONAL, INTENT(IN),   DIMENSION(XX,YY)         :: SATLYRCHK
   REAL,   OPTIONAL, INTENT(OUT),  DIMENSION(XX,YY,NTRACER) :: QSUBRT_SUB_GAIN
   REAL,   OPTIONAL, INTENT(OUT),  DIMENSION(XX,YY,NTRACER) :: QSUBRT_TR_GAIN
   REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(XX,YY,NTRACER) :: QSUBBDRYRT_TR_PASS
   REAL,   OPTIONAL, INTENT(INOUT)                          :: QSUBBDRYTRT_TR_PASS
   REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(XX,YY,NTRACER) :: SMCRT_TR_PASS
   REAL,   OPTIONAL, INTENT(INOUT),DIMENSION(XX,YY,NTRACER) :: SMCRT_SUB_PASS
   !-----local vars----
   INTEGER    :: KSUB, ISUB, KSOIL
   REAL       :: PROPOR_TR
   REAL, DIMENSION(NTRACER)       :: qqsub_tr_lose, qqsub_sub_lose
!-------------------------------------------------
   ! Local variables
   integer :: i, j ! indices for local loops
   integer :: IXX0, JYY0, index ! i, j, and direction index for steepest slope neighbor
   real    :: beta ! total head slope (m/m)
   real    :: hh ! interim variable for subsurface solution per DHSVM (dimensionless)
   real    :: gamma ! interim variable for subsurface solution per DHSVM (m3/s)
   real    :: ksat ! saturated hydraulic conductivity (m/s)
   real    :: qqsub ! net subsurface flow out of cell (m3/s?)
   real    :: waterToRoute ! total water to route from cell (m)
   ! Local arrays
   real*8, DIMENSION(XX,YY) :: qsub_tmp, QSUBDRY_tmp ! temp trackers for fluxes (m3/s)
   ! Local parameters
   real :: tmp_dist(9)

   print *, "ADCHECK: Start ROUTE_SUBSURFACE1"

   ! Initialize temp variables
   qsub_tmp = 0.
   QSUBDRY_tmp = 0.
  !----wt cal----
  if (present(WVTFLAG)) then
     qsubrt_sub_gain = 0.0
     qsubrt_tr_gain  = 0.0
     qqsub_tr_lose   = 0.0
     qqsub_sub_lose  = 0.0
  end if     
  !--------------

#ifdef HYDRO_D
   write(6,*) "call subsurface routing xx= , yy =", yy, xx
#endif

   do j=2,YY-1 ! start j loop

      do i=2,XX-1 ! start i loop

         if (i.ge.2.AND.i.le.XX-1.AND.j.ge.2.AND.j.le.YY-1) then !if grdcl chk
            ! check for boundary grid point?

            ! Set initial guess values for steepest slope neighbor based on terrain
            IXX0 = SO8RT_D(i,j,1)
            JYY0 = SO8RT_D(i,j,2)
            index = SO8RT_D(i,j,3)
            !beta = -1.0

            ! ADCHANGE: Always call routine to check for steepest elev+head slope
            ! Updates IXX0, JYY0, index, and beta for steepest slope neighbor.
            ! Note that beta will be -1 if steepest slope cannot be found.
            tmp_dist = dist(i,j,:)
            call GETSUB8(i, j, XX, YY, Z, so8rt, tmp_dist, &
                         IXX0, JYY0, index, beta)

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
               gamma = -1.0 * ( (dist(i,j,index) * ksat * soldep(i,j)) / nexp(i,j) ) * beta
               qqsub = gamma * hh

               ! AD_NOTE: Moved this water available constraint from outside qsub calc loop
               ! to inside to better account for adjustments to adjacent cells.
               ! Calculate total water to route (where dist(i,j,9) is cell area):
               waterToRoute = ABS(qqsub) / dist(i,j,9) * SUBDT
               if ( (qqsub .le. 0.0) .and. (CWATAVAIL(i,j) .lt. waterToRoute) ) THEN
                  qqsub = -1.0 * CWATAVAIL(i,j) * dist(i,j,9) / SUBDT
               endif

               !-----wt calc-----
               if (present(WVTFLAG)) then
                  print *, "ADCHECK: In ROUTE_SUBSURFACE1 tracer1"
                  if (qqsub .le. 1.0e-12 .and. SATLYRCHK(i,j) .lt. NSOIL+1) then
                     call calculate_qsubout_tr(i,j,SATLYRCHK(i,j),NSOIL,NTRACER,PMOPT,NSUB,SLDPTH,&
                                            qqsub,dist(i,j,9),SUBDT,qqsub_tr_lose,qqsub_sub_lose,  &
                                            SMCRT_TR_PASS(i,j,:), SMCRT_SUB_PASS(i,j,:))
                     do KSOIL  = 1, NSOIL
                        do KSUB = 1, NSUB(KSOIL)
                           ISUB  = sum(NSUB(1:KSOIL))-NSUB(KSOIL)+KSUB
                           qsubrt_tr_gain(ixx0,jyy0,ISUB)  = qsubrt_tr_gain(ixx0,jyy0,ISUB)  - qqsub_tr_lose(ISUB)   ! downstream pixel water gain
                           qsubrt_sub_gain(ixx0,jyy0,ISUB) = qsubrt_sub_gain(ixx0,jyy0,ISUB) - qqsub_sub_lose(ISUB)  ! downstream pixel water gain
                        end do
                     end do
!----------------------
                  else
                     qqsub_tr_lose   = 0.0
                     qqsub_sub_lose  = 0.0
                  end if
               end if
               !-----------------

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

                  !-----wt calc-----
                  if (present(WVTFLAG)) then
                     print *, "ADCHECK: In ROUTE_SUBSURFACE1 tracer2"
                     do KSOIL = 1, NSOIL
                        do KSUB = 1, NSUB(KSOIL)
                           ISUB = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
                           qsubrt_tr_gain(ixx0,jyy0, ISUB) = qsubrt_tr_gain(ixx0,jyy0,ISUB) + qqsub_tr_lose(ISUB)
                           qsubrt_sub_gain(ixx0,jyy0,ISUB) = qsubrt_sub_gain(ixx0,jyy0,ISUB) + qqsub_sub_lose(ISUB)
                           QSUBBDRYRT_TR_PASS(ixx0,jyy0,ISUB) = QSUBBDRYRT_TR_PASS(ixx0,jyy0,ISUB) + qqsub_tr_lose(ISUB)
                           QSUBBDRYTRT_TR_PASS = QSUBBDRYTRT_TR_PASS + qqsub_tr_lose(ISUB)
                        end do
                     end do
                  end if
                  !-----------------

               endif   ! end if for edge check

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

  print *, "ADCHECK: End ROUTE_SUBSURFACE1"

#ifdef MPP_LAND
   call MPP_LAND_COM_REAL(qsub,XX,YY,99)
   call MPP_LAND_COM_REAL(QSUBDRY,XX,YY,99)
   !-----wt cal-----
   if (present(WVTFLAG)) then
      call MPP_LAND_COM_REAL(qsubrt_tr_gain,XX,YY,99)
      call MPP_LAND_COM_REAL(qsubrt_sub_gain,XX,YY,99)
      call MPP_LAND_COM_REAL(QSUBBDRYRT_TR_PASS,XX,YY,99)
   end if
   !----------------
#endif

   return

end subroutine ROUTE_SUBSURFACE1

!===================================================================================================



!DJG ----------------------------------------------------------------
!DJG     SUBROUTINE ROUTE_SUBSURFACE2
!DJG ----------------------------------------------------------------

SUBROUTINE ROUTE_SUBSURFACE2(                                 &
        dist,z,qsub,sox,soy,                                   &
        latksat,soldep,XX,YY,QSUBDRY,QSUBDRYT,CWATAVAIL,   &
        SUBDT)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  Subroutine to route subsurface flow through the watershed
    !DJG ----------------------------------------------------------------
    !
    !  Called from: main.f (Noah_router_driver)
    !
    !  Returns: qsub=DQSUB   which in turn becomes SUBFLO in head calc.
    !
    !  Created:    D. Gochis                           3/27/03
    !              Adaptded from Wigmosta, 1994
    !
    !  Modified:   D. Gochis                           1/05/04
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MPP_LAND
    use module_mpp_land, only: left_id,down_id,right_id,&
        up_id,mpp_land_com_real,MPP_LAND_UB_COM, &
        MPP_LAND_LR_COM,mpp_land_com_integer
#endif
    use module_hydro_stop, only:HYDRO_stop

    IMPLICIT NONE


    !! Declare Passed variables

    INTEGER, INTENT(IN) :: XX,YY

    !! Declare passed arrays

    REAL, INTENT(IN), DIMENSION(XX,YY) :: z
    REAL, INTENT(IN), DIMENSION(XX,YY) :: sox
    REAL, INTENT(IN), DIMENSION(XX,YY) :: soy
    REAL, INTENT(IN), DIMENSION(XX,YY) :: latksat
    REAL, INTENT(IN), DIMENSION(XX,YY) :: CWATAVAIL
    REAL, INTENT(IN), DIMENSION(XX,YY) :: soldep
    REAL, INTENT(OUT), DIMENSION(XX,YY) :: qsub
    REAL, INTENT(INOUT), DIMENSION(XX,YY) :: QSUBDRY
    REAL, INTENT(INOUT) :: QSUBDRYT
    REAL, INTENT(IN) :: SUBDT
    real, intent(in), dimension(xx,yy,9) :: dist

    !!! Declare Local Variables

    REAL :: dzdx,dzdy,beta,gamma
    REAL :: qqsub,hh,ksat, gsize

    INTEGER :: i,j
    !!! Initialize variables
    REAL, PARAMETER :: nexp=1.0      ! local power law exponent
    qsub = 0.                        ! initialize flux = 0. !DJG 5 May 2014

    !yw        soldep = 2.


    ! Begin Subsurface routing

    !!! Loop to route water in x-direction
    do j=1,YY
        do i=1,XX
            ! check for boundary grid point?
            if (i.eq.XX) GOTO 998
            gsize = dist(i,j,3)

            dzdx= (z(i,j) - z(i+1,j))/gsize
            beta=sox(i,j) + dzdx + 1E-30
            if (abs(beta) .lt. 1E-20) beta=1E-20
            if (beta.lt.0) then
                !yw            hh=(1-(z(i+1,j)/soldep(i,j)))**nexp
                hh=(1-(z(i+1,j)/soldep(i+1,j)))**nexp
                ! Change later to use mean Ksat of two cells
                ksat=latksat(i+1,j)
            else
                hh=(1-(z(i,j)/soldep(i,j)))**nexp
                ksat=latksat(i,j)
            end if

            if (hh .lt. 0.) then
                print *, "hsub<0 at gridcell...", i,j,hh,z(i+1,j),z(i,j), &
                    soldep(i,j),nexp
                call hydro_stop("In ROUTE_SUBSURFACE2() - hsub<0 at gridcell")
            end if

            !Err. tan slope          gamma=-1.*((gsize*ksat*soldep(i,j))/nexp)*tan(beta)
            !AD_CHANGE: beta is already a slope so no tan (consistent with ROUTE_SUBSURFACE1)
            gamma=-1.*((gsize*ksat*soldep(i,j))/nexp)*beta
            !DJG lacks tan(beta) of original Wigmosta version          gamma=-1.*((gsize*ksat*soldep(i,j))/nexp)*beta

            qqsub = gamma * hh
            qsub(i,j) = qsub(i,j) + qqsub
            qsub(i+1,j) = qsub(i+1,j) - qqsub

            ! Boundary adjustments
#ifdef MPP_LAND
            if ((i.eq.1).AND.(beta.lt.0.).and.(left_id.lt.0)) then
#else
                if ((i.eq.1).AND.(beta.lt.0.)) then
#endif
                    qsub(i,j) = qsub(i,j) - qqsub
                    QSUBDRY(i,j) = QSUBDRY(i,j) - qqsub
                    QSUBDRYT = QSUBDRYT - qqsub
#ifdef MPP_LAND
                else if ((i.eq.(xx-1)).AND.(beta.gt.0.) &
                        .and.(right_id.lt.0) ) then
#else
                else if ((i.eq.(xx-1)).AND.(beta.gt.0.)) then
#endif
                    qsub(i+1,j) = qsub(i+1,j) + qqsub
                    QSUBDRY(i+1,j) = QSUBDRY(i+1,j) + qqsub
                    QSUBDRYT = QSUBDRYT + qqsub
                end if

                998       continue

                !! End loop to route sfc water in x-direction
        end do
    end do

#ifdef MPP_LAND
    call MPP_LAND_LR_COM(qsub,XX,YY,99)
    call MPP_LAND_LR_COM(QSUBDRY,XX,YY,99)
#endif


    !!! Loop to route water in y-direction
    do j=1,YY
        do i=1,XX
            ! check for boundary grid point?
            if (j.eq.YY) GOTO 999
            gsize = dist(i,j,1)

            dzdy= (z(i,j) - z(i,j+1))/gsize
            beta=soy(i,j) + dzdy + 1E-30
            if (abs(beta) .lt. 1E-20) beta=1E-20
            if (beta.lt.0) then
                !yw            hh=(1-(z(i,j+1)/soldep(i,j)))**nexp
                hh=(1-(z(i,j+1)/soldep(i,j+1)))**nexp
                ksat=latksat(i,j+1)
            else
                hh=(1-(z(i,j)/soldep(i,j)))**nexp
                ksat=latksat(i,j)
            end if

            if (hh .lt. 0.) GOTO 999

            !Err. tan slope          gamma=-1.*((gsize*ksat*soldep(i,j))/nexp)*tan(beta)
            gamma=-1.*((gsize*ksat*soldep(i,j))/nexp)*beta

            qqsub = gamma * hh
            qsub(i,j) = qsub(i,j) + qqsub
            qsub(i,j+1) = qsub(i,j+1) - qqsub

            ! Boundary adjustments

#ifdef MPP_LAND
            if ((j.eq.1).AND.(beta.lt.0.).and.(down_id.lt.0)) then
#else
            if ((j.eq.1).AND.(beta.lt.0.)) then
#endif
                qsub(i,j) = qsub(i,j) - qqsub
                QSUBDRY(i,j) = QSUBDRY(i,j) - qqsub
                QSUBDRYT = QSUBDRYT - qqsub
#ifdef MPP_LAND
            else if ((j.eq.(yy-1)).AND.(beta.gt.0.)  &
                .and. (up_id.lt.0) ) then
#else
            else if ((j.eq.(yy-1)).AND.(beta.gt.0.)) then
#endif
                qsub(i,j+1) = qsub(i,j+1) + qqsub
                QSUBDRY(i,j+1) = QSUBDRY(i,j+1) + qqsub
                QSUBDRYT = QSUBDRYT + qqsub
            end if

            999       continue

                !! End loop to route sfc water in y-direction
        end do
    end do

#ifdef MPP_LAND
    call MPP_LAND_UB_COM(qsub,XX,YY,99)
    call MPP_LAND_UB_COM(QSUBDRY,XX,YY,99)
#endif

    return
    !DJG------------------------------------------------------------
end subroutine ROUTE_SUBSURFACE2
!DJG------------------------------------------------------------

!------------------------------------------------------
! Huancui Hu: calculate_wt_ratio
subroutine calculate_wt_ratio (a,b,c)

  implicit none
  real,    intent(in)     :: a
  real,    intent(in)     :: b
  real,    intent(inout)  :: c

  if (b .gt. 1.0E-10) then
    c = min(1.0,max(a/b, 0.0))
  else
    c = 0.0
  end if

end subroutine
!------------------------------------------------------

!------------------------------------------------------
! Huancui Hu: calculate_wt_ratio
subroutine calculate_wt_ratio_pm (a,b,c)

  implicit none
  real,    intent(in)     :: a
  real,    intent(in)     :: b
  real,    intent(inout)  :: c
  !----local vars----
  real      :: d,b_new

  d = b-a
  if (d .gt. 1.0E-10) then
    b_new = a+d*0.5
    if (b_new .gt. 1.0E-8) then
      c = min(1.0,max(a/b_new, 0.0))
    else
      c = 0.0
    end if
  else
    c = 0.0
  end if

end subroutine
!-----------------------------------------------------

!-----------------------------------------------------
! Huancui Hu: implement_tracerflow_bottom2top
subroutine implement_tracerflow_bottom2top(I,J, NSOIL, NTRACER, SLDPTH, NSUB, &
                                     smc, smc_sub, smc_tr,        &
                                     subflo_init, subflo_sub, subflo_tr, &
                                     infxsubrt_tr_fill, infxsubrt_tr)

  implicit none
  integer, intent(in)                       :: I,J
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout),dimension(NTRACER) :: smc_tr
  real,    intent(in)                       :: subflo_init
  real,    intent(inout),dimension(NTRACER) :: subflo_sub
  real,    intent(inout),dimension(NTRACER) :: subflo_tr
  integer, intent(in)                       :: infxsubrt_tr_fill
  real,    intent(inout)                    :: infxsubrt_tr
  !----local vars----
  real                      :: subflo_sub_sum
  real                      :: subflo_tr_sum
  integer                   :: KSOIL, KSUB, ISUB
  real                      :: PROPOR_TR
  real                      :: dsmc_tmp
  real                      :: smc_residual
  !---------------------

  !initialize vars
  dsmc_tmp       =  0.0
  subflo_sub_sum =  0.0
  subflo_tr_sum  =  0.0

     ! sum up nlayer sub and tracer flow
  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      subflo_sub_sum = subflo_sub_sum + subflo_sub(ISUB)     !unit:m
      subflo_tr_sum  = subflo_tr_sum  + subflo_tr(ISUB)
    end do
  end do

  do KSOIL = NSOIL, 1, -1
    do KSUB = NSUB(KSOIL), 1,-1
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      if (subflo_sub_sum .gt. 1.0e-6) then
        smc_sub(ISUB) = smc_sub(ISUB) + subflo_sub_sum/(SLDPTH(KSOIL)/NSUB(KSOIL))
        smc_tr(ISUB)  = smc_tr(ISUB)  + subflo_tr_sum/(SLDPTH(KSOIL)/NSUB(KSOIL))
        call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
        dsmc_tmp      = smc(KSOIL) - smc_sub(ISUB)
        if (dsmc_tmp .gt. 1.0e-8) then  !downward adjustment is needed
          subflo_sub_sum  = 0.0
          subflo_tr_sum   = 0.0
        else
          if (abs(dsmc_tmp) .le. 1.0e-8) then  !meaning soil moisture matches after adjustment
            smc_sub(ISUB)  = smc(KSOIL)
            smc_tr(ISUB)   = min(max(0.0,smc_tr(ISUB)), smc_sub(ISUB))
            subflo_sub_sum = 0.0
            subflo_tr_sum  = 0.0
          else !dsmc_tmp < 0
            smc_sub(ISUB)  = smc(KSOIL)
            smc_tr(ISUB)   = smc_sub(ISUB) * PROPOR_TR
            subflo_sub_sum = 0.0 - dsmc_tmp*(SLDPTH(KSOIL)/NSUB(KSOIL))   !becomes positive
            subflo_tr_sum  = subflo_sub_sum* PROPOR_TR
          end if
        end if
      else
        if (abs(subflo_sub_sum) .le. 1.0e-6) then
          subflo_sub_sum = 0.0
          subflo_tr_sum  = 0.0
        else !subflo_sub_sum < 0
          print *, "ERROR 1 in implement_tracerflow_bottom2top"
        end if
      end if
    end do
  end do

  ! now the adjustment due to upstream flow is done
  ! may have all saturaed layers
  if (subflo_sub_sum .gt. 1.0e-6) then
    if (infxsubrt_tr_fill  .eq. 1) then  !sanity check
      infxsubrt_tr    = infxsubrt_tr + subflo_tr_sum*1000.   !unit:mm
      subflo_sub_sum  = 0.0
      subflo_tr_sum   = 0.0
    else
      print *, "ERROR 3 in implement_tracerflow_bottom2top"
      print *, "I,J,subflo_sub_sum,infxsubrt_tr_fill,subflo_sub:",I,J,subflo_sub_sum,infxsubrt_tr_fill,subflo_sub
    end if
  else
    !----do soil moisture adjustment downward-----
    call adjust_smc_downward(I,J,NSOIL, NTRACER, SLDPTH, NSUB, smc, smc_sub, smc_tr, smc_residual)
    if (abs(smc_residual)  .gt. 1.0e-6) then
      print *, "ERROR 4 in implement_tracerflow_bottom2top"
    end if
  end if

end subroutine implement_tracerflow_bottom2top
!-----------------------------------------------------

!-----------------------------------------------------
! Huancui Hu: implement_tracerflow_laminar
subroutine implement_tracerflow_laminar(I,J, NSOIL, NTRACER, SLDPTH, NSUB, &
                                     smc, smc_sub, smc_tr,        &
                                     subflo_init, subflo_sub, subflo_tr, &
                                     infxsubrt_tr_fill, subflo_excess, infxsubrt_tr)

  implicit none
  integer, intent(in)                       :: I,J
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout),dimension(NTRACER) :: smc_tr
  real,    intent(in)                       :: subflo_init
  real,    intent(inout),dimension(NTRACER) :: subflo_sub
  real,    intent(inout),dimension(NTRACER) :: subflo_tr
  integer, intent(in)                       :: infxsubrt_tr_fill
  real,    intent(in)                       :: subflo_excess
  real,    intent(inout)                    :: infxsubrt_tr
  !-----local vars------
  real                      :: subflo_sub_sum
  real                      :: subflo_tr_sum
  integer                   :: KSOIL, KSUB, ISUB
  real                      :: PROPOR_TR
  real                      :: dsmc_tmp, dwater_avail
  real                      :: smc_residual
  real, dimension(NTRACER)  :: sldpth_sub
  !---------------------

!  print *, "FLAG 1 IN implement_tracerflow_laminar"
  dsmc_tmp = 0.0
  subflo_sub_sum = 0.0       !not used
  subflo_tr_sum  = 0.0       !not used
  sldpth_sub          = 0.0
  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      sldpth_sub(ISUB) = SLDPTH(KSOIL)/NSUB(KSOIL)
    end do
  end do

  ! add the subflo_sub upstreamflow to each layer
  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      smc_sub(ISUB)  = smc_sub(ISUB) + subflo_sub(ISUB)/(SLDPTH(KSOIL)/NSUB(KSOIL))     !unit:m3/m3
      smc_tr(ISUB)   = smc_tr(ISUB)  + subflo_tr(ISUB)/(SLDPTH(KSOIL)/NSUB(KSOIL))
    end do
  end do

  ! now start to adjust soil moisture from bottom to top first
  do KSOIL = NSOIL, 1, -1
    do KSUB = NSUB(KSOIL), 1, -1
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      dsmc_tmp = smc(KSOIL) - smc_sub(ISUB)
      if (dsmc_tmp .gt. 1.0e-6) then
        !needs downward flow
      else
        if (abs(dsmc_tmp) .lt. 1.0e-6) then
          smc_sub(ISUB) = smc(KSOIL)
          smc_tr(ISUB)  = min(max(0.0, smc_tr(ISUB)), smc_sub(ISUB))
        else  !dsmc_tmp < 0
          call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
          if (ISUB .gt. 1) then
            smc_sub(ISUB)   = smc(KSOIL)
            smc_tr(ISUB)    = smc_tr(ISUB) + dsmc_tmp*PROPOR_TR
            smc_sub(ISUB-1) = smc_sub(ISUB-1) - dsmc_tmp*sldpth_sub(ISUB)/sldpth_sub(ISUB-1)
            smc_tr(ISUB-1)  = smc_tr(ISUB-1)  - dsmc_tmp*sldpth_sub(ISUB)*PROPOR_TR/sldpth_sub(ISUB-1)
          else !ISUB = 1, meaning there is saturated flow
            if (infxsubrt_tr_fill .eq. 1) then
              infxsubrt_tr  = infxsubrt_tr  + subflo_excess*PROPOR_TR*1000.  !unit:mm
              smc_sub(ISUB) = smc_sub(ISUB) - subflo_excess/sldpth_sub(ISUB)
              smc_tr(ISUB)  = smc_tr(ISUB)  - subflo_excess*PROPOR_TR/sldpth_sub(ISUB)
            else ! it is also possible if some lower layer is not filled to smc yet              
!              print *, "ERROR 1 in implement_tracerflow_laminar"
!              print *, "I,J, KSOIL, ISUB, smc(KSOIL),smc_sub(ISUB),subflo_sub(ISUB),infxsubrt_tr_fill:",I,J, KSOIL, ISUB, smc(KSOIL),smc_sub(ISUB),subflo_sub(ISUB),infxsubrt_tr_fill
            end if
          end if
        end if
      end if
    end do
  end do

!  print *, "FLAG 3 IN implement_tracerflow_laminar"

  ! do soil moisture adjustment from top to bottom then
  call adjust_smc_downward(I,J,NSOIL,NTRACER, SLDPTH, NSUB, smc, smc_sub, smc_tr, smc_residual)
  if (abs(smc_residual) .gt. 1.0e-6) then
    print *, "ERROR 2 in implement_tracerflow_laminar"
  end if

end subroutine implement_tracerflow_laminar
!-----------------------------------------------------

!-----------------------------------------------------
! Huancui Hu: implement_tracerflow_up2bottom 
subroutine implement_tracerflow_top2bottom(I,J, NSOIL, NTRACER, SLDPTH, NSUB, &
                                     smc, smc_sub, smc_tr,        &
                                     subflo_init, subflo_sub, subflo_tr, &
                                     infxsubrt_tr_fill,subflo_excess, infxsubrt_tr)

  implicit none
  integer, intent(in)                       :: I,J
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout),dimension(NTRACER) :: smc_tr
  real,    intent(in)                       :: subflo_init
  real,    intent(inout),dimension(NTRACER) :: subflo_sub
  real,    intent(inout),dimension(NTRACER) :: subflo_tr
  integer, intent(in)                       :: infxsubrt_tr_fill
  real,    intent(in)                       :: subflo_excess      !unit:m
  real,    intent(inout)                    :: infxsubrt_tr       !unit:mm
  !-----local vars------
  real                      :: subflo_sub_sum
  real                      :: subflo_tr_sum
  integer                   :: KSOIL, KSUB, ISUB
  real                      :: PROPOR_TR
  real                      :: dsmc_tmp, dwater_avail
  real                      :: smc_residual
  !---------------------

  dsmc_tmp = 0.0
  subflo_sub_sum = 0.0
  subflo_tr_sum  = 0.0


  ! sum up nlayer sub and tracer flow
  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      subflo_sub_sum = subflo_sub_sum + subflo_sub(ISUB)     !unit:m
      subflo_tr_sum  = subflo_tr_sum  + subflo_tr(ISUB)
    end do
  end do
  call calculate_wt_ratio(subflo_tr_sum, subflo_sub_sum, PROPOR_TR)

  ! add subflo_sub_sum and subflo_sub_tr to the top soil layer
  if (subflo_sub_sum .gt. 1.0e-6) then
    if (infxsubrt_tr_fill .eq. 1) then
       subflo_sub_sum  = subflo_sub_sum - subflo_excess
       subflo_tr_sum   = subflo_tr_sum  - subflo_excess*PROPOR_TR
       infxsubrt_tr    = infxsubrt_tr   + subflo_excess*PROPOR_TR*1000.0   !unit:mm
    else
    end if
    smc_sub(1)   = smc_sub(1) + subflo_sub_sum/(SLDPTH(1)/NSUB(1))  !unit:m3/m3
    smc_tr(1)    = smc_tr(1)  + subflo_tr_sum/(SLDPTH(1)/NSUB(1))
  else
    if (abs(subflo_sub_sum) .le. 1.0e-6) then
      subflo_sub_sum   = 0.0
      subflo_tr_sum    = 0.0
    else
      print *, "ERROR 1 in implement_tracerflow_up2bottom"
    end if
  end if

  ! do soil moisture adjustment from top to bottom then
  call adjust_smc_downward(I,J,NSOIL, NTRACER, SLDPTH, NSUB, smc, smc_sub, smc_tr, smc_residual)
  if (abs(smc_residual) .gt. 1.0e-6) then
    print *, "ERROR 2 in implement_tracerflow_up2bottom"
  end if

end subroutine implement_tracerflow_top2bottom
!-----------------------------------------------------

!-----------------------------------------------------
! Huancui Hu: implement_tracerflow_fill
subroutine implement_tracerflow_fill(I,J, NSOIL, NTRACER, SLDPTH, NSUB, &
                                     smc, smc_sub, smc_tr,        &
                                     subflo_init, subflo_sub, subflo_tr, &
                                     infxsubrt_tr_fill, infxsubrt_tr)

  implicit none
  integer, intent(in)                       :: I,J
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout),dimension(NTRACER) :: smc_tr
  real,    intent(in)                       :: subflo_init
  real,    intent(inout),dimension(NTRACER) :: subflo_sub
  real,    intent(inout),dimension(NTRACER) :: subflo_tr
  integer, intent(in)                       :: infxsubrt_tr_fill
  real,    intent(inout)                    :: infxsubrt_tr
  !-----local vars------
  real                      :: subflo_sub_sum
  real                      :: subflo_tr_sum
  integer                   :: KSOIL, KSUB, ISUB
  real                      :: PROPOR_TR
  real                      :: dsmc_tmp, dwater_avail
  real                      :: smc_residual
  !---------------------

  ! initialize vars
  dsmc_tmp = 0.0
  subflo_sub_sum = 0.0
  subflo_tr_sum  = 0.0

  ! sum up nlayer sub and tracer flow
  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      subflo_sub_sum = subflo_sub_sum + subflo_sub(ISUB)     !unit:m
      subflo_tr_sum  = subflo_tr_sum  + subflo_tr(ISUB)
    end do
  end do

  do KSOIL = NSOIL, 1, -1
    do KSUB = NSUB(KSOIL), 1, -1
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      if (subflo_sub_sum .gt. 1.0e-8) then
        dsmc_tmp      = smc(KSOIL) - smc_sub(ISUB)
        dwater_avail  = dsmc_tmp*SLDPTH(KSOIL)/NSUB(KSOIL)
        if (dsmc_tmp .gt. 1.0e-8) then
          if (dwater_avail .le. subflo_sub_sum ) then
            call calculate_wt_ratio(subflo_tr_sum, subflo_sub_sum, PROPOR_TR)
            subflo_sub_sum = subflo_sub_sum - dwater_avail
            subflo_tr_sum  = subflo_tr_sum  - dwater_avail*PROPOR_TR
            smc_sub(ISUB)  = smc(KSOIL)
            smc_tr(ISUB)   = smc_tr(KSOIL)  + dsmc_tmp*PROPOR_TR
          else  !gain from upstream flow is not sufficient to cover soil moisture change, need uplayer soil moiture to fill
            smc_sub(ISUB)  = smc_sub(ISUB)  + subflo_sub_sum/SLDPTH(KSOIL) !/NSUB(KSOIL)*NSUB(KSOIL) m3/m3
            smc_tr(ISUB)   = smc_tr(ISUB)   + subflo_tr_sum/SLDPTH(KSOIL)
            subflo_sub_sum = 0.0
            subflo_tr_sum  = 0.0
          end if
        else
          if (abs(dsmc_tmp) .le. 1.0e-8) then
            smc_sub(ISUB) = smc(KSOIL)
            smc_tr(ISUB)  = min(max(0.0,smc_tr(ISUB)), smc_sub(ISUB))
          else !dsmc_tmp<0,smc_sub need to lose some water
            call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
            smc_sub(ISUB)  = smc(KSOIL)
            smc_tr(ISUB)   = smc_tr(ISUB)   + dsmc_tmp*PROPOR_TR
            subflo_sub_sum = subflo_sub_sum - dwater_avail   !dwater_avail<0, subflo_sub_sum will increase
            subflo_tr_sum  = subflo_tr_sum  - dwater_avail*PROPOR_TR
          end if
        end if
      else
        if ((abs(subflo_sub_sum) .lt. 1.0e-8)) then
          subflo_sub_sum   = 0.0
          subflo_tr_sum    = 0.0
        else
          print *, "ERROR 2 in implement_tracerflow_fill"
        end if
      end if

    end do
  end do

  ! now the adjustment due to upstream flow is done
  ! may have all saturaed layers
  if (subflo_sub_sum .gt. 1.0e-6) then
    if (infxsubrt_tr_fill  .eq. 1) then  !sanity check
      infxsubrt_tr    = infxsubrt_tr + subflo_tr_sum*1000.   !unit:mm
      subflo_sub_sum  = 0.0
      subflo_tr_sum   = 0.0
    else
      print *, "ERROR 3 in implement_tracerflow_fill"
!----debug purpose----
      print *, "I, J, ISUB, subflo_sub_sum, subflo_tr_sum, smc, smc_sub, smc_tr, infxsubrt_tr_fill:"
      print *, I, J, ISUB, subflo_sub_sum, subflo_tr_sum, smc, smc_sub, smc_tr, infxsubrt_tr_fill
    end if
  else
    !----do soil moisture adjustment downward-----
    call adjust_smc_downward(I,J,NSOIL, NTRACER, SLDPTH, NSUB, smc, smc_sub, smc_tr, smc_residual)
    if (abs(smc_residual) .gt. 1.0e-6) then
      print *, "ERROR 4 in implement_tracerflow_fill"
    end if
  end if

end subroutine implement_tracerflow_fill
!------------------------------------------------------

!------------------------------------------------------
! Huancui Hu: adjust_smc_downward
subroutine adjust_smc_downward(I, J, NSOIL, NTRACER, SLDPTH, NSUB,  &
                               smc, smc_sub, smc_tr, smc_residual)

  implicit none
  integer, intent(in)                       :: I,J
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout),dimension(NTRACER) :: smc_tr
  real,    intent(out)                      :: smc_residual
  !----local vars-----
  INTEGER             :: KSOIL, KSUB, ISUB
  real                :: PROPOR_TR
  real                :: dsmc_tmp
  real, dimension(NTRACER)    :: sldpth_sub
  INTEGER, dimension(NTRACER) :: solve_flag
  real                        :: water_sub_extra, water_tr_extra

  smc_residual        = 0.0
  sldpth_sub          = 0.0
  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      sldpth_sub(ISUB) = SLDPTH(KSOIL)/NSUB(KSOIL)
    end do
  end do

  solve_flag      = 0
  water_sub_extra = 0.0
  water_tr_extra  = 0.0
  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      dsmc_tmp = smc(KSOIL) - smc_sub(ISUB)
      if (dsmc_tmp .gt. 1.0e-6) then
        !it is possible to have smc>smc_sub because smc can only be adjusted to reference soil moisture then use the lower levels for drainage
        solve_flag(ISUB) = 1
      else
        if (abs(dsmc_tmp) .le. 1.0e-6) then
          smc_sub(ISUB) = smc(KSOIL)
          smc_tr(ISUB)  = min(max(0.0,smc_tr(ISUB)), smc_sub(ISUB))
        else !dsmc_tmp <0
          call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
          smc_sub(ISUB)    = smc(KSOIL)
          smc_tr(ISUB)     = smc_tr(ISUB)  + dsmc_tmp*PROPOR_TR
          water_sub_extra  = water_sub_extra - dsmc_tmp*sldpth_sub(ISUB)              !unit:m
          water_tr_extra   = water_tr_extra  - dsmc_tmp*sldpth_sub(ISUB)*PROPOR_TR    !unit:m
        end if
      end if
    end do
  end do

  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      if (solve_flag(ISUB) .eq. 1) then
        if (water_sub_extra .lt. -1.0e-6) then
          print *, "ERROR 1 in adjust_smc_downward"
        else
          if (abs(water_sub_extra) .lt. 1.0e-6) then
            smc_sub(ISUB) = smc(KSOIL)
            smc_tr(ISUB)  = min(max(0.0,smc_tr(ISUB)), smc_sub(ISUB))
          else     !smc>smc_sub
            call calculate_wt_ratio(water_tr_extra, water_sub_extra, PROPOR_TR)
            dsmc_tmp       = smc(KSOIL) - smc_sub(ISUB)     !dsmc_tmp>0
            smc_sub(ISUB)  = smc(KSOIL)
            smc_tr(ISUB)   = smc_tr(ISUB) + dsmc_tmp*PROPOR_TR
            water_sub_extra = water_sub_extra - dsmc_tmp*sldpth_sub(ISUB)
            water_tr_extra  = water_tr_extra  - dsmc_tmp*sldpth_sub(ISUB)*PROPOR_TR
          end if
        end if
        solve_flag(ISUB)  = 0
      end if
    end do
  end do

  do KSOIL = 1, NSOIL
    do KSUB = 1, NSUB(KSOIL)
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      if (solve_flag(ISUB)  .ne. 0) then
        print *, "ERROR 2 in adjust_smc_downward"
      end if
    end do
  end do
  !-----calculate soil moisture residual-----
  call calculate_smc_residual(NSOIL, NTRACER, SLDPTH, NSUB, smc,smc_sub, smc_residual)

end subroutine adjust_smc_downward
!------------------------------------------------------  

!------------------------------------------------------  
! Huancui Hu: calculate_smc_residual
subroutine calculate_smc_residual(NSOIL, NTRACER, SLDPTH, NSUB, smc, smc_sub, smc_residual)
  implicit none
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout)                    :: smc_residual
  !----local vars-----
  integer                :: KSOIL, KSUB, ISUB

  smc_residual = 0.0
  do KSOIL = NSOIL, 1, -1
    do KSUB = NSUB(KSOIL), 1, -1
      ISUB  = sum(NSUB(1:KSOIL)) - NSUB(KSOIL) + KSUB
      smc_residual = smc_residual + abs(smc(KSOIL) - smc_sub(ISUB))*SLDPTH(KSOIL)/NSUB(KSOIL)  !unit:m
    end do
  end do

end subroutine calculate_smc_residual
!------------------------------------------------------

!------------------------------------------------------
! Huancui Hu: adjust_smc_tr
subroutine adjust_smc_tr(I,J,NSOIL, NTRACER, SLDPTH, NSUB, Ksat_ini,&
                         Ksat_now, smc, smc_sub, smc_tr,&
                         subflo_sub, subflo_tr)

  implicit none
  integer, intent(in)                       :: I,J
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  integer, intent(in)                       :: Ksat_ini
  integer, intent(in)                       :: Ksat_now     ! the level below Ksat_now is saturated
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout),dimension(NTRACER) :: smc_tr
  real,    intent(inout)                    :: subflo_sub
  real,    intent(inout)                    :: subflo_tr
  !-----local vars-----
  integer        :: K, KSUB, ISUB
  real           :: PROPOR_TR
  real           :: dsmc_tmp
  real           :: dsmc_acc
  real           :: dsmc_tr_acc
  integer,dimension(NTRACER) :: SOLVE_ORDER
  !--------------------

  if (Ksat_now .ge. Ksat_ini) then   ! only works if Ksat_now is below Kast_ini by default
    if (SUBFLO_SUB .gt. 1.0E-12) then   !  SUBFLO_SUB >0 while SUBLO<0?? correct, SUBFLO_SUB is always >0 because it only has the gain component
      call calculate_wt_ratio(subflo_tr, subflo_sub, PROPOR_TR)
!----for debug purpose----
if (I .eq. 93 .and. J .eq. 117) then
!  print *, "FLAG IN adjust_smc_tr01, subflo_tr, subflo_sub, PROPOR_TR:",subflo_tr, subflo_sub, PROPOR_TR
end if
!-------------------------

      do K = NSOIL, Ksat_ini, -1     ! not an error, do need to adjust each layer below Ksat_ini for tracer
        do KSUB = NSUB(K), 1, -1
          ISUB = sum(NSUB(1:K)-NSUB(K)+KSUB)
          dsmc_tmp = smc(K)-smc_sub(ISUB)    !dsmc_tmp can be +/-
!----for debug purpose----
if (I .eq. 93 .and. J .eq. 117) then
!  print *, "FLAG IN adjust_smc_tr02, K, ISUB, smc(K),smc_sub(ISUB),dsmc_tmp,subflo_sub:",K, ISUB, smc(K),smc_sub(ISUB),dsmc_tmp,subflo_sub
end if
!-------------------------
          if (subflo_sub .gt. 1.0E-10) then
            if (dsmc_tmp .gt. subflo_sub/SLDPTH(K)) then
              call calculate_wt_ratio(subflo_tr, subflo_sub, PROPOR_TR)
              smc_tr(ISUB)  = smc_tr(ISUB) + subflo_sub/SLDPTH(K)*PROPOR_TR
              smc_sub(ISUB) = smc_sub(ISUB)+ subflo_sub/SLDPTH(K)
              subflo_sub  = 0.0
              subflo_tr   = 0.0
            else
              if (dsmc_tmp .gt. 1.0E-10) then   !dsmc_tmp >0
                call calculate_wt_ratio(subflo_tr, subflo_sub, PROPOR_TR)
                smc_tr(ISUB)  = smc_tr(ISUB) + dsmc_tmp*PROPOR_TR
                smc_sub(ISUB) = smc(K)
                subflo_sub  = subflo_sub - dsmc_tmp*SLDPTH(K)/NSUB(K)
                subflo_tr   = subflo_tr  - dsmc_tmp*SLDPTH(K)/NSUB(K)*PROPOR_TR
              else 
                if (dsmc_tmp .lt. -1.0E-10) then
                  call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
                  smc_tr(ISUB)  = smc_tr(ISUB) + dsmc_tmp*PROPOR_TR
                  smc_sub(ISUB) = smc(K)
                  subflo_sub  = subflo_sub - dsmc_tmp*SLDPTH(K)/NSUB(K)
                  subflo_tr   = subflo_tr  - dsmc_tmp*SLDPTH(K)/NSUB(K)*PROPOR_TR
                end if
              end if
            end if
          end if
!----for debug purpose----
if (I .eq. 93 .and. J .eq. 117) then
!  print *, "FLAG IN adjust_smc_tr03, K, ISUB, smc(K),smc_sub(ISUB),dsmc_tmp,subflo_sub:",K, ISUB, smc(K),smc_sub(ISUB),dsmc_tmp,subflo_sub
end if
!-------------------------
        end do
      end do
    end if
  end if
  if (subflo_sub .gt. 1.0E-6) then
    print *, "ERROR in adjust_smc_tr"
    print *, "subflo_sub,i,j:", subflo_sub,i,j
  else
    subflo_sub = 0.0
  end if
  !-----adjust soil moisture downward/upward------
  dsmc_acc  = 0.0
  dsmc_tr_acc = 0.0
  SOLVE_ORDER = 0
  do K = Ksat_ini, NSOIL
    do KSUB = 1, NSUB(K)
      ISUB = sum(NSUB(1:K))-NSUB(K)+KSUB
      call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
      dsmc_tmp = smc_sub(ISUB) - smc(K)
      if (dsmc_tmp .lt. -1.E-12) then
        SOLVE_ORDER(ISUB) = 3
      else
        if (dsmc_tmp .gt. 1.E-12) then
          SOLVE_ORDER(ISUB) = 1
          dsmc_acc      = dsmc_acc + dsmc_tmp*SLDPTH(K)/NSUB(K)
          dsmc_tr_acc   = dsmc_tr_acc + dsmc_tmp*SLDPTH(K)/NSUB(K)*PROPOR_TR
          smc_tr(ISUB)  = smc_tr(ISUB) - dsmc_tmp*PROPOR_TR
          smc_sub(ISUB) = smc(K)
        end if
      end if
    end do
  end do
!----for debug purpose----
if (I .eq. 93 .and. J .eq. 117) then
!  print *, "FLAG IN adjust_smc_tr04, Ksat_ini, NSOIL, SOLVE_ORDER, smc,smc_sub,smc_tr:",Ksat_ini, NSOIL, SOLVE_ORDER, smc,smc_sub,smc_tr
end if
!-------------------------

  do K = Ksat_ini, NSOIL
    do KSUB = 1, NSUB(K)
      ISUB = sum(NSUB(1:K))-NSUB(K)+KSUB
      if (SOLVE_ORDER(ISUB) .eq. 3) then
        dsmc_tmp = smc_sub(ISUB) - smc(K)  !now dsmc_tmp is <0
        call calculate_wt_ratio(dsmc_tr_acc, dsmc_acc, PROPOR_TR)
        dsmc_acc = dsmc_acc + dsmc_tmp*SLDPTH(K)/NSUB(K)
        dsmc_tr_acc = dsmc_tr_acc + dsmc_tmp*SLDPTH(K)/NSUB(K)*PROPOR_TR
        smc_tr(ISUB) = smc_tr(ISUB) - dsmc_tmp*PROPOR_TR
        smc_sub(ISUB) = smc(K)
      end if 
    end do
  end do
  if (abs(dsmc_acc) .gt. 1.0E-7) then
    print *, "ERROR in adjust_smc_tr2"
    print *, "dsmc_acc, dsmc_tr_acc,I,J:"
    print *, dsmc_acc, dsmc_tr_acc, I,J
  end if

end subroutine adjust_smc_tr
!------------------------------------------------------

!------------------------------------------------------
! Huancui Hu: subroutine adjust_smc_tr2
subroutine adjust_smc_tr2(I,J,NSOIL, NTRACER, SLDPTH, NSUB,&
                          KK, smc, smc_sub, smc_tr,&
                          SUBFLO_SUB, SUBFLO_TR)
  implicit none
  integer, intent(in)                       :: I,J
  integer, intent(in)                       :: NSOIL
  integer, intent(in)                       :: NTRACER
  real,    intent(in),   dimension(NSOIL)   :: SLDPTH
  integer, intent(in),   dimension(NSOIL)   :: NSUB
  integer, intent(in)                       :: KK     ! the level adjust in bulk soil water
  real,    intent(in),   dimension(NSOIL)   :: smc
  real,    intent(inout),dimension(NTRACER) :: smc_sub
  real,    intent(inout),dimension(NTRACER) :: smc_tr
  real,    intent(inout)                    :: subflo_sub
  real,    intent(inout)                    :: subflo_tr
  !-----local vars-----
  integer        :: K, KSUB, ISUB
  real           :: PROPOR_TR
  real           :: dsmc_tmp
  real           :: dsmc_acc
  real           :: dsmc_tr_acc
  integer,dimension(NTRACER) :: SOLVE_ORDER

if (I .eq. 173 .and. J .eq. 91) then
!  print *, "FLAG IN adjust_smc_tr2001, KK, SMC, SMC_SUB, SUBFLO_SUB:",KK, SMC, SMC_SUB, SUBFLO_SUB
end if  
 
  dsmc_tmp = 0.0
  dsmc_acc = 0.0
  if (SUBFLO_SUB .gt. 1.0E-12) then   !  SUBFLO_SUB >0 while SUBLO<0?? correct, SUBFLO_SUB is always >0 because it only has the gain component
      call calculate_wt_ratio(subflo_tr, subflo_sub, PROPOR_TR)
      do K = NSOIL, 1, -1     ! not an error, do need to adjust each layer below Ksat_ini for tracer
        do KSUB = NSUB(K), 1, -1
          ISUB = sum(NSUB(1:K)-NSUB(K)+KSUB)
          dsmc_tmp = smc(K)-smc_sub(ISUB)    !theoretically dsmc_tmp can only be +
          if (dsmc_tmp .gt. subflo_sub/SLDPTH(K)) then
              call calculate_wt_ratio(subflo_tr, subflo_sub, PROPOR_TR)
              smc_tr(ISUB)  = smc_tr(ISUB) + subflo_sub/SLDPTH(K)*PROPOR_TR
              smc_sub(ISUB) = smc_sub(ISUB)+ subflo_sub/SLDPTH(K)
              subflo_sub  = 0.0
              subflo_tr   = 0.0
          else
              if (dsmc_tmp .gt. 1.0E-10) then   !dsmc_tmp >0
                call calculate_wt_ratio(subflo_tr, subflo_sub, PROPOR_TR)
                smc_tr(ISUB)  = smc_tr(ISUB) + dsmc_tmp*PROPOR_TR
                smc_sub(ISUB) = smc(K)
                subflo_sub  = subflo_sub - dsmc_tmp*SLDPTH(K)/NSUB(K)
                subflo_tr   = subflo_tr  - dsmc_tmp*SLDPTH(K)/NSUB(K)*PROPOR_TR
              else
                if (dsmc_tmp .lt. -1.0E-7) then
                  print *, "ERROR in adjust_smc_tr201,I,J, dsmc_tmp,K,smc(K),smc_sub(ISUB):",I,J,dsmc_tmp,K,smc(K),smc_sub(ISUB)
                  call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
                  smc_tr(ISUB)  = smc_tr(ISUB) + dsmc_tmp*PROPOR_TR
                  smc_sub(ISUB) = smc(K)
                  subflo_sub  = subflo_sub - dsmc_tmp*SLDPTH(K)/NSUB(K)
                  subflo_tr   = subflo_tr  - dsmc_tmp*SLDPTH(K)/NSUB(K)*PROPOR_TR
                else
                  smc_sub(ISUB) = smc(K)
                  smc_tr(ISUB)  = max(0.0,min(smc_tr(ISUB),smc_sub(ISUB)))
                end if
              end if
          end if
        end do
      end do
  end if
  ! now subflow_sub should be close to 0
  if (abs(subflo_sub) .lt. 1.0E-6) then
    subflo_sub = 0.0
    subflo_tr  = 0.0
  else
    print *, "ERROR in adjust_smc_tr202, subflo_sub,subflo_tr, I,J,K, ISUB, smc, smc_sub, smc_tr:",subflo_sub,subflo_tr, I,J,K, ISUB, smc, smc_sub, smc_tr
  end if 
              

end subroutine
!------------------------------------------------------

!------------------------------------------------------
! Huancui Hu: calculate_qsubout_tr
subroutine calculate_qsubout_tr(ii, jj, satk,nsoil,ntracer,pmopt,nsub,sldepth,&
                                qqsub_lose,grid_area,subdt,  &
                                qqsub_tr_lose,  &
                                qqsub_sub_lose, &
                                smc_tr, smc_sub)
  implicit none
  integer, intent(in)                            :: ii, jj, satk
  integer, intent(in)                            :: nsoil
  integer, intent(in)                            :: ntracer
  integer, intent(in)                            :: pmopt
  integer, intent(in),   dimension(nsoil)        :: nsub
  real,    intent(in),   dimension(nsoil)        :: sldepth
  real,    intent(in)                            :: qqsub_lose     !unit:m3/s, has to be <0
  real,    intent(in)                            :: grid_area     !unit:m2
  real,    intent(in)                            :: subdt    
  real,    intent(out),  dimension(ntracer)      :: qqsub_tr_lose
  real,    intent(out),  dimension(ntracer)      :: qqsub_sub_lose
  real,    intent(inout),dimension(ntracer)      :: smc_tr
  real,    intent(inout),dimension(ntracer)      :: smc_sub
  !-----local vars-----
  integer     :: KSUB, ISUB, K
  real        :: PROPOR_TR
  real        :: qqsub_tmp
  real        :: qqsub_sub_tmp
  !--------------------

  qqsub_tr_lose  = 0.0
  qqsub_sub_lose = 0.0
  qqsub_sub_tmp  = 0.0
  do K = satk, nsoil
    qqsub_sub_tmp = qqsub_lose*sldepth(K)/sum(sldepth(satk:nsoil))    !qqsub_sub_tmp<0
    do KSUB = 1, NSUB(K)
      ISUB = sum(NSUB(1:K)) - NSUB(K) + KSUB
      if (pmopt .eq. 2  .or. pmopt .eq. 12) then
        call calculate_wt_ratio_pm(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
      else
        call calculate_wt_ratio(smc_tr(ISUB), smc_sub(ISUB), PROPOR_TR)
      end if
      qqsub_tr_lose(ISUB)  = qqsub_tr_lose(ISUB)  + qqsub_sub_tmp/NSUB(K)*PROPOR_TR
      smc_tr(ISUB)         = smc_tr(ISUB)  + (qqsub_sub_tmp/grid_area*subdt)/sldepth(K)*PROPOR_TR
      qqsub_sub_lose(ISUB) = qqsub_sub_lose(ISUB) + qqsub_sub_tmp/NSUB(K)
      smc_sub(ISUB)        = smc_sub(ISUB) + (qqsub_sub_tmp/grid_area*subdt)/sldepth(K)
      if (smc_tr(ISUB) .lt. -1.0E-8) then
        print *, "ERROR in calculate_qsubout_tr00, smc_tr<0,ii,jj,pmopt,PROPOR_TR,smc_tr,smc_sub:",ii,jj,pmopt,PROPOR_TR,smc_tr,smc_sub
      end if
    end do
  end do
  if (abs(sum(qqsub_sub_lose(1:NTRACER)) - qqsub_lose) .gt. 1.5E-7 ) then
    if (abs(qqsub_lose) .lt. 1.0) then
      print *, "ERROR in calculate_qsubout_tr!"
      print *, "satk, nsoil, qqsub_lose, qqsub_sub_lose:"
      print *, ii,jj, satk, nsoil, qqsub_lose, qqsub_sub_lose
    else
      if (abs(sum(qqsub_sub_lose(1:NTRACER)) - qqsub_lose) .gt. 1.5E-6) then
        print *, "ERROR in calculate_qsubout_tr!"
        print *, "satk, nsoil, qqsub_lose, qqsub_sub_lose:"
        print *, ii,jj,satk, nsoil, qqsub_lose, qqsub_sub_lose
      else
!        qqsub_sub_lose = qqsub_lose
!        qqsub_tr_lose  = max(qqsub_tr_lose, qqsub_lose)      ! is negative
      end if
    end if
  else
!    qqsub_sub_lose = qqsub_lose
!    qqsub_tr_lose  = max(qqsub_tr_lose, qqsub_lose)      ! is negative
  end if
end subroutine calculate_qsubout_tr
!------------------------------------------------------

end module module_subsurfaceRouting 
