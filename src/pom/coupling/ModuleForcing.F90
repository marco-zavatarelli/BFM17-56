!
!!!!!!
!WARNING THIS IS A TEST VERSION
!(MONTHLY FREQUENCY)!!!!!!!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!ONE-DIMENSIONAL BFM-POM SYSTEM
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !INTERFACE
  MODULE  Forcing
!
!
!................................................Marco.Zavatarelli@unibo.it
!................................................G.Mussap@sincem.unibo.it
!                                                    2014
!
! -----DEFINITION AND ALLOCATION OF PARAMETERS, SCALAR AND ARRAYS USED---
! -----TO DEFINE THE PHYSICAL AND  BIOGEOCHEMICAL FORCING-----
!
! !USES:
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
  use POM,ONLY: KB,ilong
  use global_mem, ONLY:RLEN
!
!-------------------------------------------------------------------------!
!BOC
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Implicit typing is never allowed
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
  IMPLICIT NONE
!
!     -----MONTHLY SHORTWAVE RADIATION-----
!
!     ***************************************
!     ***************************************
!     ** N.B.!                             **
!     ** ALWAYS NEEDED: WHEN THE MODEL IS  **
!     ** RUN IN DIAGNOSTIC MODE PROVIDES   **
!     ** ONLY PAR TO BFM. IN PROGNOSTIC    **
!     ** CONTRIBUTES TO THE DEFINITION OF  **
!     ** THE TEMPERATURE SURFACE BOUNDARY  **
!     ** CONDITION.                        **
!     ***************************************
!     ***************************************
!
      REAL(RLEN),SAVE                     :: SWRAD1,SWRAD2
!      REAL(RLEN),SAVE                     :: SLUX1,SLUX2
!
!     ***************************************
!     ***************************************
!     ** N.B.!                             **
!     ** THE FOLLOWING SCALARS ARE USED    **
!     ** ONLY WHEN THE MODEL IS RUN IN     **
!     ** PROGNOSTIC MODE.                  **
!     ***************************************
!     ***************************************
!
!     -----MONTHLY LOSS TERM OF THE SURFACE HEAT FLUX-----

       REAL (RLEN),SAVE                     :: WTSURF1,WTSURF2
!
!     ****************************************
!     ****************************************
!     ** N.B.                               **
!     ** THE FOLLOWING ARE ALWAYS USED      **
!     **                                    **
!     ****************************************
!     ****************************************
!
!     -----PRESCRIBED T&S PROFILES-----
!
      real(RLEN),public,dimension(KB),SAVE   :: TSTAR,SSTAR
!
!     -----MONTHLY WIND STRESS-----
!
      REAL (RLEN),SAVE                       :: WSU1,WSU2,WSV1,WSV2
!
!     -----MONTHLY SURFACE SALINITY

      REAL (RLEN),SAVE                       :: SSS1,SSS2
!
!     -----MONTHLY BOTTOM OXYGEN
!
      REAL (RLEN),SAVE                       :: O2_b1,O2_b2
!
!     -----MONTHLY SURFACE AND BOTTOM NITRATE
!
      REAL (RLEN),SAVE                       :: NO3_s1,NO3_s2
      REAL (RLEN),SAVE                       :: NO3_b1,NO3_b2
!
!     -----MONTHLY SURFACE AND BOTTOM PHOSPHATE-----
!
      REAL (RLEN),SAVE                       :: PO4_s1,PO4_s2
      REAL (RLEN),SAVE                       :: PO4_b1,PO4_b2
!
!     -----MONTHLY SURFACE AMMONIA-----
!
      REAL (RLEN),SAVE                       :: NH4_s1,NH4_s2
!
!     -----MONTHLY BOTTOM PON GRADIENT-----
!
      REAL (RLEN),SAVE                       :: PON_b1,PON_b2
!
!     -----MONTHLY SURFACE SILICATE-----
!
      REAL (RLEN),SAVE                       :: SIO4_s1,SIO4_s2
!
!     -----MONTHLY PROFILES OF INORGANIC SUSPENDED MATTER-----
!
      real(RLEN),public,dimension(KB-1),SAVE :: ISM1,ISM2
!
!     -----MONTHLY PROFILES OF T & S-----
!
      real(RLEN),public,dimension(KB),SAVE :: TCLIM1,TCLIM2
      real(RLEN),public,dimension(KB),SAVE :: SCLIM1,SCLIM2
      real(RLEN),public,dimension(KB),SAVE :: WCLIM1,WCLIM2
      real(RLEN),public,dimension(KB),SAVE :: WEDDY1,WEDDY2,WEDDY3,WEDDY4
!
       REAL (RLEN),SAVE                    :: SLUX1,QCORR1,QCORR2
!
!     -----INTERPOLATORS AND COUNTERS-----
!
      INTEGER(ilong),SAVE                   ::  ICOUNTF,IDOUNTF,          &
                                                IFCHGE,IFDCHGE,           &
                                                IFDINT,IFINT
!
      REAL(RLEN),SAVE                       ::  RATIOF,RATIOD
!
contains
!
!
      SUBROUTINE FORCING_MANAGER
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use global_mem,ONLY: RLEN,ZERO,PI,ONE,NML_OPEN,NML_READ,error_msg_prn
!
      use constants, ONLY: SEC_PER_DAY
!
      use POM, ONLY: IDIAGN,                                              &
                     DTI,                                                 &
                     INTT,                                                &
                     RCP,                                                 &
                     KB,                                                  &
                     TF,                                                  &
                     SF ,                                                 &
                     WUSURF,WVSURF,                                       &
                     SWRAD,                                               &
                     WTSURF,WSSURF,                                       &
                     TSURF,SSURF,                                         &
                     TB,SB,                                               &
                     ilong
!
      use Service,ONLY: ISM, savef, &
                        PO4SURF,NO3SURF,NH4SURF,SIO4SURF,                 &
                        PO4BOTT,NO3BOTT,O2BOTT,PONBOTTgrad,WGEN,WEDDY,    &
                        wind_input,radiance_input,                        &
                        ism_input,Sprofile_input,W_input,                 &
                        Tprofile_input,surfNut_input,bottNut_input,       &
                        Sal_input,Temp_input,heat_input,SurfaceS_input
!------------------------------------------------------------------------!
!
!BOC
!
!     -----LOOP COUNTER-----
!
      integer(ilong)               :: K
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Arrays
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      INTEGER                      :: RLENGTH
!
!     -------INITIALISATION AND FIRST FORCING READING-----
!
      if(intt==INT(ONE)) then
!
        CALL opendat(WSU1,WSV1,SSS1,SLUX1,ISM1,SCLIM1,TCLIM1,WCLIM1,&
             WEDDY1,WEDDY2,SB,TB,SWRAD1,WTSURF1,QCORR1,NO3_s1,NH4_s1,&
             PO4_s1,SIO4_s1,O2_b1,NO3_b1,PO4_b1,PON_b1)
!
!         -----DAY COUNTER----
!
          IDOUNTF=1
!
!         -----MONTH COUNTER-----
!
          ICOUNTF =  1
!
!         -----TIME STEPS TO COVER ONE DAY-----
!
          IFDCHGE=INT(SEC_PER_DAY)/INT(DTI)
!
!         -----TIME STEPS TO COVER ONE MONTH----
!
          IFCHGE  = 30_ilong*IFDCHGE
!
!         -----DAY INTERPOLATOR---
!
!         ******************************************************************
!         ******************************************************************
!         **                                                              **
!         ** THE DAILY   CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE    **
!         ** CENTERED AT h 00.00  OF EACH CLIMATOLOGICAL DAY. THEREFORE   **
!         ** THE MONTH INTERPOLATOR (IFDINT) IS INITIALISED AT THE VALUE  **
!         ** CORRESPONDING TO MIDNIGHT  MINUS 1 TIMESTEP.                 **
!         **                                                              **
!         ******************************************************************
!         ******************************************************************
!
          IFDINT   = -INT(ONE)
!
!         -----MONTH INTERPOLATOR---
!
!         ******************************************************************
!         ******************************************************************
!         **                                                              **
!         ** THE MONTHLY CLIMATOLOGICAL FORCING DATA ARE ASSUMED TO BE    **
!         ** CENTERED AT DAY 15 OF EACH CLIMATOLOGICAL MONTH. THEREFORE   **
!         ** THE MONTH INTERPOLATOR (IFINT) IS INITIALISED AT THE VALUE   **
!         ** (IFCHGE/2)-1 CORRESPONDING TO DAY 15 MINUS 1 TIMESTEP.       **
!         **                                                              **
!         ******************************************************************
!         ******************************************************************
!
          IFINT   = (IFCHGE/2_ilong)-INT(ONE)
!
!         *************************************************
!         *************************************************
!         **                                             **
!         **  INITIAL READING OF THE MONTHLY FORCING     **
!         **                                             **
!         *************************************************
!         *************************************************
!
!         -----WIND STRESS-----
!
          READ (11,REC=ICOUNTF)   WSU1,WSV1
          READ (11,REC=ICOUNTF+1) WSU2,WSV2
!
!         -----SURFACE HEAT FLUX-----
!
!         **************************************************
!         **************************************************
!         **                                              **
!         ** N.B.: IF THE MODEL IS RUN IN DIAGNOSTIC MODE **
!         ** ONLY THE SOLAR RADIATION  (SWRAD#) IS USED   **
!         **                                              **
!         **************************************************
!         **************************************************
!
          select case (IDIAGN)
!
                 case(0)
!
                    READ(21,REC=ICOUNTF)   SWRAD1,WTSURF1
                    READ(21,REC=ICOUNTF+1) SWRAD2,WTSURF2
!
                 case(1)
!
!            ---------------DAILY----------------
!                    READ (17,REC=IDOUNTF)   SWRAD1
!                    READ (17,REC=IDOUNTF+1) SWRAD2
!------------------------MONTHLY--------------------------
                     READ (21,REC=ICOUNTF)   SWRAD1,WTSURF1,QCORR1
                     READ (21,REC=ICOUNTF+1) SWRAD2,WTSURF2,QCORR2
!                     READ (17,REC=ICOUNTF)   SLUX1
!                     READ (17,REC=ICOUNTF+1) SLUX2

           end select
!
!         ***********************************************************
!         ***********************************************************
!         **                                                       **
!         ** ELIMINATE SOLAR RADIATION from TOTAL HEAT FLUX        **
!         ** N.B.: THE FOLLOWING TWO LINES SHOULD BE ACTIVE IF     **
!         **       (AND ONLY IF!) IN WTSURF# THERE IS THE TOTAL    **
!         **       HEAT FLUX AND NOT ONLY THE LOSS TERMS           **
!         **                                                       **
!         ***********************************************************
!         ***********************************************************
!
!         WTSURF1 = WTSURF1-SWRAD1
!         WTSURF2 = WTSURF2-SWRAD2
!
!
!         -----CLIMATOLOGICAL T&S PROFILES-------
!
!         ***********************************************************
!         ***********************************************************
!         **                                                       **
!         ** N.B.:IF THE MODEL IS RUN IN DIAGNOSTIC MODE  THE      **
!         **      WHOLE PROFILES ARE USED. IN PROGNOSTIC MODE      **
!         **      ONLY THE SURFACE T&S VALUES ARE RETAINED AND     **
!         **      STORED IN TSURF AND SURF AS OPTION FOR THE       **
!         **      T&S SURFACE BOUNDARY CONDITION                   **
!         **                                                       **
!         ***********************************************************
!         ***********************************************************
!
          DO K = 1,KB
             READ (20,REC=(ICOUNTF-1)*KB+K)  SCLIM1(K)
             READ (15,REC=(ICOUNTF-1)*KB+K)  TCLIM1(K)
             READ (215,REC=(ICOUNTF-1)*KB+K) WCLIM1(K)
             READ (216,REC=(ICOUNTF-1)*KB+K) WEDDY1(K)
             READ (217,REC=(ICOUNTF-1)*KB+K) WEDDY2(K)
             READ (20,REC=ICOUNTF*KB+K)      SCLIM2(K)
             READ (15,REC=ICOUNTF*KB+K)      TCLIM2(K)
             READ (215,REC=ICOUNTF*KB+K)     WCLIM2(K)
             READ (216,REC=ICOUNTF*KB+K)     WEDDY3(K)
             READ (217,REC=ICOUNTF*KB+K)     WEDDY4(K)

          END DO
!
!         -----SUSPENDED INORGANIC MATTER PROFILES-----
!
          DO K = 1,KB-1
             READ (19,REC=(ICOUNTF-1)*(KB-1)+K) ISM1(K)
             READ (19,REC=ICOUNTF*(KB-1)+K)     ISM2(K)
          END DO
!
!         -----SURFACE NUTRIENTS-----
!
          READ(18,REC=ICOUNTF)    NO3_s1,NH4_s1,PO4_s1,SIO4_s1
          READ(18,REC=ICOUNTF+1)  NO3_s2,NH4_s2,PO4_s2,SIO4_s2
!
!         -----BOTTOM NUTRIENTS-----
!
          READ(300,REC=ICOUNTF)   O2_b1,NO3_b1,PO4_b1,PON_b1
          READ(300,REC=ICOUNTF+1) O2_b2,NO3_b2,PO4_b2,PON_b2
!
!
!         -----WIND STRESS CONVERTED TO POM UNITS (N/m2-->m2/s2)-----
!
          WSU1 = -WSU1*1.e-3_RLEN
          WSU2 = -WSU2*1.e-3_RLEN
          WSV1 = -WSV1*1.e-3_RLEN
          WSV2 = -WSV2*1.e-3_RLEN
!
!         -----HEAT FLUX CONVERTED TO POM UNITS(W/m2-->deg.C*m/s)-----
!
          SWRAD1  = -SWRAD1/rcp
          SWRAD2  = -SWRAD2/rcp
          WTSURF1 = -WTSURF1/rcp
          WTSURF2 = -WTSURF2/rcp
!
!         -----VERTICAL VELOCITY CONVERTED TO POM UNITS (m/s-->m/s)-----
!
          WCLIM1  = WCLIM1
          WCLIM2  = WCLIM2
!
!         -----UPDATE THE DAY COUNTER-----
!
          IDOUNTF=IDOUNTF + INT(ONE)
!
!         -----UPDATE THE MONTH COUNTER-----
!
          ICOUNTF = ICOUNTF + INT(ONE)
!
      endif
!
!     -----UPDATE INTERPOLATION COUNTERS-----
!
      IFDINT = IFDINT + INT(ONE)
      RATIOD = FLOAT(IFDINT)/FLOAT(IFDCHGE)
!
      IFINT  = IFINT + INT(ONE)
      RATIOF = FLOAT(IFINT)/FLOAT(IFCHGE)
!
!     -----INTERPOLATE WIND STRESS-----
!
      WUSURF = WSU1 + RATIOF * (WSU2-WSU1)
      WVSURF = WSV1 + RATIOF * (WSV2-WSV1)
!
!
!     -----INTERPOLATE HEAT FLUX-----
!
      select case (IDIAGN)
!
             case (INT(ZERO))
!
                  WTSURF = WTSURF1 + RATIOF * (WTSURF2-WTSURF1)
                  SWRAD  = SWRAD1  + RATIOF * (SWRAD2-SWRAD1)
!
             case (INT(ONE))
!
!               ------------DAILY----------------
!                  SWRAD  = SWRAD1  + RATIOD * (SWRAD2-SWRAD1)
!               ------------MONTHLY----------------
                  SWRAD  = SWRAD1  + RATIOF * (SWRAD2-SWRAD1)
!                  SLUX   = SLUX1   + RATIOF * (SLUX2-SLUX1)
!
      end select
!
!     -----INTERPOLATE T&S PROFILES-----
!
      TSTAR(:) = TCLIM1(:) + RATIOF * (TCLIM2(:)-TCLIM1(:))
      SSTAR(:) = SCLIM1(:) + RATIOF * (SCLIM2(:)-SCLIM1(:))
      WGEN(:)  = WCLIM1(:) + RATIOF * (WCLIM2(:)-WCLIM1(:))
      if(RATIOF.le.0.5)then
         WEDDY(:) = WEDDY1(:)
      else
         WEDDY(:) = WEDDY2(:)
      endif
!
      select case (IDIAGN)
!
             case (INT(ZERO))
!
                  TSURF = TSTAR(1)
                  SSURF = SSTAR(1)
!
             case (INT(ONE))
!
                  TF(:) = TSTAR(:)
                  SF(:) = SSTAR(:)
!
      end select
!
!     -----INTERPOLATE SUSPENDED INORGANIC MATTER-----
!
      ISM(:) = ISM1(:) + RATIOF * (ISM2(:)-ISM1(:))
!
!     -----INTERPOLATE SURFACE NUTRIENTS-----
!
      NO3SURF  = NO3_s1  + RATIOF * (NO3_s2-NO3_s1)
      NH4SURF  = NH4_s1  + RATIOF * (NH4_s2-NH4_s1)
      PO4SURF  = PO4_s1  + RATIOF * (PO4_s2-PO4_s1)
      SIO4SURF = SIO4_s1 + RATIOF * (SIO4_s2-SIO4_s1)
!
!     -----INTERPOLATE BOTTOM NUTRIENTS-----
!
      O2BOTT      = O2_b1   + RATIOF * (O2_b2-O2_b1)
      NO3BOTT     = NO3_b1  + RATIOF * (NO3_b2-NO3_b1)
      PO4BOTT     = PO4_b1  + RATIOF * (PO4_b2-PO4_b1)
      PONBOTTgrad = PON_b1  + RATIOF * (PON_b2-PON_b1)
!
!      IF (IFDINT==IFDCHGE.AND.IDIAGN==INT(ONE)) THEN
!
!        ----A DAY HAS GONE...IT IS NECESSARY TO....-----
!
!        ----UPDATE DAY COUNTER....----
!
!         IDOUNTF = IDOUNTF+1
!
!        ----RESET INTERPOLATOR-----
!
!         IFDINT = INT(ZERO)
!
!        -----SHIFT THE DAILY DATA....-----
!
!         SWRAD1=SWRAD2
!
!        -----IF 360 DAYS HAVE GONE, RESTART THE READING SEQUENCE-----
!
!         if (IDOUNTF.GT.361) then
!            IDOUNTF=2_ilong
!            READ(17,REC=1) SWRAD1
!            SWRAD1=-SWRAD1/rcp
!         endif
!
!      ENDIF
!
!
      IF (IFINT==IFCHGE) THEN
!
!        -----A MONTH HAS GONE...IT IS NECESSARY TO....-----
!
!        -----....UPDATE MONTH COUNTER....-----
!
         ICOUNTF = ICOUNTF + 1
!
         PRINT *, 'ICOUNTF', ICOUNTF
!
!        -----....RESET INTERPOLATOR....-----
!
         IFINT   = INT(ZERO)
!
!        -----....SHIFT THE MONTHLY DATA....-----
!
         WSU1       = WSU2
         WSV1       = WSV2
!         if (IDIAGN==INT(ZERO)) then
            SWRAD1     = SWRAD2
            WTSURF1    = WTSURF2
!            SLUX1      = SLUX2
!         endif
         NO3_s1      = NO3_s2
         NH4_s1      = NH4_s2
         PO4_s1      = PO4_s2
         SIO4_s1     = SIO4_s2
         NO3_b1      = NO3_b2
         O2_b1       = O2_b2
         PO4_b1      = PO4_b2
         PON_b1      = PON_b2
         ISM1(:)     = ISM2(:)
         TCLIM1(:)   = TCLIM2(:)
         SCLIM1(:)   = SCLIM2(:)
         WCLIM1(:)   = WCLIM2(:)
         WEDDY1(:)   = WEDDY3(:)
         WEDDY2(:)   = WEDDY4(:)
!
         IF (ICOUNTF.GT.13) THEN
!
!            -----IF 12 MONTHS HAVE GONE, RESTART THE READING SEQUENCE-----
!
             ICOUNTF = 2
!
             READ (11,REC=1) WSU1,WSV1
             WSU1 = -WSU1*1.e-3_RLEN
             WSV1 = -WSV1*1.e-3_RLEN
!
!             if(IDIAGN==INT(ZERO)) then
                READ(21,REC=1) SWRAD1,WTSURF1
!                READ(17,REC=1) SLUX1
!               WTSURF1=WTSURF1-SWRAD1
                WTSURF1=-WTSURF1/rcp
                SWRAD1=-SWRAD1/rcp
!             endif
!
             READ(18,REC=1)  NO3_s1,NH4_s1,PO4_s1,SIO4_s1
             READ(300,REC=1) O2_b1,NO3_b1,PO4_b1,PON_b1
!
             DO K = 1,KB
                READ (20,REC=K)  SCLIM1(K)
                READ (15,REC=K)  TCLIM1(K)
                READ (215,REC=K) WCLIM1(K)
                READ (216,REC=K) WEDDY1(K)
                READ (217,REC=K) WEDDY2(K)
             END DO
!
             DO K = 1, KB-1
                READ (19,REC=K) ISM1(K)
             END DO
!
         END IF
!
!        -----READ FOLLOWING MONTH-----
!
         READ(18,REC=ICOUNTF)  NO3_s2,NH4_s2,PO4_s2,SIO4_s2
         READ(300,REC=ICOUNTF) O2_b2,NO3_b2,PO4_b2,PON_b2
!
         READ (11,REC=ICOUNTF) WSU2,WSV2
         WSU2 = -WSU2*1.e-3_RLEN
         WSV2 = -WSV2*1.e-3_RLEN
!
!         select case (IDIAGN)
!
!               case (INT(ZERO))
!
                     READ (21,REC=ICOUNTF) SWRAD2,WTSURF2
!                    WTSURF2=WTSURF2-SWRAD2
                     WTSURF2=-WTSURF2/rcp
                     SWRAD2 =-SWRAD2/rcp
!                     READ (17,REC=ICOUNTF) SLUX2
!
!        ------------- DAILY -------------------
!                case (INT(ONE))
!
!                     READ (17,REC=IDOUNTF) SWRAD2
!                     SWRAD2=-SWRAD2/rcp
!        ----------------------------------------
!
!         end select
!
         DO K = 1,KB
             READ (20,REC=(ICOUNTF-1)*KB+K)  SCLIM2(K)
             READ (15,REC=(ICOUNTF-1)*KB+K)  TCLIM2(K)
             READ (215,REC=(ICOUNTF-1)*KB+K) WCLIM2(K)
             READ (216,REC=(ICOUNTF-1)*KB+K) WEDDY3(K)
             READ (217,REC=(ICOUNTF-1)*KB+K) WEDDY4(K)
         END DO

         DO K = 1,KB-1
             READ (19,REC=(ICOUNTF-1)*(KB-1)+K) ISM2(K)
         END DO
!
      END IF
!
      return
!
      end
 end module Forcing
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
