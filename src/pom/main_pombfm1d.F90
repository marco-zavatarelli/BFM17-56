#include"cppdefs.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! This is the "main" program of the coupled numerical model originating by 
! the direct on-line coupling of the 1D Version of the Princeton Ocean model
! "POM" and the Biological Flux Model "BFM".
! The whole modelling system is available for download from the BFM web site:
!
! bfm-community.eu
!
! This release has been finalised by Marco Zavatarelli, Giulia Mussap and 
! Nadia Pinardi. However, previous very significant past contributions were 
! provided  also by Momme Butenschoen and Marcello Vichi. 
! 
!                                              Marco Zavatarelli@unibo.it
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      PROGRAM MAIN
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use global_mem,ONLY: RLEN,ZERO,PI,ONE
!
      use constants, ONLY: SEC_PER_DAY
!
      use POM, ONLY: IDIAGN,                                              &    
                     KL1, KL2,                                            &
                     IHOTST,                                              &
                     DTI,                                                 &
                     IDAYS,                                               &
                     IEND,                                                &
                     INTT,                                                &
                     TIME,                                                &
                     TIME0,                                               &
                     H,                                                   &
                     ALAT,                                                &
                     COR,                                                 &
                     UMOL,                                                &
                     UMOLT,UMOLS,UMOLBFM,                                 &
                     SMOTH,                                               &
                     RCP,                                                 &
                     DAYI,                                                &
                     KB,                                                  &
                     NBCT,NBCS,NBCBFM,                                     &
                     NTP,                                                 &
                     TRT, SRT,                                            &
                     upperH,                                              &
                     SSRT,                                                &
                     Z,ZZ,DZ,DZZ,DZR,                                     &
                     TF,T,TB,                                             &
                     SF,S,SB,                                             &
                     RHO,                                                 & 
                     UF,U,UB,VF,V,VB,                                     &
                     Q2F,Q2,Q2B,                                          &
                     L,                                                   &
                     Q2LF,Q2L,Q2LB,                                       &
                     KM,KH,KQ,                                            & 
                     GM,GH,SM,SH,KN,SPROD,BPROD, A, C, VH, VHP,PROD,DTEF, &
                     D,DT,                                                &
                     WUSURF,WVSURF,                                       &
                     WUBOT,WVBOT,                                         &
                     SWRAD,                                               &
                     WTSURF,                                              &
                     TSURF, SSURF,                                        &
                     WSSURF,                                              &
                     WTADV, WSADV,                                        &
                     NRT_o2o,NRT_n1p,NRT_n3n,NRT_n4n
!
      use Service,ONLY: ilong, savef
!
      use Forcing,ONLY: Forcing_manager
!      use netcdf_pom, only: init_netcdf_pom,init_save_pom
!
!------------------------------------------------------------------------!
!
!BOC
!
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Implicit typing is never allowed
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      IMPLICIT NONE
!
!     -----NAMELIST READING UNIT------
!
      integer(ilong),parameter            :: namlst=10
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Scalars
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!     -----LOOP COUNTER-----
!
      integer(ilong)               :: K
!
!     -----TWICE THE TIME STEP-----
!
      REAL(RLEN)                   :: DT2
!
!     -----EARTH ROTATION ANGULAR VELOCITY-----
!
      REAL(RLEN),parameter         ::OMEGA=7.29E-5
!
!     -----INTERPOLATED CLIMATOLOGICAL T AND S PROFILES-----
!
      REAL(RLEN),dimension(KB)   :: TSTAR,SSTAR
!
!     -----FORCING COUNTERS-----
!
      INTEGER(ilong)               :: ICOUNTF,IDOUNTF,IFCHGE,IFDCHGE,&
                                      IFDINT,IFINT,IFNCHGE, IFNINT
!
!     -----MONTHLY CLIMATOLOGICAL T AND S PROFILES-----
!
      REAL(RLEN),SAVE              :: SCLIM1(KB),SCLIM2(KB),TCLIM1(KB),&
                                      TCLIM2(KB)
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Arrays
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      REAL(RLEN),SAVE              :: WSU1,WSV1,SSS1,SSS2,SLUX1,SLUX2,SWRAD1,WTSURF1,QCORR1, &
                                      QCORR2,NO3_s1,NO3_s2,NH4_s1,NH4_s2,PO4_s1,PO4_s2,SIO4_s1,SIO4_s2, &
                                      O2_b1,O2_b2,NO3_b1,NO3_b2,PO4_b1,PO4_b2,WCLIM1,WCLIM2
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  External Subroutines
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      EXTERNAL DENS,PROFQ,PROFT,PROFU,PROFV, CALCDEPTH
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     !  Intrinsic Functions
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      INTRINSIC FLOAT,MOD,NINT,SIN,IFIX
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  NAMELIST READING
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      NAMELIST /Params_POMBFM/ H,DTI,ALAT,IDIAGN,IDAYS,SMOTH,&
                               ihotst,UMOL,KL1,KL2,savef,NRT_O2o,NRT_N1p,NRT_N3n,NRT_N4n, &
                               NBCT,NBCS,NBCBFM,UMOL,UMOLT,UMOLS,UMOLBFM,NTP,TRT,SRT, &
                               UPPERH,SSRT
!
      OPEN(namlst,file='params_POMBFM.nml',status='old',action='read')
      READ(namlst,nml=Params_POMBFM)
      CLOSE(namlst)
!
!
!     -----GENERAL INITIALISATION-----
!
          Z(:)           = ZERO
          ZZ(:)          = ZERO
          DZ(:)          = ZERO
          DZZ(:)         = ZERO
          TF(:)          = ZERO
          T(:)           = ZERO
          TB(:)          = ZERO
          SF(:)          = ZERO
          S(:)           = ZERO
          SB(:)          = ZERO
          RHO(:)         = ZERO
          UF(:)          = ZERO
          U(:)           = ZERO
          UB(:)          = ZERO
          VF(:)          = ZERO
          V(:)           = ZERO
          VB(:)          = ZERO 
          Q2F(:)         = 1.E-7_RLEN
          Q2(:)          = Q2F(:)
          Q2B(:)         = Q2(:)
          Q2LF(:)        = 1.E-7_RLEN
          Q2L(:)         = Q2LF(:)
          Q2LB(:)        = Q2L(:)
          L(:)           = 1.0_RLEN
          L(1)           = ZERO
          L(KB)          = ZERO
          KM(:)          = ZERO
          KH(:)          = ZERO
          KQ(:)          = ZERO 
          GM(:)          = ZERO
          GH(:)          = ZERO
          SM(:)          = ZERO
          SH(:)          = ZERO
          KN(:)          = ZERO
          SPROD(:)       = ZERO
          BPROD(:)       = ZERO
          PROD(:)        = ZERO    
          VH(:)          = ZERO
          VHP(:)         = ZERO
          DTEF(:)        = ZERO
          D(:)           = ZERO
          DT(:)          = ZERO
          A(:)           = ZERO 
          C(:)           = ZERO
          WTADV(:)       = ZERO
          WSADV(:)       = ZERO
          WUSURF         = ZERO
          WVSURF         = ZERO
          WUBOT          = ZERO
          WVBOT          = ZERO
          WTSURF         = ZERO
          SWRAD          = ZERO 
          WSSURF         = ZERO
          TSURF          = ZERO
          SSURF          = ZERO
!
!     -----DEFINE VERTICAL COORDINATE-----
!
      call calcdepth(z,zz,dz,dzz,kb,KL1,KL2)
!     
      dz(kb)=1.e-6_RLEN
      dzr(:)=ONE/dz(:)
      dz(kb)=ZERO
!
!     -----CORIOLIS PARAMETER-----
!
!      COR = 1.e-4.0_RLEN*SIN(ALAT*2.0_RLEN*PI/360.0_RLEN)
      COR = 2.0_RLEN*OMEGA*SIN(ALAT*2.0_RLEN*PI/360.0_RLEN)
!
!     -----TWICE THE TIME STEP-----
!
      DT2  = 2.0_RLEN*DTI
!
!     -----ITERATIONS NEEDED TO CARRY OUT AN"IDAYS" SIMULATION-----
!
      iend = idays*IFIX(SEC_PER_DAY)/IFIX(dti)
!
!     -----READ  T&S INITIAL CONDITIONS (IHOTST=0) OR RESTART FILE (IHOTST=1)-----
!
      select case (ihotst)
!
             case (0)
!
                  time0=ZERO
!
                  call get_TS_IC
!       
!                 -----DEFINE INITIAL DENSITY FIELD----
!    
                  call DENS(T,S,ZZ,H,RHO,KB)
!
             case (1)
!
!                 -----READ RESTART-----
                  call get_rst
! 
       end select 
!
#ifndef POM_only
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Initialization of BFM 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       call pom_ini_bfm_1d
#endif

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Begin the time march 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      WRITE(6,*) 'Icount before time march loop = ',ICOUNTF

!      -----BFM SET-UP, INITIALISATION AND RESTART READING (IF REQUIRED)-----
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  Begin the time march 
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      DO intt = 1, IEND           
!
!         -----COMPUTE TIME IN DAYS-----
!
          time = time0+(DTI*float(intt)*dayi)       
!
!         -----TURBULENCE CLOSURE-----
        
          Q2F(:)  = Q2B(:)
          Q2LF(:) = Q2LB(:)
!
          call PROFQ(DT2)
!
!         -----DEFINE ALL FORCINGS-----
!  
          call forcing_manager

!         -----T&S COMPUTATION-----
!
     select case (IDIAGN)
!           
          case (INT(ZERO))
!
!            *************************************
!            *************************************
!            **                                 **
!            ** PROGNOSTIC MODE:                **
!            ** T&S FULLY COMPUTED BY MODEL     **
!            **                                 **
!            *************************************
!            *************************************
!
!            -----COMPUTE LATERAL ADVECTION TERM FOR T&S-----
!
             wtadv(:) = ZERO
             wsadv(:) = ZERO
!
             IF(TRT.ne.ZERO) THEN
!
                do K = 1, KB
 
                   if((-ZZ(K)*H).ge.upperH)   &
                   wtadv(K) = (tstar(k)-T(k))/(TRT*SEC_PER_DAY)
!
                enddo
!
             ENDIF
!
             IF(SRT.ne.ZERO) THEN
!
                do K = 1, KB
!
                   if((-ZZ(K)*H).ge.upperH)   &
                   wsadv(K) = (sstar(k)-S(k))/(SRT*SEC_PER_DAY)
!
                enddo
!
             ENDIF
!
!            -----COMPUTE SURFACE SALINITY FLUX-----
!
             wssurf = -(ssurf-S(1))*SSRT/SEC_PER_DAY
!
!            -----COMPUTE TEMPERATURE-----
!
             TF(:) = TB(:) + (WTADV(:) * DT2)
!
             CALL PROFTS(TF,WTSURF,SWRAD,TSURF,NBCT,DT2,NTP,UMOLT)
!
!            -----COMPUTE SALINITY-----
!
             SF(:) = SB(:) + (WSADV(:) * DT2)
!
             CALL PROFTS(SF,WSSURF,ZERO,SSURF,NBCS,DT2,NTP,UMOLS)
!!
!            -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!            !        MIXING THE TIMESTEP (ASSELIN)                            !
!            -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
             T(:)  = T(:) + 0.5_RLEN * SMOTH * (TF(:) + TB(:) - 2.0_RLEN * T(:))
             S(:)  = S(:) + 0.5_RLEN * SMOTH * (SF(:) + SB(:) - 2.0_RLEN * S(:))

     end select
!
!     -----COMPUTE VELOCITY-----
!
     UF(:) = UB(:) + DT2 * COR * V(:)
     call PROFU(DT2)
     VF(:) = VB(:) - DT2 * COR * U(:)
     call PROFV(DT2)

!
!     -----MIX TIME STEP (ASSELIN FILTER)-----
!
      Q2(:)   = Q2(:)  + 0.5_RLEN * SMOTH * (Q2F(:)  + Q2B(:)  - 2.0_RLEN * Q2(:))
      Q2L(:)  = Q2L(:) + 0.5_RLEN * SMOTH * (Q2LF(:) + Q2LB(:) - 2.0_RLEN * Q2L(:))
!
      U(:) = U(:) + 0.5_RLEN * SMOTH * (UF(:) + UB(:) - 2.0_RLEN * U(:))
      V(:) = V(:) + 0.5_RLEN * SMOTH * (VF(:) + VB(:) - 2.0_RLEN * V(:))
!
!     -----RESTORE TIME SEQUENCE-----
!
      Q2B(:)  = Q2(:)
      Q2(:)   = Q2F(:)
      Q2LB(:) = Q2L(:)
      Q2L(:)  = Q2LF(:)
!
      UB(:) = U(:)
      U(:)  = UF(:)
      VB(:) = V(:)
      V(:)  = VF(:)
!
      TB(:) = T(:)
      T(:)  = TF(:)
      SB(:) = S(:)
      S(:)  = SF(:)
!
!     -----UPDATE DENSITY-----
!
      call DENS(T,S,ZZ,H,RHO,KB)
!
!     -----PHYSICS HAS BEEN COMPUTED.....GO BFM!!!!-----
!
#ifndef POM_only
!
      call pom_bfm_1d
!
#endif
!
       END DO
!
!----- WRITING OF POM RESTART-----
!
      WRITE (71) time,                      &
                 u,ub,v,vb,                 &
                 t,tb,s,sb,                 &
                 q2,q2b,q2l,q2lb,           &
                 kh,km,kq,                  &
                 l,                         &
                 wubot,wvbot,               &
                 RHO
!
!-----BFM RESTART-----
!
#ifndef POM_only
!
!  -----WRITE RESTART -----
!
      call restart_BFM_inPOM

!
#endif
!
!
       PRINT *,'Main done'
!
     end program MAIN 
