!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!ONE-DIMENSIONAL BFM-POM SYSTEM  
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! !INTERFACE:
!
   MODULE POM
!
!................................................Marco.Zavatarelli@unibo.it
!................................................G.Mussap@sincem.unibo.it
!                                                    2014
! -----DEFINITION AND ALLOCATION OF PARAMETERS, SCALAR AND ARRAYS USED (MOSTLY)-----
! -----BY THE PHYSICAL COMPONENT OF THE SYSTEM-----
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
    use global_mem,ONLY: RLEN,ONE
    use constants, ONLY: SEC_PER_DAY
!
!-------------------------------------------------------------------------!
!
!BOC
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   IMPLICIT NONE
!
!     -----SET INTEGER PRECISION------
!
      integer,parameter                 :: ilong=selected_int_kind(12)
!
!     -----SWITCH FOR PROGNOSTIC/DIAGNOSTIC SIMULATION-----
!
!     *****************************************************
!     *****************************************************
!     **                                                 **
!     ** IDIAGN=0 PROGNOSTIC (T&S PROFILES COMPUTED)     **
!     ** IDIAGN=1 DIAGNOSTIC (T&S PROFILES PRESCRIBED)   **
!     **                                                 **
!     *****************************************************
!     *****************************************************
!
      integer(ilong)           ::IDIAGN
!
!     -----# OF SURFACE/BOTTOM LAYERS WITH LOG DISTRIBUTION-----
!
!     *****************************************************
!     *****************************************************
!     **                                                 **
!     ** KL1: SURFACE LAYERS                             **
!     ** KL2: BOTTOM LAYERS                              **
!     **                                                 **
!     *****************************************************
!     *****************************************************
!
      integer(ilong)          :: KL1, Kl2
!
!     -----SWITCH FOR COLD/HOT START-----
!
!     *****************************************************
!     *****************************************************
!     **                                                 **
!     ** IHOTST=0 "COLD" START FROM INITIAL CONDITION    **
!     ** IHOTST=1 "HOT" START FROM RESTART FILE          **
!     **                                                 **
!     ** SEE SUBROUTINE "CALCDEPTH"                      **
!     *****************************************************
!     *****************************************************
!
      integer(ilong)           ::IHOTST
!
!     -----MODEL TIME STEP-----
! 
      REAL(RLEN)               :: DTI
!
!     -----LENGTH OF THE RUN (DAYS)-----
!
      integer(ilong)           :: IDAYS
!
!     ----ITERATIONS NEEDED FOR AN "IDAYS" RUN-----
!  
      integer(ilong)           :: IEND
!
!     -----COUNTER FOR THE TIME MARCHING LOOP-----
!
      integer(ilong)           :: intt
!     
!     -----RUNNING TIME-----
!
      REAL(RLEN)               :: TIME
!
!     -----TIME AT RESTART-----
!
      REAL(RLEN)               :: TIME0
!
!     -----BOTTOM DEPTH-----
!
      REAL(RLEN)               :: H
!   
!     -----LATITUDE & LONGITUDE-----
!  
      REAL(RLEN)               :: ALAT, ALON
!    
!     -----CORIOLIS PARAMETER-----
!
      REAL(RLEN)               :: COR
!
!     -----BACKGROUND DIFFUSION FOR U, V, Q2, Q2L-----
!
      REAL(RLEN)               :: UMOL
!
!     -----BACKGROUND DIFFUSION FOR T,S & BFM TRACERS
!
      REAL(RLEN)               :: UMOLT,UMOLS,UMOLBFM
!
!     -----NUTRIENT RELAXATION TIME-----
!    
      REAL(RLEN)               :: NRT_o2o
      REAL(RLEN)               :: NRT_n1p
      REAL(RLEN)               :: NRT_n3n
      REAL(RLEN)               :: NRT_n4n
!
!     -----PARAMETER FOR THE HASSELIN FILTER----
!     -----(TIME STEP MIXING)-----
!
      REAL(RLEN)               :: SMOTH
!
!     -----SPECIFIC HEAT TIMES RHO0-----
!
      real(RLEN),parameter     :: RCP=4.187E6
!
!     -----1 DAY IN SECONDS (RECIPROCAL)-----
!
      real(RLEN), Parameter    :: DAYI=ONE/SEC_PER_DAY
!
!     -----VERTICAL LAYERS-----
!
      integer(ilong),parameter :: KB=151
!
!     -----FLAGS TO CHOOSE T,S AND BFM TRACERS SURFACE B.C. IN PROFTS-----
!
      integer(ilong)           :: NBCT, NBCS, NBCBFM
!
!     -----FLAG TO CHOOSE JERLOV WATER TYPE IN PROFTS-----
!
      integer(ilong)           :: NTP
!
!     -----T&S RELAXATION TIME (DAYS) FOR LATERAL ADVECTION-----
!
      integer(ilong)           :: TRT, SRT
!
!     -----DEPTH (m) AT WHICH LATERAL ADVECTION STARTS-----
!
      REAL(RLEN)               :: upperH 
!
!     -----RELAXATION TIME (DAYS) FOR SURFACE SALINITY FLUX-----
!
      REAL(RLEN)               :: SSRT
!
!     -----VERTICAL COORDINATE SYSTEM-----
!     
      real(RLEN),dimension(KB) :: Z,ZZ,DZ,DZZ,DZR     
!
!     -----TEMPERATURE----
!
      real(RLEN),dimension(KB) :: TF,T,TB
!
!     -----SALINITY-----
!
      real(RLEN),dimension(KB) :: SF,S,SB
!
!     -----DENSITY-----
!
      real(RLEN),dimension(KB) :: RHO
!
!     -----VELOCITY-----
!
      real(RLEN),dimension(KB) :: UF,U,UB,VF,V,VB
!
!     -----TURBULENT KINETIC ENERGY (T.K.E.X2)-----
!
      real(RLEN),dimension(KB) :: Q2F,Q2,Q2B
!
!     -----LENGTH SCALE-----
!
      real(RLEN),dimension(KB) :: L
!
!     -----(T.K.E.X2) TIMES LENGTH SCALE-----
!
      real(RLEN),dimension(KB) :: Q2LF,Q2L,Q2LB
!
!     -----VERTICAL DIFFUSION COEFFICIENTS-----
!
      real(RLEN),dimension(KB) ::KM,KH,KQ
!
!     -----SERVICE ARRAYS USED IN POM ROUTINES-----
!
      real(RLEN),public,dimension(KB) :: GM,GH,SM,SH,KN,SPROD,BPROD, A, C, VH, VHP, & 
                                         PROD,DTEF,D, DT
!
!     -----WIND STRESS----
!
      real(RLEN)               :: WUSURF,WVSURF
!
!     -----BOTTOM STRESS----
!
      real(RLEN)               :: WUBOT,WVBOT
!
!     *****************************************
!     *****************************************
!     **                                     **
!     ** N.B.                                **
!     **  WHEN THE MODEL IS RUN IN           **
!     **  DIAGNOSTIC MODE, THIS IS USED ONLY **
!     **  TO PROVIDE PAR TO THE              **
!     **  BIOGEOCHEMICAL COMPONENT           **
!     **                                     **
!     *****************************************
!     *****************************************
!
      real(RLEN)               :: SWRAD
!
!     *****************************************
!     *****************************************
!     **                                     **
!     ** N.B.                                **
!     ** THE FOLLOWING SCALARS ARE USED ONLY **
!     ** WHEN THE MODEL IS RUN IN PROGNOSTIC **
!     ** MODE                                **
!     **                                     **
!     *****************************************
!     *****************************************
!
!     -----LOSS TERM OF THE SURFACE HEAT FLUX-----
!
      real(RLEN)               :: WTSURF
! 
!     ------PRESCRIBED SURFACE T & S-----
!
!     *****************************************
!     *****************************************
!     **                                     **
!     ** N.B.                                **
!     ** TO BE USED (IF DESIRED) AS SURFACE  **
!     ** T & S SURFACE BOUNDARY CONDITION    **
!     **                                     **
!     *****************************************
!     *****************************************
!
      real(RLEN)               :: TSURF, SSURF
!
!     -----SURFACE SALINITY FLUX-----
!
      real(RLEN)               :: WSSURF
!
!     -----LATERAL ADVECTION FLUX FOR T & S-----
!
      real(RLEN), dimension(KB):: WTADV, WSADV
!              
      end module POM

! EOC
