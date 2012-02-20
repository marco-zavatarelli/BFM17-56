MODULE trcini
   !!======================================================================
   !!                         ***  MODULE trcini  ***
   !! TOP :   Manage the passive tracer initialization
   !!         This is a special version of trcini to be used with the BFM
   !!======================================================================
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_init  :   Initialization for passive tracer
   !!   top_alloc :   allocate the TOP arrays
   !!----------------------------------------------------------------------
   USE oce_trc
   USE trc
   USE trcrst
   USE trcnam          ! Namelist read
   USE trcini_cfc      ! CFC      initialisation
   USE trcini_lobster  ! LOBSTER  initialisation
   USE trcini_pisces   ! PISCES   initialisation
   USE trcini_c14b     ! C14 bomb initialisation
   USE trcini_my_trc   ! MY_TRC   initialisation
   USE trcdta   
   USE daymod
   USE zpshde          ! partial step: hor. derivative   (zps_hde routine)
   USE prtctl_trc      ! Print control passive tracers (prt_ctl_trc_init routine)
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_init   ! called by opa

    !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2011)
   !! $Id: trcini.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE trc_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_init  ***
      !!
      !! ** Purpose :   Initialization of the passive tracer fields 
      !!
      !! ** Method  : - read namelist
      !!              - control the consistancy 
      !!              - compute specific initialisations
      !!              - set initial tracer fields (either read restart 
      !!                or read data or analytical formulation
      !!---------------------------------------------------------------------
      INTEGER ::   jk, jn    ! dummy loop indices
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_init : initial set up of the passive tracers'
      IF(lwp) WRITE(numout,*) '           Using the Biogeochemical Flux Model (BFM'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      CALL top_alloc()              ! allocate TOP arrays

      !                             ! masked grid volume
      DO jk = 1, jpk
         cvol(:,:,jk) = e1t(:,:) * e2t(:,:) * fse3t(:,:,jk) * tmask(:,:,jk) 
      END DO

      !                             ! total volume of the ocean
#if ! defined key_degrad
      areatot = glob_sum( cvol(:,:,:) )
#else
      areatot = glob_sum( cvol(:,:,:) * facvol(:,:,:) )  ! degrad option: reduction by facvol
#endif

      IF( ln_dm2dc )      &
         &       CALL ctl_stop( ' The diurnal cycle is not compatible with PISCES or LOBSTER or BFM ' )

      IF( nn_cla == 1 )   &
         &       CALL ctl_stop( ' Cross Land Advection not yet implemented with passive tracer ; nn_cla must be 0' )

      CALL trc_nam                  ! read passive tracers namelists
      CALL trc_ini_bfm              ! Initialize BFM tracers

      IF( lk_offline )  THEN
           neuler = 0                  ! Set time-step indicator at nit000 (euler)
           CALL day_init               ! set calendar
      ENDIF
      trb(:,:,:,:) = trn(:,:,:,:)
      ! 
 
      tra(:,:,:,:) = 0._wp
      
   END SUBROUTINE trc_init


   SUBROUTINE top_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE top_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!----------------------------------------------------------------------
      USE trcadv        , ONLY:   trc_adv_alloc          ! TOP-related alloc routines...
      USE trc           , ONLY:   trc_alloc
      USE trcnxtbfm     , ONLY:   trc_nxt_alloc
      USE trczdf        , ONLY:   trc_zdf_alloc
      USE trdmod_trc_oce, ONLY:   trd_mod_trc_oce_alloc
#if ! defined key_iomput
      USE trcdia        , ONLY:   trc_dia_alloc
#endif
#if defined key_trcdmp 
      USE trcdmp        , ONLY:   trc_dmp_alloc
#endif
#if defined key_dtatrc
      USE trcdta        , ONLY:   trc_dta_alloc
#endif
#if defined key_trdmld_trc   ||   defined key_esopa
      USE trdmld_trc    , ONLY:   trd_mld_trc_alloc
#endif
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        trc_adv_alloc()          ! Start of TOP-related alloc routines...
      ierr = ierr + trc_alloc    ()
      ierr = ierr + trc_nxt_alloc()
      ierr = ierr + trc_zdf_alloc()
      ierr = ierr + trd_mod_trc_oce_alloc()
#if ! defined key_iomput
      ierr = ierr + trc_dia_alloc()
#endif
#if defined key_trcdmp 
      ierr = ierr + trc_dmp_alloc()
#endif
#if defined key_dtatrc
      ierr = ierr + trc_dta_alloc()
#endif
#if defined key_trdmld_trc   ||   defined key_esopa
      ierr = ierr + trd_mld_trc_alloc()
#endif
      !
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'top_alloc : unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE top_alloc

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_init                      ! Dummy routine   
   END SUBROUTINE trc_init
#endif

   !!======================================================================
END MODULE trcini
