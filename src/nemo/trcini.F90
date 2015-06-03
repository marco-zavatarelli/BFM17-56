MODULE trcini
   !!======================================================================
   !!                         ***  MODULE trcini  ***
   !! TOP :   Manage the passive tracer initialization
   !!======================================================================
   !! History :   -   ! 1991-03 (O. Marti)  original code
   !!            1.0  ! 2005-03 (O. Aumont, A. El Moussaoui) F90
   !!            2.0  ! 2005-10 (C. Ethe, G. Madec) revised architecture
   !!            4.0  ! 2011-01 (A. R. Porter, STFC Daresbury) dynamical allocation
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_init  :   Initialization for passive tracer
   !!   top_alloc :   allocate the TOP arrays
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables
   USE trcrst          ! passive tracers restart
   USE trcnam          ! Namelist read
   USE trcini_cfc      ! CFC      initialisation
   USE trcini_pisces   ! PISCES   initialisation
   USE trcini_c14b     ! C14 bomb initialisation
   USE trcini_my_trc   ! MY_TRC   initialisation
   USE trcdta          ! initialisation from files
   USE daymod          ! calendar manager
   USE zpshde          ! partial step: hor. derivative   (zps_hde routine)
   USE prtctl_trc      ! Print control passive tracers (prt_ctl_trc_init routine)
   USE trcsub          ! variables to substep passive tracers
   USE lib_mpp         ! distribued memory computing library
   USE sbc_oce
   ! BFM
   USE mem_PAR, ONLY : LightPeriodFlag
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_init   ! called by opa

    !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2011)
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
      IF( nn_timing == 1 )   CALL timing_start('trc_init')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_init : initial set up of the passive tracers'
      IF(lwp) WRITE(numout,*) '           Using the Biogeochemical Flux Model (BFM)'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      CALL top_alloc()              ! allocate TOP arrays

      IF( nn_cla == 1 )   &
         &       CALL ctl_stop( ' Cross Land Advection not yet implemented with passive tracer ; nn_cla must be 0' )

      CALL trc_nam                  ! read passive tracers namelists
      !                             ! masked grid volume
      DO jk = 1, jpk
         cvol(:,:,jk) = e1e2t(:,:) * fse3t(:,:,jk) * tmask(:,:,jk)
      END DO
      IF( lk_degrad ) cvol(:,:,:) = cvol(:,:,:) * facvol(:,:,:)      ! degrad option: reduction by facvol
      !                                                              ! total volume of the ocean 
      areatot = glob_sum( cvol(:,:,:) )


      IF ( MOD(nitend,nn_dttrc) /= 0 ) &
                 CALL ctl_stop( 'The length of the experiment (nn_itend) is not a multiple of the tracer step freq. (nn_dttrc)')

      CALL trc_ini_bfm              ! Initialize BFM tracers

      IF( lk_offline )  THEN
           neuler = 0                  ! Set time-step indicator at nit000 (euler)
           CALL day_init               ! set calendar
      ENDIF
      ! Initialize temporary NEMO arrays for tracer transport
      trn(:,:,:,:) = 0._wp
      tra(:,:,:,:) = 0._wp
      trb(:,:,:,:) = trn(:,:,:,:)
      ! Note: there is no initial computation of partial step because it is done
      ! runtime in trc_trp_bfm

      !
      IF( nn_dttrc /= 1 )        CALL trc_sub_ini      ! Initialize variables for substepping passive tracers
      !
      ! check consistency of light paramterizations
      if ( ln_dm2dc .AND. LightPeriodFlag .NE. 1) &
                 CALL ctl_stop( ' The diurnal cycle (ln_dm2dc) is not compatible with the BFM LightPeriodFlag = 1.' )

      IF( nn_timing == 1 )   CALL timing_stop('trc_init')
      
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
      USE trdtrc_oce    , ONLY:   trd_trc_oce_alloc
#if defined key_trdmxl_trc 
      USE trdmxl_trc    , ONLY:   trd_mxl_trc_alloc
#endif
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        trc_adv_alloc()          ! Start of TOP-related alloc routines...
      ierr = ierr + trc_alloc    ()
      ierr = ierr + trc_nxt_alloc()
      ierr = ierr + trc_zdf_alloc()
      ierr = ierr + trd_trc_oce_alloc()
#if defined key_trdmxl_trc 
      ierr = ierr + trd_mxl_trc_alloc()
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
