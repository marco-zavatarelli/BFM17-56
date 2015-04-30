MODULE trcnxtbfm
!   !!======================================================================
!   !!                       ***  MODULE  trcnxt  ***
!   !! Ocean passive tracers:  time stepping on passives tracers
!   !!======================================================================
#if defined key_top 
!   !!----------------------------------------------------------------------
!   !!   trc_nxt     : time stepping on passive tracers
!   !!----------------------------------------------------------------------
!   !! * Modules used
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl_trc      ! Print control for debbuging
   USE trd_oce
   USE trdtra
   USE tranxt
   USE trcnam_trp
   USE c1d, only: lk_c1d
# if defined key_agrif
   USE agrif_top_interp
# endif

   IMPLICIT NONE
   PRIVATE

!   !! * Routine accessibility
   PUBLIC trc_nxt_bfm      ! routine called by trc_trp_bfm.F90
   PUBLIC trc_nxt_alloc    ! routine called by nemogcm.F90

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   r2dt

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_nxt_alloc()
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_nxt_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( r2dt(jpk), STAT=trc_nxt_alloc )
      !
      IF( trc_nxt_alloc /= 0 )   CALL ctl_warn('trc_nxt_alloc : failed to allocate array')
      !
   END FUNCTION trc_nxt_alloc


   SUBROUTINE trc_nxt_bfm( kt,m )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trcnxt  ***
      !!
      !! ** Purpose :   Compute the passive tracers fields at the 
      !!      next time-step from their temporal trends and swap the fields.
      !!      This is a modified version for the BFM to check that
      !!      Euler integration schemes are used 
      !!
      !! ** Method  :   Apply lateral boundary conditions on (ua,va) through
      !!      call to lbc_lnk routine
      !!   default:
      !!      arrays swap
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!         (trb) = (trn)
      !!
      !! ** Action  : - update trb, trn
      !!----------------------------------------------------------------------
      !! * Arguments

!   !! * Routine accessibility
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      INTEGER, INTENT( in ) ::   m          ! index of the BFM state variable
      !! * Local declarations
      INTEGER  ::   jk,jn   ! dummy loop indices
      REAL(wp) ::   zfact   ! temporary scalar
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  ztrdt
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_nxt_bfm')
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nxt_bfm : time stepping on BFM tracer',m
      ENDIF


      jn = 1 ! BFM uses only one tracer as a working array

      ! 0. Lateral boundary conditions on tra (T-point, unchanged sign)
      ! ---------------------------------============
      CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )

      ! Check that the BFM compliant integrations are used
      IF(.NOT.lk_c1d) THEN
        IF( ln_trczdf_exp .OR. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) THEN
           WRITE(numout,*)
           WRITE(numout,*) 'The following must all be set to FALSE with the BFM:'
           WRITE(numout,*) 'ln_trczdf_exp: ',ln_trczdf_exp
           WRITE(numout,*) 'ln_trcadv_cen2: ',ln_trcadv_cen2
           WRITE(numout,*) 'ln_trcadv_tvd: ',ln_trcadv_tvd
           CALL ctl_stop( ' BFM can only be run with implicit integration and Euler stepping ' )
        ENDIF
      ENDIF

      r2dt(:) =     rdttrc(:)           ! ( Euler step ) 

#if defined key_agrif
      ! Update tracer at AGRIF zoom boundaries
      IF( .NOT.Agrif_Root() )    CALL Agrif_Update_Trc( kt )      ! children only
#endif

!  case of smolar scheme or muscl (the only ones allowed with BFM)
      DO jk = 1, jpk
            trb(:,:,jk,jn) = tra(:,:,jk,jn)
            trn(:,:,jk,jn) = tra(:,:,jk,jn)
            tra(:,:,jk,jn) = 0.
      END DO


   END SUBROUTINE trc_nxt_bfm

!
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nxt_bfm( kt,m )
      INTEGER, INTENT(in) :: kt,m
      WRITE(*,*) 'trc_nxt_bfm: You should not have seen this print! error?', kt
   END SUBROUTINE trc_nxt_bfm
#endif
!   !!======================================================================
END MODULE trcnxtbfm
