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
   USE trcnam_trp      ! pasive tracers transport
   USE trdmod_oce
   USE trdtra
   USE tranxt
   USE prtctl_trc      ! Print control for debbuging
#if defined key_agrif
   USE agrif_top_update
   USE agrif_top_interp
#endif
   IMPLICIT NONE
   PRIVATE

!   !! * Routine accessibility
   PUBLIC trc_nxt_bfm      ! routine called by trc_trp_bfm.F90
   PUBLIC trc_nxt_alloc    ! routine called by nemogcm.F90

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   r2dt

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnxt.F90 2690 2011-03-15 15:27:46Z gm $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
!
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
      !!      This is a modified version for the BFM to include
      !!      open boundary conditions
      !!
      !! ** Method  :   Apply lateral boundary conditions on (ua,va) through
      !!      call to lbc_lnk routine
      !!   default:
      !!      arrays swap
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!         (trb) = (trn)
      !!
      !!   For Arakawa or TVD Scheme :
      !!      A Asselin time filter applied on now tracers (trn) to avoid
      !!      the divergence of two consecutive time-steps and tr arrays
      !!      to prepare the next time_step:
      !!         (trb) = (trn) + atfp [ (trb) + (tra) - 2 (trn) ]
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!
      !!
      !! ** Action  : - update trb, trn
      !!
      !! History :
      !!   7.0  !  91-11  (G. Madec)  Original code
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  95-02  (M. Levy)   passive tracers
      !!        !  96-02  (G. Madec & M. Imbard)  opa release 8.0
      !!   8.0  !  96-04  (A. Weaver)  Euler forward step
      !!   8.2  !  99-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
      !!   9.0  !  04-03  (C. Ethe) passive tracers
      !!----------------------------------------------------------------------
      !! * Arguments
      USE oce_trc         ! ocean dynamics and tracers variables
      USE trc             ! ocean passive tracers variables
      USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
      USE trcnam_trp      ! pasive tracers transport
      USE prtctl_trc      ! Print control for debbuging
#ifdef key_obcbfm
      USE obctrc_bfm
#endif

!   !! * Routine accessibility
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      INTEGER, INTENT( in ) ::   m          ! index of the BFM state variable
      !! * Local declarations
      INTEGER  ::   jk,jn   ! dummy loop indices
      REAL(wp) ::   zfact   ! temporary scalar
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  ztrdt
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

       IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nxt_bfm : time stepping on BFM tracer',m
       ENDIF


      jn = 1 ! BFM uses only one tracer as a working array

      ! 0. Lateral boundary conditions on tra (T-point, unchanged sign)
      ! ---------------------------------============
      CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )

      ! set time step size (Euler/Leapfrog)
      IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
         IF( neuler == 0 .AND. kt ==  nit000) THEN  ;  r2dt(:) =     rdttrc(:)   ! at nit000             (Euler)
         ELSEIF( kt <= nit000 + 1 )           THEN  ;  r2dt(:) = 2.* rdttrc(:)   ! at nit000 or nit000+1 (Leapfrog)
         ENDIF
      ELSE
         r2dt(:) =     rdttrc(:)           ! ( Euler step ) 
      ENDIF

      ! trends computation initialisation
      IF( l_trdtrc )  THEN
         ALLOCATE( ztrdt(jpi,jpj,jpk,1) )  !* store now fields before applying the Asselin filter
         ztrdt(:,:,:,:)  = trn(:,:,:,:)
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpk                                   ! Horizontal slab
         !                                             ! ===============
         ! 1. Leap-frog scheme (only in explicit case, otherwise the 
         ! -------------------  time stepping is already done in trczdf)
         IF( ln_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) THEN
            zfact = 2. * rdttra(jk) * FLOAT(nn_dttrc)
            IF( neuler == 0 .AND. kt == nit000 ) zfact = rdttra(jk) * FLOAT(nn_dttrc)
            tra(:,:,jk,jn) = ( trb(:,:,jk,jn) + zfact * tra(:,:,jk,jn) ) * tmask(:,:,jk)
         ENDIF

      END DO

#if defined key_obcbfm
         ! Update tracers on open boundaries.
      CALL obc_trc_bfm(kt,m)
#endif

#if defined key_agrif
      ! Update tracer at AGRIF zoom boundaries
      IF( .NOT.Agrif_Root() )    CALL Agrif_Update_Trc( kt )      ! children only
#endif

      DO jk = 1, jpk

         ! 2. Time filter and swap of arrays
         ! ---------------------------------
         IF ( ln_trcadv_cen2 .OR. ln_trcadv_tvd  ) THEN

            IF( neuler == 0 .AND. kt == nit000 ) THEN
               trb(:,:,jk,jn) = trn(:,:,jk,jn)
               trn(:,:,jk,jn) = tra(:,:,jk,jn)
               tra(:,:,jk,jn) = 0.
            ELSE
               trb(:,:,jk,jn) = atfp  * ( trb(:,:,jk,jn) + tra(:,:,jk,jn) ) + atfp1 * trn(:,:,jk,jn)
               trn(:,:,jk,jn) = tra(:,:,jk,jn)
               tra(:,:,jk,jn) = 0.
            ENDIF

         ELSE
!  case of smolar scheme or muscl
            trb(:,:,jk,jn) = tra(:,:,jk,jn)
            trn(:,:,jk,jn) = tra(:,:,jk,jn)
            tra(:,:,jk,jn) = 0.
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! trends computation
      IF( l_trdtrc ) THEN                                      ! trends
         DO jk = 1, jpkm1
               zfact = 1.e0 / r2dt(jk)
               ztrdt(:,:,jk,jn) = ( trb(:,:,jk,jn) - ztrdt(:,:,jk,jn) ) * zfact
               CALL trd_tra( kt, 'TRC', jn, jptra_trd_atf, ztrdt )
         END DO
         DEALLOCATE( ztrdt )
      END IF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nxt')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF


   END SUBROUTINE trc_nxt_bfm
!
#else
!   !!----------------------------------------------------------------------
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
