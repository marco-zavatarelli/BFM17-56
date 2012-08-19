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

#if defined key_mfs && defined key_my_trc
      ! Set boundary conditions for MED equal to first interion line of the domain
      ! TOM: zero gradient condition -> to be put into obc_trc_bfm
      CALL obc_trc_bfm_tom(kt,m)
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

   !!----------------------------------------------------------------------------
   !! Subroutine to set BFM OBC zero gradient condition at the open boundary
   !! TOM: to be placed in obctrc_bfm 
   !!----------------------------------------------------------------------------
#ifdef key_obc
   SUBROUTINE obc_trc_bfm_tom(kt,m)
   !!----------------------------------------------------------------------------
   ! Purpose : Used to set OBC equal to internal points as OBC is not yet present 
   ! This Subroutine is basedon the old structure of bfm obc treatment
   ! NOTe that all boudaries are asumed to be solid 
   ! Tomas Lovato
   !!----------------------------------------------------------------------------
   ! include file for substituting indeces for open boundaries
#  include "obc_vectopt_loop_substitute.h90"

      USE oce             ! ocean dynamics and tracers variables
      USE dom_oce         ! ocean space and time domain variables
      USE phycst          ! physical constants
      USE obc_oce         ! ocean open boundary conditions
!      USE trcobc_oce_bfm      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE in_out_manager  ! I/O manager
!      USE daymod,ONLY     :nday,nmonth
      use global_mem, only:LOGUNIT
      USE par_kind,ONLY   :wp
      USE api_bfm
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,m

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices

      ! check consitency of indeces
!     if ( kt == 1 ) then
!      WRITE(numout,*) 'Enter OBC tom for BFM tracer: ',m
!      WRITE(numout,*) 'Check indexes for ', narea
!      WRITE(numout,*) 'Flag: ',lp_obc_west,' West  :',fs_niw0, fs_niw1
!      WRITE(numout,*) 'Flag: ',lp_obc_east,' East  :',fs_nie0, fs_nie1
!      WRITE(numout,*) 'Flag: ',lp_obc_north,' North :',fs_njn0, fs_njn1
!      WRITE(numout,*) 'Flag: ',lp_obc_south,' South :',fs_njs0, fs_njs1
!     endif
      ! Set Western boundary value
     if ( lp_obc_west ) then
         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                if ( twmsk(ji,jk) == 1 ) then
                   ! set equal to first inner domain cell (no gradient)
                   if ( tra(ji+1,jj,jk,1) /= 0 ) then
                      tra(ji,jj,jk,1) = tra(ji+1,jj,jk,1)
                   ! set equal to the one of the upper level
                   elseif ( tra(ji+1,jj,jk,1) == 0 ) then
                      tra(ji,jj,jk,1) = tra(ji,jj,jk-1,1)
!         if (kt == 1) WRITE(numout,*) 'Step :',kt,' Fault at ',narea,' East.',ji,jj,jk,' Used  ', tra(ji,jj,jk,1)
                   endif
                endif
!                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1. - twmsk(jj,jk)) + &
!                                 tfow(jj,jk)*twmsk(jj,jk)
             END DO
         END DO
      END DO
     endif
      ! Set Eastern boundary value
     if ( lp_obc_east ) then
      DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
         DO jk = 1, jpkm1
             DO jj = 1, jpj
                if ( temsk(ji,jk) == 1 ) then
                   ! set equal to first inner domain cell (no gradient)
                   if ( tra(ji-1,jj,jk,1) /= 0) then
                      tra(ji,jj,jk,1) = tra(ji-1,jj,jk,1)
                   ! set equal to the one of the upper level
                   elseif ( tra(ji-1,jj,jk,1) == 0 ) then
                      tra(ji,jj,jk,1) = tra(ji,jj,jk-1,1)
!         if (kt == 1) WRITE(numout,*) 'Step :',kt,' Fault at ',narea,' East.',ji,jj,jk,' Used  ', tra(ji,jj,jk,1)
                   endif
                endif 
!                tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1. - temsk(jj,jk)) + &
!                                    trfoe(jj,jk)*temsk(jj,jk)
             END DO
         END DO
      END DO
     endif
   ! Set northern boundary value
     if ( lp_obc_north ) then
      DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
          DO jk = 1, jpkm1
              DO ji = 1, jpi
                 if ( tnmsk(ji,jk) == 1 ) then 
                   ! set equal to first inner domain cell (no gradient)
                    if ( tra(ji,jj-1,jk,1) /= 0) then
                       tra(ji,jj,jk,1) = tra(ji,jj-1,jk,1)
                   ! set equal to the one of the upper level
                    elseif ( tra(ji,jj-1,jk,1) == 0 ) then
                       tra(ji,jj,jk,1) = tra(ji,jj,jk-1,1)
!      if (kt == 1) WRITE(numout,*) 'Step :',kt,' Fault at ',narea,' North.',ji,jj,jk,' Used  ', tra(ji,jj,jk,1)
                    endif
                 endif
!                 tra(ji,jj,jk,1)= tra(ji,jj,jk,1) * (1.-tnmsk(ji,jk)) + &
!                               tnmsk(ji,jk) * trfon(ji,jk)
              END DO
          END DO
      END DO
     endif
   ! Set southtern boundary value
     if (lp_obc_south) then
      DO jj = fs_njs0, fs_njs1  ! Vector opt.
          DO jk = 1, jpkm1
             DO ji = 1, jpi
                if ( tsmsk(ji,jk) == 1 ) then
                   ! set equal to first inner domain cell (no gradient)
                   if ( tra(ji,jj+1,jk,1) /= 0 )  then
                      tra(ji,jj,jk,1) = tra(ji,jj+1,jk,1)
                   ! set equal to the one of the upper level
                   elseif ( tra(ji,jj+1,jk,1) == 0 ) then 
                      tra(ji,jj,jk,1) = tra(ji,jj,jk-1,1)
!         if (kt == 1) WRITE(numout,*) 'Step :',kt,' Fault at ',narea,' South.',ji,jj,jk,' Used  ', tra(ji,jj,jk,1)
                   endif
                endif
!                tra(ji,jj,jk,1)= tra(ji,jj,jk,1) * (1.-tsmsk(ji,jk)) + &
!                              tsmsk(ji,jk) * trfos(ji,jk)
             END DO
         END DO
      END DO
    endif

      !-------------------------------------------------------------
      ! Based on MFS-MERSEA implementation 2 ocean corners (S-W ans N-W).
      ! Tracer data in these points are replaced using inner values (numerical!)
      !-------------------------------------------------------------

        DO jk = 1, jpkm1
           do jj=njs0,njs1
              do ji=niw0,niw1
                    tra(ji,jj,jk,1) = tra(ji+1,jj,jk,1)
!                 if(jk==1)write(numout,*)'MFS CHECK OBC SOUTH: JI wanted=2, actual=',ji
!                 if(jk==1)write(numout,*)'MFS CHECK: JJwanted=2, actual=',jj
                 ! ta(ji,jj,jk)     =  ta(ji+1,jj,jk)
              enddo
           enddo
           do jj=njn0p1,njn1p1
                 do ji=niw0,niw1
                    tra(ji,jj,jk,1) = tra(ji+1,jj,jk,1)
!                 if(jk==1)write(numout,*)'MFS CHECK OBC NORTH: JI wanted=2, actual=',ji
!                 if(jk==1)write(numout,*)'MFS CHECK: JJ wanted=272, actual=',jj
                 ! ta(ji,jj,jk) =  ta(ji+1,jj,jk)
                 enddo
           enddo
        END DO
      ! Link lbc values to global domain state vector
      IF( lk_mpp ) THEN                  !!bug ???
         IF( kt >= nit000+3 .AND. ln_rstart ) THEN
            CALL lbc_lnk( trb(:,:,:,1), 'T', 1. )
         END IF
         CALL lbc_lnk( tra(:,:,:,1), 'T', 1. )
      ENDIF


!   if (kt == 1) WRITE(numout,*) 'OBC tom done.'

   END SUBROUTINE obc_trc_bfm_tom
#endif
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
