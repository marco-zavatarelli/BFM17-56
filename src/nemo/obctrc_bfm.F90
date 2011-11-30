#include"cppdefs.h"
MODULE obctrc_bfm
   !!=================================================================================
   !!                       ***  MODULE  obctra  ***
   !! Ocean tracers:   Radiation of tracers on each open boundary
   !!=================================================================================

#if defined key_obcbfm
   !!---------------------------------------------------------------------------------
   !!   'key_obcbfm'   :                                      Open Boundary Conditions
   !!---------------------------------------------------------------------------------
   !!   obc_trc        : call the subroutine for each open boundary
   !!   obc_trc_east   : radiation of the east open boundary tracers
   !!   obc_trc_west   : radiation of the west open boundary tracers
   !!   obc_trc_north  : radiation of the north open boundary tracers
   !!   obc_trc_south  : radiation of the south open boundary tracers
   !!   Key changed by Tomas Lovato - 10/21/11 
   !!----------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! ???
   USE lbclnk          ! ???
   USE in_out_manager  ! I/O manager

   USE trc
   USE par_kind,ONLY   :wp
   USE trcobc_oce_bfm      ! ocean open boundary conditions
   IMPLICIT NONE

   !! * Accessibility
   PUBLIC obc_trc_bfm     ! routine called in trc_nxt_bfm.F90

   !! * Module variables
   INTEGER ::      & ! ... boundary space indices
      nib   = 1,   & ! nib   = boundary point
      nibm  = 2,   & ! nibm  = 1st interior point
      nibm2 = 3,   & ! nibm2 = 2nd interior point
                     ! ... boundary time indices
      nit   = 1,   & ! nit    = now
      nitm  = 2,   & ! nitm   = before
      nitm2 = 3      ! nitm2  = before-before
#if defined  key_mfs
   INTEGER :: jk,jj,ji      !  loop index
#endif

   REAL(wp) ::     &
      rtaue  , rtauw  , rtaun  , rtaus  ,  &  ! Boundary restoring coefficient
      rtauein, rtauwin, rtaunin, rtausin      ! Boundary restoring coefficient for inflow

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005)
   !! $Id: obctra.F90 1152 2008-06-26 14:11:13Z rblod $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_trc_bfm( kt,m )
      !!-------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_trc_bfm  ***
      !!
      !! ** Purpose :   Compute tracer fields  along the open boundaries.
      !!      This routine is called by the trcnxt_bfm.F routine
      !!        The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
      !!      and/or lp_obc_south allow the user to determine which boundary is an
      !!      open one (must be done in the param_obc.h90 file).
      !!
      !! Reference :
      !!   Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
      !!
      !!  History :
      !!        !  95-03 (J.-M. Molines) Original, SPEM
      !!        !  97-07 (G. Madec, J.-M. Molines) addition
      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) F90
      !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! ???
   USE lbclnk          ! ???
   USE in_out_manager  ! I/O manager


      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,m
      !! Local variables
      REAL(wp)              ::   zrdt

      ! 0. Local constant initialization

      IF( kt == nit000 .OR. ln_rstart) THEN
         ! ... Boundary restoring coefficient (NOTE, discarded the vertical dependence of rdttrc!)
         zrdt =  rdt * FLOAT(nn_dttrc)
         rtaue = 2. * zrdt / rn_trdpeob
         rtauw = 2. * zrdt / rn_trdpwob
         rtaun = 2. * zrdt / rn_trdpnob
         rtaus = 2. * zrdt / rn_trdpsob
         ! ... Boundary restoring coefficient for inflow ( all boundaries)
         rtauein = 2. * rdt / rn_trdpein
         rtauwin = 2. * rdt / rn_trdpwin
         rtaunin = 2. * rdt / rn_trdpnin
         rtausin = 2. * rdt / rn_trdpsin
      END IF

      IF( lp_obc_east  )   CALL obc_trc_bfm_east ( kt ,m)    ! East open boundary

      IF( lp_obc_west  )   CALL obc_trc_bfm_west ( kt ,m)    ! West open boundary

      IF( lp_obc_north )   CALL obc_trc_bfm_north( kt ,m)    ! North open boundary

      IF( lp_obc_south )   CALL obc_trc_bfm_south( kt ,m)    ! South open boundary

#if defined  key_mfs
      !-------------------------------------------------------------
      ! If MFS-MERSEA implementation 2 ocean corners (S-W ans N-W).
      ! Values at these points are replaced with inner values (numerical!)
      !-------------------------------------------------------------
        DO jk = 1, jpkm1
           do jj=njs0,njs1
              do ji=niw0,niw1
                 tra(ji,jj,jk,1)     =  tra(ji+1,jj,jk,1)
              enddo
           enddo
           do jj=njn0p1,njn1p1
                 do ji=niw0,niw1
                    tra(ji,jj,jk,1) =  tra(ji+1,jj,jk,1)
                 enddo
           enddo
        END DO
#endif

      IF( lk_mpp ) THEN                  !!bug ???
         IF( kt >= nit000+3 .AND. ln_rstart ) THEN
            CALL lbc_lnk( trb(:,:,:,1), 'T', 1. )
         END IF
         CALL lbc_lnk( tra(:,:,:,1), 'T', 1. )
      ENDIF

  END SUBROUTINE obc_trc_bfm

   SUBROUTINE obc_trc_bfm_east ( kt ,m)
      !!------------------------------------------------------------------------------
      !!                ***  SUBROUTINE obc_trc_east  ***
      !!
      !! ** Purpose :
      !!      Apply the radiation algorithm on east OBC tracers ta, sa using the
      !!      phase velocities calculated in obc_rad_east subroutine in obcrad.F90 module
      !!      If the logical lfbceast is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines)
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!------------------------------------------------------------------------------
      USE oce             ! ocean dynamics and tracers variables
      USE dom_oce         ! ocean space and time domain variables
      USE phycst          ! physical constants
      USE obc_oce         ! ocean open boundary conditions
      USE trcobc_oce_bfm      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE lbclnk          ! ???
      USE in_out_manager  ! I/O manager
      USE trc
      USE daymod,ONLY     :nday,nmonth
      use global_mem, only:LOGUNIT
      USE par_kind,ONLY   :wp
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,m

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
      INTEGER ::   &
      nlecto,   &  ! switch for the first read
      ntobc1,   &  ! first record used
      ntobc2,   &  ! second record used
      ntobc        ! number of time steps in OBC files

      ! 0. Interpolation of OBC data
      ! --------------------------------------------------------
      ! Time interpolation should be done here to avoid array
      ! duplication. Not yet implemented
      ! --------------------------------------------------------
      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
         DO jk = 1, jpkm1
            DO jj = nje0p1, nje1m1
               ij = jj -1 + njmpp
               trfoe(jj,jk) =  tredta(ij,jk,1,m) * temsk(jj,jk)
             END DO
         END DO
      ELSE
         CALL ctl_stop('Usage of external open boundary data not implemented with BFM',' routine obc_trc_bfm_east')
      ENDIF

      ! 1. First three time steps and more if ltrfbceast is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF(  (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. ltrfbceast ) THEN
         DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1. - temsk(jj,jk)) + &
                                 trfoe(jj,jk)*temsk(jj,jk)
               END DO
            END DO
         END DO

      ELSE

      ! 2. Beyond the fourth time step if ltrfbceast is .FALSE.
      ! -----------------------------------------------------

         ! Tracer radiation/advection
         ! ----------------------------------
         !
         !            nibm2      nibm      nib
         !              |   nibm  |   nib///|///
         !              |    |    |    |////|///
         !  jj   line --v----f----v----f----v---
         !              |    |    |    |////|///
         !                   |         |///   //
         !  jj   line   T    u    T    u/// T //
         !                   |         |///   //
         !              |    |    |    |////|///
         !  jj-1 line --v----f----v----f----v---
         !              |    |    |    |////|///
         !                jpieob-1    jpieob / ///
         !              |         |         |
         !           jpieob-1    jpieob     jpieob+1
         !
         ! ... radiative conditions + relaxation toward a climatology
         !     the phase velocity is the normal current velocity

         DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
         ! ... average of the normal velocity (normalized and bounded for stability)
                  z05cx = zrdt * ( 0.5 * ( uebnd(jj,jk,nibm,nit) +  &
                          uebnd(jj-1,jk,nibm,nit) ) ) / e1t(ji-1,jj)
                  z05cx = min( z05cx, 1. )
         ! ... z05cx=< 0, inflow  zin=0, ztau=rtauein
         !           > 0, outflow zin=1, ztau=rtaue
                  zin = sign( 1., z05cx )
                  zin = 0.5*( zin + abs(zin) )
         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
                  ztau = (1.-zin ) * rtauein  + zin * rtaue
                  z05cx = z05cx * zin
         ! ... update with radiative or climatological
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1. - temsk(jj,jk)) +           &
                                 temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
                                 * trebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
                                 * trebnd(jj,jk,nibm,nit ) + ztau * trfoe (jj,jk) ) &
                                 / (1. + z05cx)
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE obc_trc_bfm_east


   SUBROUTINE obc_trc_bfm_west ( kt ,m)
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_tra_west  ***
      !!
      !! ** Purpose :
      !!      Apply the advection/radiation algorithm on west OBC tracers
      !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines)
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      USE oce             ! ocean dynamics and tracers variables
      USE dom_oce         ! ocean space and time domain variables
      USE phycst          ! physical constants
      USE obc_oce         ! ocean open boundary conditions
      USE trcobc_oce_bfm      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE lbclnk          ! ???
      USE in_out_manager  ! I/O manager
      USE trc
      USE daymod,ONLY     :nday,nmonth
      USE par_kind,ONLY   :wp
      use global_mem, only:LOGUNIT

      INTEGER, INTENT( in ) ::   kt,m

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
            INTEGER ::   &
      nlecto,   &  ! switch for the first read
      ntobc1,   &  ! first record used
      ntobc2,   &  ! second record used
      ntobc        ! number of time steps in OBC files

      !!------------------------------------------------------------------------------


      ! 0. Interpolation of OBC data
      ! --------------------------------------------------------
      ! Time interpolation should be done here to avoid array
      ! duplication. Not yet implemented
      ! --------------------------------------------------------

      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
         DO jk = 1, jpkm1
            DO jj = njw0p1, njw1m1
               ij = jj -1 + njmpp
                  trfow(jj,jk) =  trwdta(ij,jk,1,m) * twmsk(jj,jk)
             END DO
         END DO
      ELSE
         CALL ctl_stop('Usage of external open boundary data not implemented with BFM',' routine obc_trc_bfm_west')
      ENDIF

      ! 1. First three time steps and more if ltrfbcwest is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. ltrfbcwest ) THEN

         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1. - twmsk(jj,jk)) + &
                                 trfow(jj,jk)*twmsk(jj,jk)
               END DO
            END DO
         END DO

      ELSE

      ! 2. Beyond the fourth time step if ltrfbcwest is .FALSE.
      ! -----------------------------------------------------

         ! Tracer radiation/advection
         ! ----------------------------------
         !
         !          nib       nibm     nibm2
         !     nib///|   nibm  |  nibm2  |
         !   ///|////|    |    |    |    |
         !   ---v----f----v----f----v----f-- jj   line
         !   ///|////|    |    |    |    |
         !   //   ///|         |         |
         !   // T ///u    T    u    T    u   jj   line
         !   //   ///|         |         |
         !   ///|////|    |    |    |    |
         !   ---v----f----v----f----v----f-- jj-1 line
         !   ///|////|    |    |    |    |
         !         jpiwob    jpiwob+1    jpiwob+2
         !      |         |         |
         !    jpiwob    jpiwob+1   jpiwob+2
         !
         ! ... radiative conditions + relaxation toward a climatology
         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
         ! ... average of the normal velocity (normalized and bounded for stability)
                  z05cx = zrdt * (  0.5 * ( uwbnd(jj,jk,nibm,nit) + &
                          uwbnd(jj-1,jk,nibm,nit) ) ) / e1t(ji+1,jj)
                  z05cx = max( z05cx, -1. )
         ! ... z05cx > 0, inflow  zin=0, ztau=rtauwin
         !           < 0, outflow zin=1, ztau=rtauw
                  zin = sign( 1., -1.* z05cx )
                  zin = 0.5*( zin + abs(zin) )
                  ztau = (1.-zin )*rtauwin + zin * rtauw
                  z05cx = z05cx * zin
         ! ... update with radiative or climatological
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1. - twmsk(jj,jk)) +           &
                                 twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
                                 * trwbnd(jj,jk,nib ,nitm,m) - 2.*z05cx              &
                                 * trwbnd(jj,jk,nibm,nit,m) + ztau * trfow (jj,jk) ) &
                                 / (1. - z05cx)
               END DO
            END DO
         END DO

      END IF
   END SUBROUTINE obc_trc_bfm_west


   SUBROUTINE obc_trc_bfm_north ( kt,m )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_trc_north  ***
      !!
      !! ** Purpose :
      !!      Apply the radiation algorithm on north OBC tracers ta, sa using the
      !!      phase velocities calculated in obc_rad_north subroutine in obcrad.F90 module
      !!      If the logical lfbcnorth is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines)
      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!------------------------------------------------------------------------------
      USE oce             ! ocean dynamics and tracers variables
      USE dom_oce         ! ocean space and time domain variables
      USE phycst          ! physical constants
      USE obc_oce         ! ocean open boundary conditions
      USE trcobc_oce_bfm      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE lbclnk          ! ???
      USE in_out_manager  ! I/O manager
      USE trc
      USE daymod,ONLY     :nday,nmonth
      USE par_kind,ONLY   :wp
      use global_mem, only:LOGUNIT

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,m
             INTEGER ::   &
      nlecto,   &  ! switch for the first read
      ntobc1,   &  ! first record used
      ntobc2,   &  ! second record used
      ntobc        ! number of time steps in OBC files


      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
      !!------------------------------------------------------------------------------

      ! 0. Interpolation of OBC data
      ! --------------------------------------------------------
      ! Time interpolation should be done here to avoid array
      ! duplication. Not yet implemented
      ! --------------------------------------------------------

      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
         DO jk = 1, jpkm1
            DO ji = nin0p1, nin1m1
               ii = ji -1 + nimpp
                  trfon(ji,jk) = trndta(ii,jk,1,m) * tnmsk(ji,jk)
            END DO
         END DO
     ELSE
         CALL ctl_stop('Usage of external open boundary data not implemented with BFM',' routine obc_trc_bfm_north')
     ENDIF

      ! 1. First three time steps and more if ltrfbcnorth is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. ltrfbcnorth ) THEN

         DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  tra(ji,jj,jk,1)= tra(ji,jj,jk,1) * (1.-tnmsk(ji,jk)) + &
                                tnmsk(ji,jk) * trfon(ji,jk)
               END DO
            END DO
         END DO

      ELSE

      ! 2. Beyond the fourth time step if ltrfbcnorth is .FALSE.
      ! -------------------------------------------------------

         ! Tracer radiation/advection
         ! ----------------------------------
         !
         !           ji-1   ji   ji   ji +1
         !             |
         !    nib //// u // T // u // T //   jpjnob + 1
         !        /////|//////////////////
         !    nib  ----f----v----f----v---   jpjnob
         !             |         |
         !      nibm-- u -- T -- u -- T --   jpjnob
         !             |         |
         !   nibm  ----f----v----f----v---  jpjnob-1
         !             |         |
         !     nibm2-- u -- T -- T -- T --  jpjnob-1
         !             |         |
         !   nibm2 ----f----v----f----v---  jpjnob-2
         !             |         |
         !
         ! ... radiative conditions + relaxation toward a climatology
         DO jj = fs_njn0+1, fs_njn1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
         ! ... j-phase speed ratio (from averaged of vtnbnd)
         !        (bounded by 1)
         ! ... average of the normal velocity (normalized and bounded for stability)
                  z05cx = zrdt * (  0.5 * ( vnbnd(jj,jk,nibm,nit) + &
                          vnbnd(jj-1,jk,nibm,nit) ) ) / e1t(ji+1,jj)
                  z05cx = min( z05cx, 1. )
         ! ... z05cx=< 0, inflow  zin=0, ztau=rtaunin
         !           > 0, outflow zin=1, ztau=rtaun
                  zin = sign( 1., z05cx )
                  zin = 0.5*( zin + abs(zin) )
         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
                  ztau = (1.-zin ) * rtaunin + zin * rtaun
                  z05cx = z05cx * zin
         ! ... update (ta,sa) with radiative or climatological (t, s)
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1.-tnmsk(ji,jk)) +             &
                                 tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
                                 * trnbnd(ji,jk,nib ,nitm,m) + 2.*z05cx              &
                                 * trnbnd(ji,jk,nibm,nit,m) + ztau * trfon (ji,jk) ) &
                                 / (1. + z05cx)
               END DO
            END DO
         END DO

      END IF

   END SUBROUTINE obc_trc_bfm_north


   SUBROUTINE obc_trc_bfm_south ( kt ,m)
      !!------------------------------------------------------------------------------
      !!                ***  SUBROUTINE obc_tra_south  ***
      !!
      !! ** Purpose :
      !!      Apply the /advectionradiation algorithm on south OBC tracers
      !!      If the logical lfbcsouth is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines)
      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
      !!------------------------------------------------------------------------------
      USE oce             ! ocean dynamics and tracers variables
      USE dom_oce         ! ocean space and time domain variables
      USE phycst          ! physical constants
      USE obc_oce         ! ocean open boundary conditions
      USE trcobc_oce_bfm      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE lbclnk          ! ???
      USE in_out_manager  ! I/O manager
      USE trc
      USE daymod,ONLY     :nday,nmonth
      USE par_kind,ONLY   :wp
      use global_mem, only:LOGUNIT


      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,m

      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
             INTEGER ::   &
      nlecto,   &  ! switch for the first read
      ntobc1,   &  ! first record used
      ntobc2,   &  ! second record used
      ntobc        ! number of time steps in OBC files

      !!------------------------------------------------------------------------------


      ! 0. Interpolation of OBC data
      ! --------------------------------------------------------
      ! Time interpolation should be done here to avoid array
      ! duplication. Not yet implemented
      ! --------------------------------------------------------

      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
         DO jk = 1, jpkm1
           DO ji = nis0p1, nis1m1
              ii = ji -1 + nimpp
                 trfos(ji,jk) =   trsdta(ii,jk,1,m) * tsmsk(ji,jk)
           END DO
         END DO
      ELSE
         CALL ctl_stop('Usage of external open boundary data not implemented with BFM',' routine obc_trc_bfm_north')
      ENDIF

      ! 1. First three time steps and more if lfbcsouth is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF(  (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. ltrfbcsouth) THEN

         DO jj = fs_njs0, fs_njs1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  tra(ji,jj,jk,1)= tra(ji,jj,jk,1) * (1.-tsmsk(ji,jk)) + &
                                tsmsk(ji,jk) * trfos(ji,jk)
               END DO
            END DO
         END DO

      ELSE

      ! 2. Beyond the fourth time step if ltrfbcsouth is .FALSE.
      ! -------------------------------------------------------

         ! Temperature and salinity radiation
         ! ----------------------------------
         !
         !           ji-1   ji   ji   ji +1
         !             |         |
         !   nibm2 ----f----v----f----v---   jpjsob+2
         !             |         |
         !   nibm2 --  u -- T -- u -- T --   jpjsob+2
         !             |         |
         !   nibm  ----f----v----f----v---   jpjsob+1
         !             |         |
         !    nibm --  u -- T -- T -- T --   jpjsob+1
         !             |         |
         !   nib  -----f----v----f----v---   jpjsob
         !       //////|/////////|////////
         !    nib //// u // T // u // T //   jpjsob
         !
         !... radiative conditions + relaxation toward a climatology

         DO jj = fs_njs0, fs_njs1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
         ! ... average of the normal velocity (normalized and bounded for stability)
                  z05cx = ( 0.5 * ( vsbnd(ji,jk) + vsbnd(ji-1,jk) ) ) / e2t(ji,jj+1)
                  z05cx = max( z05cx, -1. )
         !... z05cx > 0, inflow  zin=0, ztau=rtausin
         !          < 0, outflow zin=1, ztau=rtaus
                  zin = sign( 1., -1.* z05cx )
                  zin = 0.5*( zin + abs(zin) )
                  ztau = (1.-zin ) * rtausin + zin * rtaus
                  z05cx = z05cx * zin
         !... update with radiative or climatological
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1.-tsmsk(ji,jk)) +             &
                                 tsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
                                 * trsbnd(ji,jk,nib ,nitm,m) - 2.*z05cx              &
                                 * trsbnd(ji,jk,nibm,nit,m) + ztau * trfos (ji,jk) ) &
                                 / (1. - z05cx)
               END DO
            END DO
         END DO

      END IF
   END SUBROUTINE obc_trc_bfm_south

#endif
END MODULE obctrc_bfm

