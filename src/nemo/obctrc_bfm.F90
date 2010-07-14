#include"cppdefs.h"
#include "obc_vectopt_loop_substitute.h90"

!!---------------------------------------------------------------------------------
!! this routine is compiled only if key_obc is defined
!!---------------------------------------------------------------------------------
#ifdef key_obc
   SUBROUTINE obctrc_bfm( kt,m )
      !!-------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_trc  ***
      !!                    
      !! ** Purpose :   Compute tracer fields (t,s) along the open boundaries.
      !!      This routine is called by the tranxt.F routine and updates ta,sa
      !!      which are the actual temperature and salinity fields.
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
   USE trcobc_oce      ! ocean open boundary conditions
   USE lib_mpp         ! ???
   USE lbclnk          ! ???
   USE in_out_manager  ! I/O manager
   USE trc
   USE daymod,ONLY     :nday,nmonth
   USE par_kind,ONLY   :wp

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt,m
      INTEGER ::   ji, jj, jk
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   itimo, iman, imois
      INTEGER ::   i15
      REAL(wp) ::   zxy
      INTEGER ::   &
      nlecto,   &  ! switch for the first read
      ntobc1,   &  ! first record used
      ntobc2,   &  ! second record used
      ntobc        ! number of time steps in OBC files 
      REAL(wp), DIMENSION(:), ALLOCATABLE :: tcobc      ! time_counter variable of BCs


      ! 0. Local constant initialization

      IF( kt == nittrc000 .OR. ln_rstart) THEN
         ! ... Boundary restoring coefficient
         rtaue = 2. * rdt / rdpeob
         rtauw = 2. * rdt / rdpwob
         rtaun = 2. * rdt / rdpnob
         rtaus = 2. * rdt / rdpsob
         ! ... Boundary restoring coefficient for inflow ( all boundaries)
         rtauein = 2. * rdt / rdpein 
         rtauwin = 2. * rdt / rdpwin
         rtaunin = 2. * rdt / rdpnin
         rtausin = 2. * rdt / rdpsin 
      END IF
!GELSOMINA
!      IF( ntrcobc_dta == 0 )   THEN
!         itimo = 1
!         zxy   = 0.
!      ELSE
!         IF( ntobc == 1 )   THEN
!            itimo = 1
!         ELSE IF( ntobc == 12 )   THEN      !   BC are monthly   
!            ! we assume we have climatology in that case
!            iman  = 12
!            i15   = nday / 16
!            imois = nmonth + i15 - 1
!            IF( imois == 0 )   imois = iman
!            itimo = imois
!         ELSE
!            iman  = ntobc
!            itimo = FLOOR( kt*rdt / (tcobc(2)-tcobc(1)) )
!            isrel = kt*rdt
!         ENDIF
!      ENDIF
!GELSOMINA

      IF( lp_obc_east  )   CALL obc_trc_bfm_east ( kt ,m)    ! East open boundary

      IF( lp_obc_west  )   CALL obc_trc_bfm_west ( kt ,m)    ! West open boundary

      IF( lp_obc_north )   CALL obc_trc_bfm_north( kt ,m)    ! North open boundary

      IF( lp_obc_south )   CALL obc_trc_bfm_south( kt ,m)    ! South open boundary

#if defined  key_mfs
      !-------------------------------------------------------------
      ! If MFS-MERSEA implementation 2 ocean corners (S-W ans N-W).
      ! Temperature and salinity data in these points are replaced
      ! using inner values (numerical!)
      !-------------------------------------------------------------
        DO jk = 1, jpkm1
!ANTO Enrico+ 
           do jj=njs0,njs1
              do ji=niw0,niw1
!                 if(jk==1)write(numout,*)'MFS CHECK OBC SOUTH: JI wanted=2, actual=',ji
!                 if(jk==1)write(numout,*)'MFS CHECK: JJwanted=2, actual=',jj
                 tra(ji,jj,jk,1)     =  tra(ji+1,jj,jk,1)
              enddo   
           enddo   
           do jj=njn0p1,njn1p1
                 do ji=niw0,niw1
!                 if(jk==1)write(numout,*)'MFS CHECK OBC NORTH: JI wanted=2, actual=',ji
!                 if(jk==1)write(numout,*)'MFS CHECK: JJ wanted=252, actual=',jj
                    tra(ji,jj,jk,1) =  tra(ji+1,jj,jk,1)
                 enddo
           enddo
!ANTO -
        END DO
#endif

      IF( lk_mpp ) THEN                  !!bug ???
         IF( kt >= nit000+3 .AND. ln_rstart ) THEN
            CALL lbc_lnk( trb(:,:,:,1), 'T', 1. )
         END IF
         CALL lbc_lnk( tra(:,:,:,1), 'T', 1. )
      ENDIF
  END SUBROUTINE obctrc_bfm

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
      USE trcobc_oce      ! ocean open boundary conditions
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

     REAL(wp), DIMENSION(:), ALLOCATABLE :: tcobc      ! time_counter variable of BCs

      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if ltrfbceast is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------
! Linear interpolation of BCs to current time step
!      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
!         zxy = 0.

!ATTENZIONE PER ORA OK MA DA CAMBIARE QUANDO SI AVRANNO I DATI SULLE OBC
!VARIABILI NON DEFINITE ERRORE SIG FAULT.In obcdta.F90 sono definite e 'lette' nei file delle OBC!!!!
!      IF( ntobc == 12 .OR. ntrobcdta ==1)   THEN
!         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!      ELSE 
!         zxy = (tcobc(ntobc1)-FLOAT(isrel))/(tcobc(ntobc1)-tcobc(ntobc2))
!      ENDIF

! Per evitare la duplicazione degli array tredta
!       DO m = 1,NO_D3_BOX_STATES
!         IF (D3STATEOBC(m)/= 0.) THEN
         IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
               DO jk = 1, jpkm1
                 DO jj = nje0p1, nje1m1
               ij = jj -1 + njmpp
               trfoe(jj,jk) =  tredta(ij,jk,1,m) * temsk(jj,jk)
                END DO
              END DO

          ELSE
!               DO jk = 1, jpkm1
!                 DO jj = nje0p1, nje1m1
!               ij = jj -1 + njmpp
!               trfoe(jj,jk) =  ( zxy * tredta(ij,jk,2,m) + &
!                  &           (1.-zxy) * tredta(ij,jk,1,m) ) * temsk(jj,jk)
!                END DO
!              END DO
       ENDIF
      
      IF(  (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. lp_obc_east ) THEN
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

         ! Temperature and salinity radiation
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
         !     the phase velocity is taken as the phase velocity of the tangen-
         !     tial velocity (here vn), which have been saved in (u_cxebnd,v_cxebnd)
         ! ... (jpjedp1, jpjefm1), jpieob+1
         DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
         ! ... i-phase speed ratio (from averaged of v_cxebnd)
                  z05cx = ( 0.5 * ( v_cxebnd(jj,jk) + v_cxebnd(jj-1,jk) ) ) / e1t(ji-1,jj)
                  z05cx = min( z05cx, 1. )
         ! ... z05cx=< 0, inflow  zin=0, ztau=1    
         !           > 0, outflow zin=1, ztau=rtaue
                  zin = sign( 1., z05cx )
                  zin = 0.5*( zin + abs(zin) )
         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
                  ztau = (1.-zin ) * rtauein  + zin * rtaue
                  z05cx = z05cx * zin
         ! ... update ( ta, sa ) with radiative or climatological (t, s)
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1. - temsk(jj,jk)) +           & 
                                 temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
                                 * trebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
                                 * trebnd(jj,jk,nibm,nit ) + ztau * trfoe (jj,jk) ) &
                                 / (1. + z05cx)
               END DO
            END DO
         END DO

      END IF
!     END IF
!     END DO
   END SUBROUTINE obc_trc_bfm_east


   SUBROUTINE obc_trc_bfm_west ( kt ,m)
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_tra_west  ***
      !!           
      !! ** Purpose :
      !!      Apply the radiation algorithm on west OBC tracers ta, sa using the 
      !!      phase velocities calculated in obc_rad_west subroutine in obcrad.F90 module
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
      USE trcobc_oce      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE lbclnk          ! ???
      USE in_out_manager  ! I/O manager
      USE trc
      USE obctrc
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

       REAL(wp), DIMENSION(:), ALLOCATABLE :: tcobc      ! time_counter variable of BCs

      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if ltrfbcwest is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------
! Linear interpolation of BCs to current time step
!      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
!         zxy = 0.
!ATTENZIONE PER ORA OK MA DA CAMBIARE QUANDO SI AVRANNO I DATI SULLE OBC
!      IF( ntobc == 12 .OR. ntrcobc_dta == 1)   THEN
!        zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!      ELSE
!         zxy = (tcobc(ntobc1)-FLOAT(isrel))/(tcobc(ntobc1)-tcobc(ntobc2))
!      ENDIF

! Per evitare la duplicazione degli array trwdta
!       DO m = 1,NO_D3_BOX_STATES
!         IF (D3STATEOBC(m)/= 0.) THEN
       IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
         DO jk = 1, jpkm1
            DO jj = njw0p1, njw1m1
               ij = jj -1 + njmpp
                  trfow(jj,jk) =  trwdta(ij,jk,1,m) * twmsk(jj,jk)
            END DO
         END DO
       ELSE

!         DO jk = 1, jpkm1
!            DO jj = njw0p1, njw1m1
!               ij = jj -1 + njmpp
!               trfow(jj,jk) =  ( zxy * trwdta(ij,jk,2,m) + &
!                  &           (1.-zxy) * trwdta(ij,jk,1,m) ) * twmsk(jj,jk)
!            END DO
!         END DO
       ENDIF

      IF( (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. lp_obc_west ) THEN

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
          
         ! Temperature and salinity radiation
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
         ! ... the phase velocity is taken as the phase velocity of the tangen-
         ! ... tial velocity (here vn), which have been saved in (v_cxwbnd)
         DO ji = fs_niw0, fs_niw1 ! Vector opt.
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
         ! ... i-phase speed ratio (from averaged of v_cxwbnd)
                  z05cx = (  0.5 * ( v_cxwbnd(jj,jk) + v_cxwbnd(jj-1,jk) ) ) / e1t(ji+1,jj)
                  z05cx = max( z05cx, -1. )
         ! ... z05cx > 0, inflow  zin=0, ztau=1    
         !           < 0, outflow zin=1, ztau=rtauw
                  zin = sign( 1., -1.* z05cx )
                  zin = 0.5*( zin + abs(zin) )
                  ztau = (1.-zin )*rtauwin + zin * rtauw
                  z05cx = z05cx * zin
         ! ... update (ta,sa) with radiative or climatological (t, s)
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1. - twmsk(jj,jk)) +           &
                                 twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
                                 * trwbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
                                 * trwbnd(jj,jk,nibm,nit ) + ztau * trfow (jj,jk) ) &
                                 / (1. - z05cx)
               END DO
            END DO
         END DO

      END IF
!     END IF
!    END DO
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
      USE trcobc_oce      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE lbclnk          ! ???
      USE in_out_manager  ! I/O manager
      USE trc
      USE obctrc
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

   REAL(wp), DIMENSION(:), ALLOCATABLE :: tcobc      ! time_counter variable of BCs


      !! * Local declaration
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   z05cx, ztau, zin
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if ltrfbcnorth is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------
! Linear interpolation of BCs to current time step
!      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
!         zxy = 0.
!ATTENZIONE PER ORA OK MA DA CAMBIARE QUANDO SI AVRANNO I DATI SULLE OBC
!      IF( ntobc == 12 .OR. ntrcobc_dta == 1)   THEN
!        zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!    ELSE
!         zxy = (tcobc(ntobc1)-FLOAT(isrel))/(tcobc(ntobc1)-tcobc(ntobc2))
!      ENDIF

! Per evitare la duplicazione degli array trndta
!       DO m = 1,NO_D3_BOX_STATES
!         IF (D3STATEOBC(m)/= 0.) THEN

      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
         DO jk = 1, jpkm1
            DO ji = nin0p1, nin1m1
               ii = ji -1 + nimpp
                  trfon(ji,jk) = trndta(ii,jk,1,m) * tnmsk(ji,jk)
            END DO
         END DO
    
         ELSE 
     ENDIF

      IF( (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. lp_obc_north ) THEN

         DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                   
              IF(ji.eq.2) tra(ji,jj,jk,1)=tra(ji+1,jj,jk,1)
              IF(ji.ge.4.and.ji.le.6) tra(ji,jj,jk,1)=tra(ji-1,jj,jk,1)

                  tra(ji,jj,jk,1)= tra(ji,jj,jk,1) * (1.-tnmsk(ji,jk)) + &
                                tnmsk(ji,jk) * trfon(ji,jk)
               END DO
            END DO
         END DO

      ELSE

      ! 2. Beyond the fourth time step if ltrfbcnorth is .FALSE.
      ! -------------------------------------------------------
          
         ! Temperature and salinity radiation
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
         ! ... the phase velocity is taken as the normal phase velocity of the tangen-
         ! ... tial velocity (here un), which has been saved in (u_cynbnd)
         ! ... jpjnob+1,(jpindp1, jpinfm1)
         DO jj = fs_njn0+1, fs_njn1+1 ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
         ! ... j-phase speed ratio (from averaged of vtnbnd)
         !        (bounded by 1)
                  z05cx = ( 0.5 * ( u_cynbnd(ji,jk) + u_cynbnd(ji-1,jk) ) ) / e2t(ji,jj-1)
                  z05cx = min( z05cx, 1. )
         ! ... z05cx=< 0, inflow  zin=0, ztau=1    
         !           > 0, outflow zin=1, ztau=rtaun
                  zin = sign( 1., z05cx )
                  zin = 0.5*( zin + abs(zin) )
         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
                  ztau = (1.-zin ) * rtaunin + zin * rtaun
                  z05cx = z05cx * zin
         ! ... update (ta,sa) with radiative or climatological (t, s)
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1.-tnmsk(ji,jk)) +             &
                                 tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
                                 * trnbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
                                 * trnbnd(ji,jk,nibm,nit ) + ztau * trfon (ji,jk) ) &
                                 / (1. + z05cx)
               END DO
            END DO
         END DO

      END IF
!     END IF
!    ENDDO
   END SUBROUTINE obc_trc_bfm_north


   SUBROUTINE obc_trc_bfm_south ( kt ,m)
      !!------------------------------------------------------------------------------
      !!                ***  SUBROUTINE obc_tra_south  ***
      !!     
      !! ** Purpose :
      !!      Apply the radiation algorithm on south OBC tracers ta, sa using the 
      !!      phase velocities calculated in obc_rad_south subroutine in obcrad.F90 module
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
      USE trcobc_oce      ! ocean open boundary conditions
      USE lib_mpp         ! ???
      USE lbclnk          ! ???
      USE in_out_manager  ! I/O manager
      USE trc
      USE obctrc
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

   REAL(wp), DIMENSION(:), ALLOCATABLE :: tcobc      ! time_counter variable of BCs

      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbcsouth is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------
! Linear interpolation of BCs to current time step
!      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
!         zxy = 0.
!ATTENZIONE PER ORA OK MA DA CAMBIARE QUANDO SI AVRANNO I DATI SULLE OBC
!      IF( ntobc == 12 .OR. ntrcobc_dta == 1)   THEN
!         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!      ELSE
!         zxy = (tcobc(ntobc1)-FLOAT(isrel))/(tcobc(ntobc1)-tcobc(ntobc2))
!      ENDIF

! Per evitare la duplicazione degli array trsdta
!       DO m = 1,NO_D3_BOX_STATES
!         IF (D3STATEOBC(m)/= 0.) THEN

      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN
         DO jk = 1, jpkm1
           DO ji = nis0p1, nis1m1
              ii = ji -1 + nimpp
                 trfos(ji,jk) =   trsdta(ii,jk,1,m) * tsmsk(ji,jk)
           END DO
         END DO
         ELSE
        ENDIF

      IF(  (kt < nit000+3 .AND. .NOT.ln_rstart) .OR. lp_obc_south) THEN

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
         !... the phase velocity is taken as the phase velocity of the tangen-
         !... tial velocity (here un), which has been saved in (u_cysbnd)
         !... jpjsob,(jpisdp1, jpisfm1)
         DO jj = fs_njs0, fs_njs1  ! Vector opt.
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
         !... j-phase speed ratio (from averaged of u_cysbnd)
         !       (bounded by 1)
                  z05cx = ( 0.5 * ( u_cysbnd(ji,jk) + u_cysbnd(ji-1,jk) ) ) / e2t(ji,jj+1)
                  z05cx = max( z05cx, -1. )
         !... z05cx > 0, inflow  zin=0, ztau=1
         !          < 0, outflow zin=1, ztau=rtaus
                  zin = sign( 1., -1.* z05cx )
                  zin = 0.5*( zin + abs(zin) )
                  ztau = (1.-zin ) + zin * rtaus
                  z05cx = z05cx * zin
         !... update (ta,sa) with radiative or climatological (t, s)
                  tra(ji,jj,jk,1) = tra(ji,jj,jk,1) * (1.-tsmsk(ji,jk)) +             &
                                 tsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
                                 * trsbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
                                 * trsbnd(ji,jk,nibm,nit ) + ztau * trfos (ji,jk) ) &
                                 / (1. - z05cx)
               END DO
            END DO
         END DO

      END IF   
   END SUBROUTINE obc_trc_bfm_south

#endif
