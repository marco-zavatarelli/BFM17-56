#include"cppdefs.h"
SUBROUTINE trc_trp_bfm( kstp )
!!----------------------------------------------------------------------
!!                     ***  ROUTINE trc_trp_bfm  ***
!!
!! ** Purpose : Management of passive tracers transport
!!
!! ** Method  :
!!              Compute the passive tracers trends
!!              Update the passive tracers
!!
!! History :
!!   9.0  !  04-03  (C. Ethe)  Original
!!   9.0  !  10-06  (M. Vichi, CMCC-INGV) Adapted to BFM
!!----------------------------------------------------------------------
!! * Modules used
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean passive tracers variables
   USE trcnam_trp      ! passive tracers transport namelist variables
   USE par_trc, ONLY: lk_trc_c1d
   USE trabbl          ! bottom boundary layer               (trc_bbl routine)
   USE trcbbl          ! bottom boundary layer               (trc_bbl routine)
   USE zdfkpp          ! KPP non-local tracer fluxes         (trc_kpp routine)
   USE trcdmp          ! internal damping                    (trc_dmp routine)
   USE trcldf          ! lateral mixing                      (trc_ldf routine)
   USE trcadv          ! advection                           (trc_adv routine)
   USE trczdf          ! vertical diffusion                  (trc_zdf routine)
   USE trcnxt          ! time-stepping                       (trc_nxt routine)
   USE trcrad          ! positivity                          (trc_rad routine)
   USE trcsbc          ! surface boundary condition          (trc_sbc routine)  
   USE zpshde          ! partial step: hor. derivative       (zps_hde routine)

#if defined key_agrif
   USE agrif_top_sponge ! tracers sponges
   USE agrif_top_update ! tracers updates
#endif

   !! * BFM Modules used
   USE constants, ONLY: SEC_PER_DAY
   USE global_mem
   USE mem,       ONLY: NO_D3_BOX_STATES,D3STATETYPE, &
                        D3SOURCE,D3STATE,NO_BOXES
#ifdef EXPLICIT_SINK
   USE mem,       ONLY: D3SINK
#endif
   use mem, only: ppO3c
   USE mem_param, ONLY: CalcTransportFlag, CalcConservationFlag, p_small
   USE api_bfm
   USE trcnxtbfm       ! time-stepping                       (trc_nxt_bfm routine)
   USE print_functions

   IMPLICIT NONE

   !! * Substitutions
#  include "top_substitute.h90"

   !! * Arguments
      INTEGER, INTENT( in ) ::  kstp  ! ocean time-step index
   !! ---------------------------------------------------------------------
      integer :: m,k,n
      real(RLEN) :: dummy(NO_BOXES)
      ! local variables for clipping and statistics
      integer :: nneg(NO_D3_BOX_STATES)
      real(RLEN) :: total(NO_D3_BOX_STATES),mass
      real(RLEN) :: zmean, zmin, zmax, zdrift, ztraf

   !-----------------------------------------------------------------------
   ! Exit if transport is not computed. Time integration is carried out
   ! with an ODE solver in trcbfm.F90
   ! The same is done for benthic variables if active
   !-----------------------------------------------------------------------
      if (.NOT.CalcTransportFlag) return

      if ( nn_timing == 1 )   CALL timing_start('trc_trp_bfm')

   !-----------------------------------------------------------------------
   ! Read Open boundary conditions data (only if transport is computed)
   ! Not yet implemented with NEMO 3.6 BDY
   !-----------------------------------------------------------------------

#ifdef DEBUG
   CALL prxy( LOGUNIT, 'tmask',tmask(:,:,1), jpi, 1, jpj, 1, ZERO)
   CALL prxy( LOGUNIT, 'umask',umask(:,:,1), jpi, 1, jpj, 1, ZERO)
   CALL prxy( LOGUNIT, 'vmask',vmask(:,:,1), jpi, 1, jpj, 1, ZERO)
#endif

   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! BFM tracers, loop over number of state variables
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !-----------------------------------------------------------------------
      DO m = 1,NO_D3_BOX_STATES
         IF (D3STATETYPE(m)>=ALLTRANSPORT) THEN
            ! remap the biological states and trends to 3D arrays
            trn(:,:,:,1) = unpack(D3STATE(m,:),SEAmask,ZEROS)
            IF (ln_top_euler) THEN
               trb(:,:,:,1) = trn(:,:,:,1)
            ELSE
               trb(:,:,:,1) = unpack(D3STATEB(m,:),SEAmask,ZEROS)
            END IF

            dummy(:) = ZERO
            ! sum all the rates
#ifndef EXPLICIT_SINK
            dummy(:) = D3SOURCE(m,:)
#else 
            do k=1,NO_D3_BOX_STATES
               do n=1,NO_BOXES
                  dummy(n) = dummy(n) + D3SOURCE(m,k,n) - D3SINK(m,k,n)
               end do
            end do
#endif
            ! Get biogeochemical trends for every variable
            tra(:,:,:,1) = unpack(dummy,SEAmask,ZEROS)
#ifdef DEBUG
            LEVEL2 'trc_trp_bfm at',kstp
            CALL prxy( LOGUNIT, 'trn for '//trim(var_names(m)),trn(:,:,1,1), jpi, 1, jpj, 1, ZERO)
            CALL prxy( LOGUNIT, 'tra:BIO',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif
            ! Add surface, coastal, and open boundary forcing
            IF (.NOT.CalcConservationFlag) THEN
               ! NOTE: these routines do not conserve mass, because non-dynamical volume may not be used;
               ! thus, excluded for mass conservation checkings
                CALL trc_bc_bfm( kstp, m )
            END IF
            ! Add Vertical Sinking 
            CALL trc_set_bfm( kstp, m)      ! set other boundary conditions and compute sinking

            IF( .NOT. lk_c1d ) THEN ! Compute and integrate physical trends for 3D case
               IF( lk_trabbl )      CALL trc_bbl( kstp )     ! advective (and/or diffusive) bottom boundary layer scheme
               !ST: still no defined for BFM
               !IF( lk_trcdmp )     CALL trc_dmp( kstp )     ! internal damping trends
               !IF( ln_trcdmp_clo ) CALL trc_dmp_clo( kstp ) ! internal damping trends on closed seas only
                                   CALL trc_adv( kstp )     ! horizontal & vertical advection
#ifdef DEBUG
                                   CALL prxy( LOGUNIT, 'tra:ADV',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif
               !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
               ! Partial steps:  now horizontal gradient of passive
               ! tracers at the bottom ocean level  gtru and gtrv
               ! are computed for each BFM tracer here, 
               ! as they are needed in lateral mixing 
               !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
               IF( ln_zps  .AND. .NOT. ln_isfcav)        &
                 &                 CALL zps_hde    ( kstp, jptra, trn, gtru, gtrv )
               IF( ln_zps .AND.        ln_isfcav)        &  ! case with ice shelf cavities
                 &                 CALL zps_hde_isf( kstp, jptra, trn, pgtu=gtru, pgtv=gtrv, pgtui=gtrui, pgtvi=gtrvi )

                                   CALL trc_ldf( kstp )     ! lateral mixing
#ifdef DEBUG
                                   CALL prxy( LOGUNIT, 'tra:LDF',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif

               IF( .NOT. lk_offline .AND. lk_zdfkpp )    &
                 &                 CALL trc_kpp( kstp )     ! KPP non-local tracer fluxes

#if defined key_agrif
               IF(.NOT. Agrif_Root()) CALL Agrif_Sponge_trc ! tracers sponge
#endif

                                   CALL trc_zdf( kstp )     ! vertical mixing and after tracer fields
#ifdef DEBUG
                                   CALL prxy( LOGUNIT, 'trn:ZDF',trn(:,:,1,1), jpi, 1, jpj, 1, ZERO)
                                   CALL prxy( LOGUNIT, 'tra:ZDF',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif

                                   ! Apply lateral boundary conditions (and special open boundary)
                                   CALL trc_nxt_bfm( kstp, m ) 
#ifdef DEBUG
                                   CALL prxy( LOGUNIT, 'trn:NXT',trn(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif
#if defined key_agrif
               IF( .NOT. Agrif_Root())   CALL Agrif_Update_Trc( kstp )   ! Update tracer at AGRIF zoom boundaries : children only
#endif
            ELSE ! Compute and integrate physical trends for 1D case
               IF( .NOT. lk_offline .AND. lk_zdfkpp )    &
                 &                 CALL trc_kpp( kstp )     ! KPP non-local tracer fluxes
                                   CALL trc_zdf( kstp )     ! vertical mixing and after tracer fields
                                   CALL trc_nxt_bfm( kstp, m ) 
            END IF ! 3D and 1D cases

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Remap the biochemical variables from 3D
            ! to 1D (apply land-sea mask)
            ! Values have been updated in trcnxt
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            D3STATE(m,:) = pack(trn(:,:,:,1),SEAmask)
            IF (ln_top_euler) THEN
               D3STATEB(m,:) = D3STATE(m,:)
            ELSE
               D3STATEB(m,:) = pack(trb(:,:,:,1),SEAmask)        ! Leap-frog scheme 
            END IF

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! compute global statistics
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ztraf = glob_sum( trn(:,:,:,1) * cvol(:,:,:) )
            zmin  = MINVAL( trn(:,:,:,1), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
            zmax  = MAXVAL( trn(:,:,:,1), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
            IF( lk_mpp ) THEN
               CALL mpp_min( zmin )      ! min over the global domain
               CALL mpp_max( zmax )      ! max over the global domain
            END IF
            zmean  = ztraf / areatot
            zdrift = ( ( ztraf - D3STATE_tot(m) ) / ( D3STATE_tot(m) + 1.e-12_RLEN )  ) * 100._RLEN
         ELSE ! non-trasported variables
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Time integration of benthic (non-transported) variables  
            ! is done, if active, with an ODE solver in trcbfm.F90
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Compute values from D3STATE to fill statistics
            zmin = MINVAL(D3STATE(m,:))
            zmax = MAXVAL(D3STATE(m,:))
            zmean  = SUM(D3STATE(m,:)) / NO_BOXES
            zdrift = 0.0_RLEN
         END IF ! transported
 
         ! Print statistics into bfm.log file 
         IF ( lwp .AND. ( kstp < 100 .OR. MOD(kstp,50) == 0 ) ) THEN
           IF(m==1) WRITE(LOGUNIT,*) 'Statistics on tracer at step: ' , kstp
           WRITE(LOGUNIT,9000) m, trim(var_names(stPelStateS+m-1)), zmean, zmin, zmax, zdrift
           IF(m==NO_D3_BOX_STATES) WRITE(LOGUNIT,*)
         ENDIF

      END DO ! over BFM state vars

      IF (ln_trcrad) THEN
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Clip negative concentrations (FIX!)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         total(:) = ZERO
         DO m = 1,NO_D3_BOX_STATES
            nneg(m) = 0
            dummy(:) = ZERO
            DO k = 1,NO_BOXES
               if (D3STATE(m,k)<ZERO) THEN
                  dummy(k) = D3STATE(m,k)
                  mass = p_small - D3STATE(m,k) ! store the mass loss
                  total(m) = total(m) + mass
                  D3STATE(m,k) = p_small
                  nneg(m) = nneg(m)+1
               END IF
            END DO
            nneg(m) = min(nneg(m),NO_BOXES-1)
            if (nneg(m)>ZERO) then
               LEVEL2 'Negative concentration at ',kstp,'for ',trim(var_names(m))
               LEVEL3 'number of negative grid points ',nneg(m)
               LEVEL3 'largest negative value ',minval(dummy(:))
               LEVEL3 'total mass added ',total(m),trim(var_units(m))
            end if
         END DO
      END IF

      if ( nn_timing == 1 )   CALL timing_stop('trc_trp_bfm')
9000  FORMAT(' STAT tracer nb :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, &
      &      '    max :',e18.10,'    drift :',e18.10, ' %')

END SUBROUTINE trc_trp_bfm

