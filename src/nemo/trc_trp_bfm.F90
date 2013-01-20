#include"cppdefs.h"
SUBROUTINE trc_trp_bfm( kt )
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

   !! * BFM Modules used
   USE constants, ONLY: SEC_PER_DAY
   USE global_mem
   USE mem,       ONLY: NO_D3_BOX_STATES,D3STATETYPE, &
                        D3SOURCE,D3STATE,NO_BOXES
#ifndef ONESOURCE
   USE mem,       ONLY: D3SINK
#endif
   use mem, only: ppO3c
   USE mem_param, ONLY: CalcTransportFlag, CalcConservationFlag, p_small
   USE api_bfm
   USE trcnxtbfm       ! time-stepping                       (trc_nxt_bfm routine)
   USE print_functions

   IMPLICIT NONE

   !! * Substitutions (fse3t for partial-steps and z-coords)
#  include "domzgr_substitute.h90"

   !! * Arguments
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
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
   !-----------------------------------------------------------------------

#if defined key_obcbfm
      CALL trcobc_dta_bfm( kt )    ! OBC for BFM
#endif

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
#if defined D1SOURCE
            dummy(:) = D3SOURCE(m,:)
#else
            do k=1,NO_D3_BOX_STATES
               do n=1,NO_BOXES
#ifdef ONESOURCE
                  dummy(n) = dummy(n) + D3SOURCE(m,k,n)
#else
                  dummy(n) = dummy(n) + D3SOURCE(m,k,n) - D3SINK(m,k,n)
#endif
               end do
            end do
#endif

            ! Add biogeochemical trends
            tra(:,:,:,1) = unpack(dummy,SEAmask,ZEROS)
#ifdef DEBUG
            LEVEL2 'trc_trp_bfm at',kt
            CALL prxy( LOGUNIT, 'trn for '//trim(var_names(m)),trn(:,:,1,1), jpi, 1, jpj, 1, ZERO)
            CALL prxy( LOGUNIT, 'tra:BIO',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif
            IF (.NOT.CalcConservationFlag) THEN
               ! NOTE: these routines do not conserve mass,
               ! because non-dynamical volume is used;
               ! thus, excluded for mass conservation checkings)
                CALL trc_sbc_bfm( kt,m )   ! surface boundary condition including rivers
            END IF
            CALL trc_set_bfm( kt, m)      ! set other boundary conditions and compute sinking

            IF( lk_trabbl )        CALL trc_bbl( kt )            ! advective (and/or diffusive) bottom boundary layer scheme
!MAV: still no defined for BFM
!            IF( lk_trcdmp     )   CALL trc_dmp( kt )            ! internal damping trends
             CALL trc_adv( kt )                                  ! horizontal & vertical advection

#ifdef DEBUG
            CALL prxy( LOGUNIT, 'tra:ADV',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif

            IF( ln_zps .AND. .NOT. lk_c1d ) &
               &     CALL zps_hde( kt, jptra, trb, gtru, gtrv )  ! Partial steps: now horizontal gradient
                                                                 ! gtru and gtrv are computed for each tracer
                                                                 ! therefore it is moved here
            CALL trc_ldf( kt )                                   ! lateral mixing
#ifdef DEBUG
            CALL prxy( LOGUNIT, 'tra:LDF',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif
            IF( .NOT. lk_offline .AND. lk_zdfkpp )    &
               &                 CALL trc_kpp( kt )              ! KPP non-local tracer fluxes
            CALL trc_zdf( kt )                                   ! vertical mixing and after tracer fields
#ifdef DEBUG
            CALL prxy( LOGUNIT, 'trn:ZDF',trn(:,:,1,1), jpi, 1, jpj, 1, ZERO)
            CALL prxy( LOGUNIT, 'tra:ZDF',tra(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif

            IF (ln_top_euler) THEN
               CALL trc_nxt_bfm( kt, m )                         ! tracer fields at next time step
            ELSE
               CALL trc_nxt( kt )                                ! tracer fields at next time step
            END IF

            ! Remap the biochemical variables from 3D
            ! to 1D (apply land-sea mask)
            ! Values have been updated in trcnxt
            D3STATE(m,:) = pack(trn(:,:,:,1),SEAmask)
            IF (ln_top_euler) THEN
               D3STATEB(m,:) = D3STATE(m,:)
            ELSE
               D3STATEB(m,:) = pack(trb(:,:,:,1),SEAmask)        ! Leap-frog scheme 
            END IF
#ifdef DEBUG
            CALL prxy( LOGUNIT, 'trn:NXT',trn(:,:,1,1), jpi, 1, jpj, 1, ZERO)
#endif
         END IF ! transported

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
         IF(lwp) WRITE(LOGUNIT,9000) m, trim(var_names(stPelStateS+m-1)), zmean, zmin, zmax, zdrift

      END DO ! over state vars
      
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
               LEVEL2 'Negative concentration at ',kt,'for ',trim(var_names(m))
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

