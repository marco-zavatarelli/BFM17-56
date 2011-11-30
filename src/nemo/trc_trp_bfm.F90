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
   USE trcnxtbfm       ! time-stepping                       (trc_nxt_bfm routine)
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

   IMPLICIT NONE

   !! * Substitutions (fse3t for partial-steps and z-coords)
#  include "domzgr_substitute.h90"

   !! * Arguments
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
   !! ---------------------------------------------------------------------
      integer :: m,k,n
      real(RLEN) :: dummy(NO_BOXES)
      ! local variables for clipping
      integer :: nneg(NO_D3_BOX_STATES)
      real(RLEN) :: total(NO_D3_BOX_STATES),mass

   !
   ! Exit if transport is not computed. Time integration is carried out
   ! with an ODE solver in trcbfm.F90
   ! The same is done for benthic variables if active
   !
      if (.NOT.CalcTransportFlag) return

   !-----------------------------------------------------------------------
   ! Read Open boundary conditions data (only if transport is computed)
   !-----------------------------------------------------------------------

#if defined key_obcbfm
      CALL trcobc_dta_bfm( kt )    ! OBC for BFM
#endif

   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! BFM tracers, loop over number of state variables
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !-----------------------------------------------------------------------
      DO m = 1,NO_D3_BOX_STATES
         IF (D3STATETYPE(m)>=ALLTRANSPORT) THEN
            ! remap the biological states and trends to 3D arrays
#ifdef USEPACK
            IF ( ln_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) THEN
               ! Leap-frog scheme (only in explicit case)
               trb(:,:,:,1) = unpack(D3STATEB(m,:),SEAmask,ZEROS)
               trn(:,:,:,1) = unpack(D3STATE(m,:),SEAmask,ZEROS)
            ELSE
               trb(:,:,:,1) = unpack(D3STATE(m,:),SEAmask,ZEROS)
               trn = trb
            END IF
#else
            IF ( ln_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) THEN
               DO n = 1,NO_BOXES
                  ! Leap-frog scheme (only in explicit case)
                  trb(iwet(n),jwet(n),kwet(n),1) = D3STATEB(m,n)
                  trn(iwet(n),jwet(n),kwet(n),1) = D3STATE(m,n)
               END DO
            ELSE
               DO n = 1,NO_BOXES
                  trb(iwet(n),jwet(n),kwet(n),1) = D3STATE(m,n)
                  trn = trb
               END DO
            END IF
#endif

#if defined D1SOURCE
            dummy(:) = D3SOURCE(m,:)
#else
            ! sum all the rates (loop is faster than intrinsic sum)
            dummy(:) = ZERO
            do k=1,NO_D3_BOX_STATES
               do n=1,NO_BOXES
#if defined ONESOURCE
                  dummy(n) = dummy(n) + D3SOURCE(m,k,n)
#else
                  dummy(n) = dummy(n) + D3SOURCE(m,k,n) - D3SINK(m,k,n)
#endif
               end do
            end do
#endif

#ifdef USEPACK
            tra(:,:,:,1) = unpack(dummy,SEAmask,ZEROS)
#else
            DO n = 1,NO_BOXES
               tra(iwet(n),jwet(n),kwet(n),1) = dummy(n)
            END DO
#endif

            IF (.NOT.CalcConservationFlag) THEN
               ! NOTE: these routines do not conserve mass,
               ! because non-dynamical volume is used;
               ! thus, excluded for mass conservation checkings)
                CALL trc_sbc_bfm( kt,m )   ! surface boundary condition including rivers
            END IF
            CALL trc_set_bfm( kt, m)      ! set other boundary conditions and compute sinking

# if defined key_trcbbc
       CALL ctl_stop( '  Bottom heat flux not yet implemented with passive tracer         ' &
           &          '  Check in trc_trp_bfm routine ' )
# endif
            IF( lk_trabbl )        CALL trc_bbl( kt )            ! advective (and/or diffusive) bottom boundary layer scheme
            !                                                      ! bottom boundary condition
!MAV: still no defined for BFM
!            IF( lk_trcdmp     )   CALL trc_dmp( kt )            ! internal damping trends
             CALL trc_adv( kt )                                  ! horizontal & vertical advection

            IF( nn_cla == 1   )   THEN
               WRITE(ctmp1,*) ' Cross Land Advection not yet implemented with passive tracer nn_cla = ',nn_cla
               CALL ctl_stop(ctmp1)
            ENDIF

            IF( ln_zps .AND. .NOT. lk_trc_c1d ) &
               &                     CALL zps_hde( kt, jptra, trb, gtru, gtrv )  ! Partial steps: now horizontal gradient
                                                                              ! gtru and gtrv are computed for each tracer
            CALL trc_ldf( kt )                                   ! lateral mixing

            CALL trc_zdf( kt )                                   ! vertical mixing and after tracer fields

            CALL trc_nxt_bfm( kt,m )            ! tracer fields at next time step

            !CALL trc_rad_bfm( kt )        ! Correct artificial negative concentrations for isopycnal scheme
            !                                                                 ! of passive tracers at the bottom ocean level
            ! Remap the biochemical variables from 3D
            ! to 1D (apply land-sea mask)
            ! Values have been updated in trcnxt
#ifdef USEPACK
            D3STATE(m,:) = pack(trn(:,:,:,1),SEAmask)
            IF ( ln_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) &
               D3STATEB(m,:) = pack(trb(:,:,:,1),SEAmask)        ! Leap-frog scheme (only in explicit case)
#else
            DO n = 1,NO_BOXES
               D3STATE(m,n) = trn(iwet(n),jwet(n),kwet(n),1)
            END DO
            IF ( ln_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) then
              DO n = 1,NO_BOXES
                 D3STATEB(m,n) = trb(iwet(n),jwet(n),kwet(n),1)
              END DO
            END IF
#endif
         END IF ! transported
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

END SUBROUTINE trc_trp_bfm

