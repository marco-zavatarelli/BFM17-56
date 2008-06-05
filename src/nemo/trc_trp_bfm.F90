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

   USE trctrp_lec      ! passive tracers transport

   USE trcbbl          ! bottom boundary layer               (trc_bbl routine)
   USE trcdmp          ! internal damping                    (trc_dmp routine)

   USE trcldf_bilapg   ! lateral mixing               (trc_ldf_bilapg routine)
   USE trcldf_bilap    ! lateral mixing                (trc_ldf_bilap routine)
   USE trcldf_iso      ! lateral mixing                  (trc_ldf_iso routine)
   USE trcldf_iso_zps  ! lateral mixing              (trc_ldf_iso_zps routine)
   USE trcldf_lap      ! lateral mixing                  (trc_ldf_lap routine)

   USE trcnxt          ! time-stepping                       (trc_nxt routine)
   USE trcrad          ! positivity                          (trc_rad routine)

   USE trcadv_cen2     ! 2nd order centered advection   (trc_adv_cen2 routine)
   USE trcadv_muscl    ! MUSCL advection               (trc_adv_muscl routine)
   USE trcadv_muscl2   ! MUSCL2 advection             (trc_adv_muscl2 routine)
   USE trcadv_tvd      ! TVD advection                   (trc_adv_tvd routine)
   USE trcadv_smolar   ! SMOLAR advection             (trc_adv_smolar routine)

   USE trczdf_exp      ! vertical diffusion              (trc_zdf_exp routine)
   USE trczdf_imp      ! vertical diffusion              (trc_zdf_exp routine)
   USE trczdf_iso      ! vertical diffusion              (trc_zdf_exp routine)
   USE trczdf_iso_vopt ! vertical diffusion              (trc_zdf_exp routine)
   USE trcsbc          ! surface boundary condition          (trc_sbc routine)

   USE zpshde_trc      ! partial step: hor. derivative   (zps_hde_trc routine)

   !! * BFM Modules used
   USE constants, ONLY: SEC_PER_DAY
   USE global_mem
   USE mem,       ONLY: NO_D3_BOX_STATES,D3STATETYPE, &
                        D3SOURCE,D3STATE,D3SINK,NO_BOXES
   use mem, only: ppO3c
   USE mem_param, ONLY: CalcTransportFlag, CalcConservationFlag
   USE api_bfm

   IMPLICIT NONE

   !! * Substitutions (fse3t for partial-steps and z-coords)
#  include "domzgr_substitute.h90"

   !! * Arguments
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
   !! ---------------------------------------------------------------------
      integer :: m,k
      real(RLEN) :: dummy(NO_BOXES)

   !
   ! Exit if transport is not computed. Time integration is carried out 
   ! with an ODE solver in trcbfm.F90
   ! The same is done for benthic variables if active 
   !
      if (.NOT.CalcTransportFlag) return

   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! BFM tracers, loop over number of state variables
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !-----------------------------------------------------------------------

      DO m = 1,NO_D3_BOX_STATES
         IF (D3STATETYPE(m)>=ALLTRANSPORT) THEN
            ! remap the biological states and trends to 3D arrays
            IF ( l_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) THEN
               ! Leap-frog scheme (only in explicit case)
               trb(:,:,:,1) = unpack(D3STATEB(m,:),SEAmask,ZEROS)
               trn(:,:,:,1) = unpack(D3STATE(m,:),SEAmask,ZEROS)
            ELSE
               trb(:,:,:,1) = unpack(D3STATE(m,:),SEAmask,ZEROS)
               trn = trb
            END IF
            ! sum all the rates (loop is faster than intrinsic sum)
            dummy = ZERO
            do k=1,NO_D3_BOX_STATES
               dummy = dummy + D3SOURCE(m,k,:)-D3SINK(m,k,:)
            end do
            tra(:,:,:,1) = unpack(dummy,SEAmask,ZEROS)

            IF (.NOT.CalcConservationFlag) &
               CALL trc_sbc( kt )         ! surface boundary condition (only dilution here, 
                                          ! NOTE: does not conserve mass,
                                          ! because non-dynamical volume is
                                          ! used; thus, excluded for mass conservation
                                          ! checkings)
            CALL trc_set_bfm( kt, m)      ! set other boundary conditions and compute sinking

# if defined key_trcbbc
       CALL ctl_stop( '  Bottom heat flux not yet implemented with passive tracer         ' &
           &          '  Check in trc_trp_bfm routine ' )
# endif
            !                                                      ! bottom boundary condition
            IF( lk_trcbbl_dif    )   CALL trc_bbl_dif( kt )                ! diffusive bottom boundary layer scheme
            IF( lk_trcbbl_adv    )   CALL trc_bbl_adv( kt )                ! advective (and/or diffusive) bottom boundary layer scheme

!MAV: still no defined for BFM
!            IF( lk_trcdmp        )   CALL trc_dmp( kt )            ! internal damping trends

            !                                                      ! horizontal & vertical advection
            IF( ln_trcadv_cen2   )   CALL trc_adv_cen2  ( kt )             ! 2nd order centered scheme
            IF( ln_trcadv_muscl  )   CALL trc_adv_muscl ( kt )             ! MUSCL scheme
            IF( ln_trcadv_muscl2 )   CALL trc_adv_muscl2( kt )             ! MUSCL2 scheme
            IF( ln_trcadv_tvd    )   CALL trc_adv_tvd   ( kt )             ! TVD scheme
            IF( ln_trcadv_smolar )   CALL trc_adv_smolar( kt )             ! SMOLARKIEWICZ scheme


            IF( n_cla == 1   )   THEN
               WRITE(ctmp1,*) ' Cross Land Advection not yet implemented with passive tracer n_cla = ',n_cla
               CALL ctl_stop(ctmp1)
            ENDIF

            !                                                      ! lateral mixing
            IF( l_trcldf_bilapg  )   CALL trc_ldf_bilapg ( kt )            ! s-coord. horizontal bilaplacian
            IF( l_trcldf_bilap   )   CALL trc_ldf_bilap  ( kt )            ! iso-level bilaplacian
            IF( l_trcldf_iso     )   CALL trc_ldf_iso    ( kt )            ! iso-neutral laplacian
            IF( l_trcldf_iso_zps )   CALL trc_ldf_iso_zps( kt )            ! partial step iso-neutral laplacian
            IF( l_trcldf_lap     )   CALL trc_ldf_lap    ( kt )            ! iso-level laplacian

            !                                                      ! vertical diffusion
            IF( l_trczdf_exp     )   CALL trc_zdf_exp( kt )                ! explicit time stepping (time splitting scheme)
            IF( l_trczdf_imp     )   CALL trc_zdf_imp( kt )                ! implicit time stepping (euler backward)
            IF( l_trczdf_iso     )   CALL trc_zdf_iso( kt )                ! isopycnal
            IF( l_trczdf_iso_vo  )   CALL trc_zdf_iso_vopt( kt )           ! vector opt. isopycnal

            CALL trc_nxt( kt )            ! tracer fields at next time step

!            CALL trc_rad( kt )            ! Correct artificial negative concentrations for isopycnal scheme            
                                           ! Nothing is done for the BFM yet

            IF( ln_zps .AND. .NOT. lk_trccfg_1d ) &
               &                     CALL zps_hde_trc( kt, trb, gtru, gtrv )  ! Partial steps: now horizontal gradient
            !                                                                 ! of passive tracers at the bottom ocean level

            ! Remap the biochemical variables from 3D
            ! to 1D (apply land-sea mask)
            ! Values have been updated in trcnxt
            D3STATE(m,:) = pack(trn(:,:,:,1),SEAmask)
            IF ( l_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) &
               D3STATEB(m,:) = pack(trb(:,:,:,1),SEAmask)        ! Leap-frog scheme (only in explicit case)
         END IF ! transported
      END DO ! over state vars

END SUBROUTINE trc_trp_bfm

