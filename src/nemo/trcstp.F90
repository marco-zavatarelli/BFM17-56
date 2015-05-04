MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE api_bfm,    ONLY: bio_calc
   USe global_mem, ONLY: LOGUNIT
   USE iom
   USE in_out_manager
   USE trcsub
   USE trcbc,      ONLY: trc_bc_read

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step
   
   !!----------------------------------------------------------------------
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !! Additional software compliant with GPL
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_stp( kt )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      CHARACTER (len=25)    ::  charout
      !!-------------------------------------------------------------------
#ifdef DEBUG
      write(LOGUNIT,*) 'Timestep: ',kt, nit000, nn_dttrc
#endif

      IF( nn_timing == 1 )   CALL timing_start('trc_stp')

      !---------------------------------------------
      ! Check the main BFM flag
      !---------------------------------------------
      IF (bio_calc) THEN

          IF( (nn_dttrc /= 1 ) .AND. (kt == nit000) ) THEN
             ! this is now done in trcini.F90. may not be necessary here
             !CALL trc_sub_ini                              ! Initialize variables for substepping passive tracers
             CALL trc_bc_read( kt )                        ! Read initial Boundary Conditions
          ENDIF

          IF ( nn_dttrc /= 1 )     CALL trc_sub_stp( kt )                    ! Averaging physical variables for sub-stepping

          !---------------------------------------------
          ! Proceed only every nn_dttrc
          !---------------------------------------------   
          IF ( MOD( kt , nn_dttrc ) == 0 ) THEN  
 
                                   CALL trc_bc_read ( kt )      ! read/update Boundary Conditions
                                   CALL trc_bfm( kt )           ! main call to BFM
                                   CALL trc_trp_bfm( kt )       ! transport of BFM tracers
                                   CALL trc_dia_bfm( kt )       ! diagnostic output for BFM
             IF ( nn_dttrc /= 1 )  CALL trc_sub_reset( kt )     ! reset physical variables after sub-stepping

          ENDIF 

      END IF

      IF( nn_timing == 1 )   CALL timing_stop('trc_stp')
      FLUSH(LOGUNIT)

   END SUBROUTINE trc_stp

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif

   !!======================================================================
END MODULE trcstp
