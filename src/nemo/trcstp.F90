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
   USE api_bfm, ONLY: bio_calc
#ifdef key_obcbfm
   USE trcobcdta_bfm
#endif
   USE iom
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step
   
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcstp.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
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

      !---------------------------------------------
      ! Check the main BFM flag
      !---------------------------------------------
      IF (bio_calc) THEN

                             CALL trc_bfm( kt )           ! main call to BFM

                             CALL trc_trp_bfm( kt )       ! transport of BFM tracers

                             CALL trc_dia_bfm( kt )       ! diagnostic output for BFM
      END IF

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
