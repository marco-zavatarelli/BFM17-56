SUBROUTINE trc_sbc_bfm ( kt, m )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sbc_bfm  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (concentration/dilution effect) and add it to the general 
      !!       trend of tracer equations.
      !!
      !! ** Method :
      !!      * concentration/dilution effect:
      !!            The surface freshwater flux modify the ocean volume
      !!         and thus the concentration of a tracer as :
      !!            tra = tra + emp * trn / e3t   for k=1
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation ) given in kg/m2/s is divided
      !!         by 1000 kg/m3 (density of plain water) to obtain m/s.
      !!
      !!         Runoff is separated into water flux (included in emps)
      !!         and mass flux
      !!
      !! ** Action  : - Update the 1st level of tra with the trend associated
      !!                with the tracer surface boundary condition 
      !!
      !! History :
      !!   8.2  !  98-10  (G. Madec, G. Roullet, M. Imbard)  Original code
      !!   8.2  !  01-02  (D. Ludicone)  sea ice and free surface
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !! Adapted to usage with the BFM by M. Vichi (CMCC-INGV)
   !! Added runoff input of tracer concentrations
   !!----------------------------------------------------------------------

   ! NEMO
   USE oce_trc              ! ocean dynamics and active tracers variables
   USE sbcrnf               ! contains river runoff (kg/m2/s) and river mask
   USE trc                  ! ocean  passive tracers variables
   USE fldread

   ! BFM
   use api_bfm
   use mem
   use global_mem,only:LOGUNIT
   use sbc_oce, only: ln_rnf
 
   ! substitutions
#  include "top_substitute.h90"

   !! * Arguments
   INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
   integer, intent(IN)     ::  m      ! BFM variable index

   !! * Local declarations
   INTEGER  ::   ji, jj, jn           ! dummy loop indices
   REAL(wp) ::   ztra, zsrau, zse3t   ! temporary scalars
   !!----------------------------------------------------------------------


   ! 0. initialization
   zsrau = 1. / rau0
   IF( .NOT. ln_sco )  zse3t = 1. / fse3t(1,1,1)

    ! Concentration and dilution effect on tra
    DO jj = 2, jpj
       DO ji = fs_2, fs_jpim1   ! vector opt.
           IF ( ln_sco ) zse3t = 1. / fse3t(ji,jj,1)
           ! concent./dilut. effect
           ztra = zsrau * zse3t * tmask(ji,jj,1) * emps(ji,jj) * trn(ji,jj,1,1) 
           ! add the trend to the general tracer trend
           tra(ji,jj,1,1) = tra(ji,jj,1,1) + ztra
        END DO
    END DO
	
	IF (ln_rnf) THEN
       ! Add mass from prescribed river concentration
	   allocate(rtmp2d(jpi,jpj))
       rtmp2d(:,:) = unpack(PELRIVER(m,:),SEAmask(:,:,1),ZEROS(:,:,1))
       tra(:,:,1,1) = tra(:,:,1,1) - zsrau*sf_rnf(1)%fnow(:,:,1)*rtmp2d(:,:)/fse3t(:,:,1) 
	   deallocate(rtmp2d)
	END IF

   END SUBROUTINE trc_sbc_bfm

