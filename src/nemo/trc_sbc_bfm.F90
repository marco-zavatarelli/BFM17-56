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
      !!         Runoff is separated
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
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean  passive tracers variables
   USE prtctl_trc          ! Print control for debbuging
   use flx_rnf, only: sriver ! gridpoints with rivers
   ! BFM
   use api_bfm
   use mem
   ! substitutions
#  include "passivetrc_substitute.h90"

   !! * Arguments
   INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
   integer, intent(IN)     ::  m   ! BFM variable index

   !! * Local declarations
   INTEGER  ::   ji, jj, jn           ! dummy loop indices
   REAL(wp) ::   ztra, zsrau, zse3t   ! temporary scalars
   CHARACTER (len=22) :: charout
   !!----------------------------------------------------------------------

   if ( kt == nittrc000 ) .AND. ( m == 1 ) THEN
      if (lwp) WRITE(numout,*)
      if (lwp) WRITE(numout,*) 'trc_set_bfm : BFM tracers surface boundary condition'
      if (lwp) WRITE(numout,*) '            : Initialise river mask and concentrations'
      if (lwp) WRITE(numout,*) '~~~~~~~ '
#ifdef FLUXES
! MAV: to be completed
      !-------------------------------------------------------
      ! Prepares the array containing the 1D mask with
      ! the location of the river grid points
      ! RIVmask is not used yet in this NEMO implementation
      ! but the OPA variables runoff and the new 
      !-------------------------------------------------------
      allocate(RIVmask(NO_BOXES_XY)) 
      allocate(btmp1D(NO_BOXES_XY))
      btmp1D = pack(runoff,SRFmask(:,:,1))
      where (btmp1d>ZERO)
        RIVmask = ONE
      elsewhere
        RIVmask = ZERO
      end where
      deallocate(btmp1D)
      !-------------------------------------------------------
      ! Prepares the array containing the 2D mask 
      !-------------------------------------------------------
      allocate(rmask(jpi,jpj)) 
      where (runoff>ZERO)
        rmask = ONE
      elsewhere
        rmask = ZERO
      end where
#endif
      !-------------------------------------------------------
      ! Fill-in the river concentration
      ! MAV: the strategy is to assign the initial values for 
      ! selected variables
      ! In the future it might be used the initial value close 
      ! to the river
      !-------------------------------------------------------
      allocate(RIVconcentration(NO_BOX_STATES)) 
      RIVconcentration(:) = ZERO
      RIVconcentration(ppO2o) = D3STATE(ppO2o,1)
      RIVconcentration(ppN1p) = D3STATE(ppN1p,1)
      RIVconcentration(ppN3n) = D3STATE(ppN3n,1) 
      RIVconcentration(ppN4n) = D3STATE(ppN4n,1)
      RIVconcentration(ppN5s) = D3STATE(ppN5s,1)
#ifdef INCLUDE_PELCO2
      RIVconcentration(ppO3c) = D3STATE(ppO3c,1)
      RIVconcentration(ppO3h) = D3STATE(ppO3h,1)
#endif
   end if

   ! 0. initialization
   zsrau = 1. / rauw
      IF( .NOT. ln_sco )  zse3t = 1. / fse3t(1,1,1)

      DO jn = 1, jptra
         ! 1. Concentration dilution effect on tra
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1   ! vector opt.
               IF( ln_sco ) zse3t = 1. / fse3t(ji,jj,1)
               ! concent./dilut. effect
               ztra = zsrau * zse3t * tmask(ji,jj,1) *  &
                      ((emps(ji,jj)-runoff(ji,jj)) * trn(ji,jj,1,jn) + &         ! precipitation, evaporation
                       runoff(ji,jj) * (trn(ji,jj,1,jn) - RIVconcentration(m)) ) ! river
               
               ! add the trend to the general tracer trend
               tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ztra
            END DO
         END DO
         
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sbc_bfm')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_sbc_bfm

