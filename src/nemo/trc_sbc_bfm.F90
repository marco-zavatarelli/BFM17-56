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
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !! Adapted to usage with the BFM by M. Vichi (CMCC-INGV)
   !! Added runoff input of tracer concentrations
   !! Added surface input of tracer concentrations
   !!----------------------------------------------------------------------

   ! NEMO
   USE oce_trc              ! ocean dynamics and active tracers variables
   USE sbcrnf               ! contains river runoff (kg/m2/s) and river mask
   USE trc                  ! ocean  passive tracers variables
   USE fldread
   USE trcbc

   ! BFM
   use api_bfm
   use mem
   use global_mem,only:LOGUNIT
   use sbc_oce, only: ln_rnf
   use constants, only: SEC_PER_DAY
 
   ! substitutions
#  include "top_substitute.h90"

   !! * Arguments
   INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
   integer, intent(IN)     ::  m      ! BFM variable index

   !! * Local declarations
   INTEGER  ::   ji, jj, jn             ! dummy loop indices
   REAL(wp) ::   ztra, zsrau, zse3t     ! temporary scalars
   TYPE(FLD), DIMENSION(1) ::   sf_dta  ! temporary array of information on the field to read
   !!----------------------------------------------------------------------


   ! initialization of density and scale factor
   zsrau = 1. / rau0
   IF( .NOT. ln_sco )  zse3t = 1. / fse3t(1,1,1)

    ! read and add surface input flux if needed
    IF (ln_trc_sbc(m)) THEN
       if (lwp) write(numout,*) 'BFM reading SBC data for variable:',m
       jn = n_trc_indsbc(m)
       sf_dta(1) = sf_trcsbc(jn)
       CALL fld_read( kt, 1, sf_dta )
       ! return the info (needed because fld_read is stupid!)
       sf_trcsbc(jn) = sf_dta(1) 
       DO jj = 2, jpj
          DO ji = fs_2, fs_jpim1   ! vector opt.
             IF ( ln_sco ) zse3t = 1. / fse3t(ji,jj,1)
             ! MAV: units in input files are assumed to be 1/day
             tra(ji,jj,1,1) = tra(ji,jj,1,1) + rf_trsfac(jn) * sf_trcsbc(jn)%fnow(ji,jj,1) &
                              * zse3t / SEC_PER_DAY 
          END DO
       END DO
    END IF

    ! Add mass from prescribed river concentration
    ! MAV: needs to be checked as for surface boundary conditions
    IF (ln_rnf .AND. ln_trc_cbc(m)) THEN
       jn = n_trc_indcbc(m)
       sf_dta = sf_trccbc(jn)
       CALL fld_read( kt, 1, sf_dta )
       ! return the info (needed because fld_read is stupid!)
       sf_trccbc(jn) = sf_dta(1) 
       DO jj = 2, jpj
          DO ji = fs_2, fs_jpim1   ! vector opt.
             IF ( ln_sco ) zse3t = 1. / fse3t(ji,jj,1)
             tra(ji,jj,1,1) = tra(ji,jj,1,1) - zsrau * sf_rnf(1)%fnow(ji,jj,1) &
                              * sf_trccbc(jn)%fnow(ji,jj,1) * zse3t
          END DO
       END DO
    END IF

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
       

if (lwp) write(numout,*) 'exit trc_sbc_bfm',m

   END SUBROUTINE trc_sbc_bfm

