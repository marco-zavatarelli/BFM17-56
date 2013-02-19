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
   !!         Runoff is treated separately according to the choice 
   !!         of online/offline coupling and river load
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
   USE iom                  !  print control for debugging
   ! BFM
   use api_bfm
   use mem
   use global_mem,only:LOGUNIT
   use sbc_oce, only: ln_rnf
   use constants, only: SEC_PER_DAY
   use print_functions
 
   ! substitutions
#  include "top_substitute.h90"

   !! * Arguments
   INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
   integer, intent(IN)     ::  m      ! BFM variable index

   !! * Local declarations
   INTEGER  ::   ji, jj, jn             ! dummy loop indices
   REAL(wp) ::   ztra, zsrau, zse3t     ! temporary scalars
   TYPE(FLD), DIMENSION(1) ::   sf_dta  ! temporary array of information on the field to read
   REAL(wp), POINTER, DIMENSION(:,:  ) :: zsfx,field
   CHARACTER (len=25) :: charout
   !!---------------------------------------------------------------------
   !
   IF( nn_timing == 1 )  CALL timing_start('trc_sbc_bfm')
   !
   ! Allocate temporary workspace
   CALL wrk_alloc( jpi, jpj,      zsfx  )
#ifdef DEBUG
   CALL wrk_alloc( jpi, jpj,      field  )
   field = 0.0_wp
#endif
   !!----------------------------------------------------------------------
   ! initialization of density and scale factor
   zsrau = 1._wp / rau0

   IF( kt == nittrc000 ) THEN
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_sbc_bfm : Surface boundary conditions for BFM variables'
      IF(lwp) WRITE(numout,*) '~~~~~~~ '
   ENDIF

   ! Coupling online : 
   ! 1) constant volume: river runoff is added to the horizontal divergence (hdivn) in the subroutine sbc_rnf_div 
   !    one only consider the concentration/dilution effect due to evaporation minus precipitation + 
   !    freezing/melting of sea-ice
   ! 2) variable volume : all freshwater fluxes are added to the volume change (no additional dilution)
   ! Coupling offline : runoff are in emp which contains E-P-R
   ! If there is no input of biogeochemical variables associated to the river, it is assumed that
   ! there is no dilution associated with the river runoff (the freshwater has the same concentration of the sea)
   zsfx(:,:) = emp(:,:)                                    ! - standard case: on/offline coupling with const vol
   IF( .NOT. lk_offline .AND. lk_vvl )  zsfx(:,:) = 0._wp  ! - online coupling with variable volume
   IF( ln_rnf .AND. .NOT. ln_trc_cbc(m) )  &               ! - remove river dilution effect in the absence
      zsfx(:,:) = zsfx(:,:) + rnf(:,:)                     ! of a river load

   ! Concentration and dilution effect on tra 
   DO jj = 2, jpj
      DO ji = fs_2, fs_jpim1   ! vector opt.
           zse3t = 1. / fse3t(ji,jj,1)
           ztra = zsrau * zse3t * tmask(ji,jj,1) * zsfx(ji,jj) * trn(ji,jj,1,1) 
           ! add the trend to the general tracer trend
           tra(ji,jj,1,1) = tra(ji,jj,1,1) + ztra
      END DO
   END DO

    ! read and add surface input flux if needed
    IF (ln_trc_sbc(m)) THEN
       if (lwp) write(numout,*) '   BFM: reading SBC data for variable ', &
                   TRIM( var_names(stPelStateS+m-1) ),' number:',m
       jn = n_trc_indsbc(m)
       sf_dta(1) = sf_trcsbc(jn)
       CALL fld_read( kt, 1, sf_dta )
       ! return the info (needed because fld_read is stupid!)
       sf_trcsbc(jn) = sf_dta(1) 
       DO jj = 2, jpj
          DO ji = fs_2, fs_jpim1   ! vector opt.
             zse3t = 1. / fse3t(ji,jj,1)
             ! The units in BFM input files are 1/day
             tra(ji,jj,1,1) = tra(ji,jj,1,1) + rf_trsfac(jn) * sf_trcsbc(jn)%fnow(ji,jj,1) &
                              * zse3t / SEC_PER_DAY 
          END DO
       END DO
    END IF

    ! Add mass from prescribed river concentration if river values are given
    IF (ln_rnf .AND. ln_trc_cbc(m)) THEN 
          if (lwp) write(numout,*) '   BFM: reading CBC data for variable ', &
                   TRIM( var_names(stPelStateS+m-1) ),' number:',m
          jn = n_trc_indcbc(m)
          sf_dta = sf_trccbc(jn)
          CALL fld_read( kt, 1, sf_dta )
          ! return the info (needed because fld_read is stupid!)
          sf_trccbc(jn) = sf_dta(1) 
          DO jj = 2, jpj
             DO ji = fs_2, fs_jpim1   ! vector opt.
                zse3t = 1. / (e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,1))
                ! Add river loads assuming an istantaneous mixing in the cell volume 
                ! The units in BFM input files are 1/day
                ztra = rn_rfact * rf_trcfac(jn) * sf_trccbc(jn)%fnow(ji,jj,1) * zse3t / SEC_PER_DAY
                tra(ji,jj,1,1) = tra(ji,jj,1,1) + ztra
#ifdef DEBUG
                field(ji,jj) = ztra
#endif
             END DO
          END DO
    END IF

#ifdef DEBUG
    charout = TRIM( var_names(stPelStateS+m-1) )
    WRITE(LOGUNIT,*) ''//charout//' trends in trc_sbc'
    WRITE(LOGUNIT,*)
    WRITE(LOGUNIT,*)'  level = 1'
    CALL prxy( LOGUNIT, 'trend at level = 1',field(:,:), jpi, 1, jpj, 1, ZERO)
    CALL prxy( LOGUNIT, 'trn at level = 1',trn(:,:,1,1), jpi, 1, jpj, 1, ZERO)
    CALL wrk_dealloc( jpi, jpj,      field  )       
#endif

   CALL wrk_dealloc( jpi, jpj,      zsfx  )       
    IF( nn_timing == 1 )  CALL timing_stop('trc_sbc_bfm')

   END SUBROUTINE trc_sbc_bfm

