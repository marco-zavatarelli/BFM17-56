MODULE trcdta
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! TOP :  reads passive tracer data 
   !!=====================================================================
   !! History :   1.0  !  2002-04  (O. Aumont)  original code
   !!              -   !  2004-03  (C. Ethe)  module
   !!              -   !  2005-03  (O. Aumont, A. El Moussaoui) F90
   !!            3.4   !  2010-11  (C. Ethe, G. Madec)  use of fldread + dynamical allocation 
   !!----------------------------------------------------------------------
#if  defined key_top 
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP model 
   !!----------------------------------------------------------------------
   !!   trc_dta    : read and time interpolated passive tracer data
   !!----------------------------------------------------------------------
   USE par_trc       !  passive tracers parameters
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables
   USE iom           !  I/O manager
   USE lib_mpp       !  MPP library
   USE fldread       !  read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_dta         ! called in trcini.F90 and trcdmp.F90
   PUBLIC   trc_dta_init    ! called in trcini.F90 

   INTEGER  , PARAMETER, PUBLIC                        :: MAXTRC=100  ! maximum number of tracers 
   INTEGER  , SAVE, PUBLIC                             :: nb_trcdta   ! number of tracers to be initialised with data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_index ! indice of tracer which is initialised with data
   INTEGER  , SAVE, PUBLIC                             :: ntra        ! MAX( 1, nb_trcdta ) to avoid compilation error with bounds checking
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trfac    ! multiplicative factor for tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcdta   ! structure of input SST (file informations, fields read)

! MAV: to be removed. Just for compilation and waiting for the upgrade
LOGICAL , PUBLIC, PARAMETER :: lk_dtatrc = .FALSE. !: temperature data flag

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: trcdta.F90 2977 2011-10-22 13:46:41Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_dta_init(ntrc)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dta_init  ***
      !!                    
      !! ** Purpose :   initialisation of passive tracer input data 
      !! 
      !! ** Method  : - Read namtsd namelist
      !!              - allocates passive tracer data structure 
      !!----------------------------------------------------------------------
      !
      INTEGER,INTENT(IN) :: ntrc
      INTEGER            :: jl, jn                   ! dummy loop indicies
      INTEGER            :: ierr0, ierr1, ierr2, ierr3       ! temporary integers
      CHARACTER(len=100) :: clndta, clntrc
      REAL(wp)           :: zfact
      !
      CHARACTER(len=100) :: cn_dir
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) :: slf_i     ! array of namelist informations on the fields to read
      TYPE(FLD_N), DIMENSION(MAXTRC) :: sn_trcdta
      REAL(wp)   , DIMENSION(MAXTRC) :: rn_trfac    ! multiplicative factor for tracer values
      !!
      NAMELIST/namtrc_dta/ sn_trcdta, cn_dir, rn_trfac 
      !!----------------------------------------------------------------------
      !
      !  Initialisation
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0  
      ! Compute the number of tracers to be initialised with data
      ALLOCATE( n_trc_index(ntrc), slf_i(ntrc), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_nam: unable to allocate n_trc_index' )   ;   RETURN
      ENDIF
      nb_trcdta      = 0
      n_trc_index(:) = 0
      DO jn = 1, ntrc
         IF( ln_trc_ini(jn) ) THEN
             nb_trcdta       = nb_trcdta + 1 
             n_trc_index(jn) = nb_trcdta 
         ENDIF
      ENDDO
      !
      ntra = MAX( 1, nb_trcdta )   ! To avoid compilation error with bounds checking
      IF( lwp ) WRITE(numout,*) ' '
      IF( lwp ) WRITE(numout,*) ' number of passive tracers to be initialized by data :', nb_trcdta
      IF( lwp ) WRITE(numout,*) ' '
      !                         ! allocate the arrays (if necessary)
      !
      cn_dir  = './'            ! directory in which the model is executed
      DO jn = 1, ntrc
         WRITE( clndta,'("TR_",I1)' ) jn
         clndta = TRIM( clndta )
         !                 !  file      ! frequency ! variable  ! time intep !  clim   ! 'yearly' or ! weights  ! rotation !
         !                 !  name      !  (hours)  !  name     !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs    !
         sn_trcdta(jn)  = FLD_N( clndta ,   -1      , clndta    ,  .false.   , .true.  ,  'monthly'  , ''       , ''       )
         !
         rn_trfac(jn) = 1._wp
      END DO
      !
      REWIND( numnat )               ! read nattrc
      READ  ( numnat, namtrc_dta )

      IF( lwp ) THEN
         DO jn = 1, ntrc
            IF( ln_trc_ini(jn) )  THEN    ! open input file only if ln_trc_ini(jn) is true
               clndta = TRIM( sn_trcdta(jn)%clvar ) 
               zfact  = rn_trfac(jn)
               WRITE(numout,*) ' read an initial file for passive tracer number :', jn, ' name : ', clndta, & 
               &               ' multiplicative factor : ', zfact
            ENDIF
         END DO
      ENDIF
      !
      IF( nb_trcdta > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trcdta(nb_trcdta), rf_trfac(nb_trcdta), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_dta_ini: unable to allocate  sf_trcdta structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, ntrc
            IF( ln_trc_ini(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_index(jn)
               slf_i(jl)    = sn_trcdta(jn)
               rf_trfac(jl) = rn_trfac(jn)
                                            ALLOCATE( sf_trcdta(jl)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
               IF( sn_trcdta(jn)%ln_tint )  ALLOCATE( sf_trcdta(jl)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_dta : unable to allocate passive tracer data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         ENDDO
         !                         ! fill sf_trcdta with slf_i and control print
         CALL fld_fill( sf_trcdta, slf_i, cn_dir, 'trc_dta', 'Passive tracer data', 'namtrc' )
         !
      ENDIF
      !
   END SUBROUTINE trc_dta_init


   SUBROUTINE trc_dta( kt, sf_dta, zrf_trfac )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dta  ***
      !!                    
      !! ** Purpose :   provides passive tracer data at kt
      !! 
      !! ** Method  : - call fldread routine
      !!              - s- or mixed z-s coordinate: vertical interpolation on model mesh
      !!              - ln_trcdmp=F: deallocates the data structure as they are not used
      !!
      !! ** Action  :   ptrc   passive tracer data on medl mesh and interpolated at time-step kt
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kt         ! ocean time-step
      TYPE(FLD), DIMENSION(1)   , INTENT(inout) ::   sf_dta     ! array of information on the field to read
      REAL(wp)                  , INTENT(in   ) ::   zrf_trfac  ! multiplication factor
      !
      INTEGER ::   ji, jj, jk, jl, jkk, ik    ! dummy loop indices
      REAL(wp)::   zl, zi
      REAL(wp), DIMENSION(jpk) ::  ztp                ! 1D workspace
      CHARACTER(len=100) :: clndta
      !!----------------------------------------------------------------------
      !
         !
         CALL fld_read( kt, 1, sf_dta )      !==   read data at kt time step   ==!
         !
         !
         IF( ln_sco ) THEN                   !==   s- or mixed s-zps-coordinate   ==!
            !
            IF( kt == nit000 .AND. lwp )THEN
               WRITE(numout,*)
               WRITE(numout,*) 'trc_dta: interpolates passive tracer data onto the s- or mixed s-z-coordinate mesh'
            ENDIF
            !
               DO jj = 1, jpj                         ! vertical interpolation of T & S
                  DO ji = 1, jpi
                     DO jk = 1, jpk                        ! determines the intepolated T-S profiles at each (i,j) points
                        zl = fsdept_0(ji,jj,jk)
                        IF(     zl < gdept_0(1  ) ) THEN          ! above the first level of data
                           ztp(jk) =  sf_dta(1)%fnow(ji,jj,1)
                        ELSEIF( zl > gdept_0(jpk) ) THEN          ! below the last level of data
                           ztp(jk) =  sf_dta(1)%fnow(ji,jj,jpkm1)
                        ELSE                                      ! inbetween : vertical interpolation between jkk & jkk+1
                           DO jkk = 1, jpkm1                                  ! when  gdept(jkk) < zl < gdept(jkk+1)
                              IF( (zl-gdept_0(jkk)) * (zl-gdept_0(jkk+1)) <= 0._wp ) THEN
                                 zi = ( zl - gdept_0(jkk) ) / (gdept_0(jkk+1)-gdept_0(jkk))
                                 ztp(jk) = sf_dta(1)%fnow(ji,jj,jkk) + ( sf_dta(1)%fnow(ji,jj,jkk+1) - &
                                           sf_dta(1)%fnow(ji,jj,jkk) ) * zi 
                              ENDIF
                           END DO
                        ENDIF
                     END DO
                     DO jk = 1, jpkm1
                        sf_dta(1)%fnow(ji,jj,jk) = ztp(jk) * tmask(ji,jj,jk)     ! mask required for mixed zps-s-coord
                     END DO
                     sf_dta(1)%fnow(ji,jj,jpk) = 0._wp
                  END DO
               END DO
            ! 
         ELSE                                !==   z- or zps- coordinate   ==!
            !                             
               sf_dta(1)%fnow(:,:,:) = sf_dta(1)%fnow(:,:,:) * tmask(:,:,:)    ! Mask
               !
               IF( ln_zps ) THEN                      ! zps-coordinate (partial steps) interpolation at the last ocean level
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        ik = mbkt(ji,jj) 
                        IF( ik > 1 ) THEN
                           zl = ( gdept_0(ik) - fsdept_0(ji,jj,ik) ) / ( gdept_0(ik) - gdept_0(ik-1) )
                           sf_dta(1)%fnow(ji,jj,ik) = (1.-zl) * sf_dta(1)%fnow(ji,jj,ik) + zl * sf_dta(1)%fnow(ji,jj,ik-1)
                        ENDIF
                     END DO
                  END DO
               ENDIF
            !
         ENDIF
         !
         sf_dta(1)%fnow(:,:,:) = sf_dta(1)%fnow(:,:,:) * zrf_trfac   !  multiplicative factor
         !
         IF( lwp .AND. kt == nit000 ) THEN
               clndta = TRIM( sf_dta(1)%clvar ) 
               WRITE(numout,*) ''//clndta//' data '
               WRITE(numout,*)
               WRITE(numout,*)'  level = 1'
               CALL prihre( sf_dta(1)%fnow(:,:,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
               WRITE(numout,*)'  level = ', jpk/2
               CALL prihre( sf_dta(1)%fnow(:,:,jpk/2), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
               WRITE(numout,*)'  level = ', jpkm1
               CALL prihre( sf_dta(1)%fnow(:,:,jpkm1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
               WRITE(numout,*)
         ENDIF
      ! 
   END SUBROUTINE trc_dta
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                              NO 3D passive tracer data
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_dta( kt, sf_dta, zrf_trfac )        ! Empty routine
      WRITE(*,*) 'trc_dta: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dta
#endif
   !!======================================================================
END MODULE trcdta
