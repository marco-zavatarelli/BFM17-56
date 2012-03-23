MODULE trcbc
   !!======================================================================
   !!                     ***  MODULE  trcdta  ***
   !! TOP :  module for passive tracer boundary conditions
   !!=====================================================================
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

   USE mem, ONLY: NO_D3_BOX_STATES

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_bc_init    ! called in trcini.F90 

   INTEGER  , SAVE, PUBLIC                             :: nb_trcobc   ! number of tracers with open BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trcsbc   ! number of tracers with surface BC
   INTEGER  , SAVE, PUBLIC                             :: nb_trccbc   ! number of tracers with coastal BC
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indobc ! index of tracer with OBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indsbc ! index of tracer with SBC data
   INTEGER  , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: n_trc_indcbc ! index of tracer with CBC data
   INTEGER  , SAVE, PUBLIC                             :: ntra_obc     ! MAX( 1, nb_trcxxx ) to avoid compilation error with bounds checking
   INTEGER  , SAVE, PUBLIC                             :: ntra_sbc     ! MAX( 1, nb_trcxxx ) to avoid compilation error with bounds checking
   INTEGER  , SAVE, PUBLIC                             :: ntra_cbc     ! MAX( 1, nb_trcxxx ) to avoid compilation error with bounds checking
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trofac   ! multiplicative factor for OBCtracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcobc   ! structure of data input OBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trsfac   ! multiplicative factor for SBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trcsbc   ! structure of data input SBC (file informations, fields read)
   REAL(wp) , SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: rf_trcfac   ! multiplicative factor for CBC tracer values
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: sf_trccbc   ! structure of data input CBC (file informations, fields read)

!MAV: it is also needed to define the arrays that have the same functionalities as ln_trc_ini in the same module (trc)
! ln_trc_obc, ln_trc_sbc, ln_trc_cbc

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

   SUBROUTINE trc_bc_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_bc_init  ***
      !!                    
      !! ** Purpose :   initialisation of passive tracer BC data 
      !! 
      !! ** Method  : - Read namtsd namelist
      !!              - allocates passive tracer BC data structure 
      !!----------------------------------------------------------------------
      !
      INTEGER            :: jl, jn                   ! dummy loop indicies
      INTEGER            :: ierr0, ierr1, ierr2, ierr3       ! temporary integers
      CHARACTER(len=100) :: clndta, clntrc
      !
      CHARACTER(len=100) :: cn_dir
      TYPE(FLD_N), DIMENSION(NO_D3_BOX_STATES) :: slf_i     ! array of namelist informations on the fields to read
      TYPE(FLD_N), DIMENSION(NO_D3_BOX_STATES) :: sn_trcobc ! open
      TYPE(FLD_N), DIMENSION(NO_D3_BOX_STATES) :: sn_trcsbc ! surface
      TYPE(FLD_N), DIMENSION(NO_D3_BOX_STATES) :: sn_trccbc ! coastal
      REAL(wp)   , DIMENSION(NO_D3_BOX_STATES) :: rn_trofac    ! multiplicative factor for tracer values
      REAL(wp)   , DIMENSION(NO_D3_BOX_STATES) :: rn_trsfac    ! multiplicative factor for tracer values
      REAL(wp)   , DIMENSION(NO_D3_BOX_STATES) :: rn_trcfac    ! multiplicative factor for tracer values
      !!
      NAMELIST/namtrc_bc/ cn_dir, sn_trcobc, rn_trofac, sn_trcsbc, rn_trsfac, sn_trccbc, rn_trcfac 
      !!----------------------------------------------------------------------
      !
      !  Initialisation
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0  

      ! Compute the number of tracers to be initialised with open, surface and boundary data
      ALLOCATE( n_trc_indobc(NO_D3_BOX_STATES), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_init: unable to allocate n_trc_indobc' )   ;   RETURN
      ENDIF
      nb_trcobc      = 0
      n_trc_indobc(:) = 0
      !
      ALLOCATE( n_trc_indsbc(NO_D3_BOX_STATES), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_init: unable to allocate n_trc_indsbc' )   ;   RETURN
      ENDIF
      nb_trcsbc      = 0
      n_trc_indsbc(:) = 0
      !
      ALLOCATE( n_trc_indcbc(NO_D3_BOX_STATES), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'trc_bc_init: unable to allocate n_trc_indcbc' )   ;   RETURN
      ENDIF
      nb_trccbc      = 0
      n_trc_indcbc(:) = 0
      !
      DO jn = 1, NO_D3_BOX_STATES
         IF( ln_trc_obc(jn) ) THEN
             nb_trcobc       = nb_trcobc + 1 
             n_trc_indobc(jn) = nb_trcobc 
         ENDIF
         IF( ln_trc_sbc(jn) ) THEN
             nb_trcsbc       = nb_trcsbc + 1
             n_trc_indsbc(jn) = nb_trcsbc
         ENDIF
         IF( ln_trc_cbc(jn) ) THEN
             nb_trccbc       = nb_trccbc + 1
             n_trc_indcbc(jn) = nb_trccbc
         ENDIF
      ENDDO
      ntra_obc = MAX( 1, nb_trcobc )   ! To avoid compilation error with bounds checking
      IF( lwp ) WRITE(numout,*) ' '
      IF( lwp ) WRITE(numout,*) ' number of passive tracers to be initialized with open boundary data :', nb_trcobc
      IF( lwp ) WRITE(numout,*) ' '
      ntra_sbc = MAX( 1, nb_trcsbc )   ! To avoid compilation error with bounds checking
      IF( lwp ) WRITE(numout,*) ' '
      IF( lwp ) WRITE(numout,*) ' number of passive tracers to be initialized with surface boundary data :', nb_trcsbc
      IF( lwp ) WRITE(numout,*) ' '
      ntra_cbc = MAX( 1, nb_trccbc )   ! To avoid compilation error with bounds checking
      IF( lwp ) WRITE(numout,*) ' '
      IF( lwp ) WRITE(numout,*) ' number of passive tracers to be initialized with coastal boundary data :', nb_trccbc
      IF( lwp ) WRITE(numout,*) ' '

      ! Initialize the namelists with default values
      cn_dir  = './'            ! directory in which the model is executed
      DO jn = 1, NO_D3_BOX_STATES
         WRITE( clndta,'("TR_",I1)' ) jn
         clndta = TRIM( clndta )
         !                 !  file      ! frequency ! variable  ! time inter !  clim   ! 'yearly' or ! weights  ! rotation !
         !                 !  name      !  (hours)  !  name     !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs    !
         sn_trcobc(jn)  = FLD_N( clndta ,   -1      , clndta    ,  .false.   , .true.  ,  'monthly'  , ''       , ''       )
         sn_trcsbc(jn)  = FLD_N( clndta ,   -1      , clndta    ,  .false.   , .true.  ,  'monthly'  , ''       , ''       )
         sn_trccbc(jn)  = FLD_N( clndta ,   -1      , clndta    ,  .false.   , .true.  ,  'monthly'  , ''       , ''       )
         rn_trofac(jn) = 1._wp
         rn_trsfac(jn) = 1._wp
         rn_trcfac(jn) = 1._wp
      END DO
      !
      REWIND( numnat )               ! read nattrc
      READ  ( numnat, namtrc_bc )

      ! print some information for each 
      IF( lwp ) THEN
         DO jn = 1, NO_D3_BOX_STATES
            IF( ln_trc_obc(jn) )  THEN    
               clndta = TRIM( sn_trcobc(jn)%clvar ) 
               IF(lwp) WRITE(numout,*) ' preparing to read OBC data file for passive tracer number :', jn, ' name : ', clndta, & 
               &               ' multiplicative factor : ', rn_trofac(jn)
            ENDIF
            IF( ln_trc_sbc(jn) )  THEN    
               clndta = TRIM( sn_trcsbc(jn)%clvar ) 
               IF(lwp) WRITE(numout,*) ' preparing to read SBC data file for passive tracer number :', jn, ' name : ', clndta, & 
               &               ' multiplicative factor : ', rn_trsfac(jn)
            ENDIF
            IF( ln_trc_cbc(jn) )  THEN    
               clndta = TRIM( sn_trccbc(jn)%clvar ) 
               IF(lwp) WRITE(numout,*) ' preparing to read CBC data file for passive tracer number :', jn, ' name : ', clndta, & 
               &               ' multiplicative factor : ', rn_trcfac(jn)
            ENDIF
         END DO
      ENDIF
      !
      ! The following code is written this way to reduce memory usage and repeated for each boundary data
!MAV: note that this is just a placeholder and the dimensions must be changed according to 
!     what will be done with BDY. A new structure will probably need to be included
      IF( nb_trcobc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trcobc(nb_trcobc), rf_trofac(nb_trcobc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_init: unable to allocate  sf_trcdta structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, NO_D3_BOX_STATES
            IF( ln_trc_obc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indobc(jn)
               slf_i(jl)    = sn_trcobc(jn)
               rf_trofac(jl) = rn_trofac(jn)
                                            ALLOCATE( sf_trcobc(jl)%fnow(jpi,jpj,jpk)   , STAT=ierr2 )
               IF( sn_trcobc(jn)%ln_tint )  ALLOCATE( sf_trcobc(jl)%fdta(jpi,jpj,jpk,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc_init : unable to allocate passive tracer OBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         ENDDO
         !                         ! fill sf_trcdta with slf_i and control print
         CALL fld_fill( sf_trcobc, slf_i, cn_dir, 'trc_obc', 'Passive tracer OBC data', 'namtrc_bc' )
         !
      ENDIF
      !
      IF( nb_trcsbc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trcsbc(nb_trcsbc), rf_trsfac(nb_trcsbc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_bc_init: unable to allocate  sf_trcsbc structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, NO_D3_BOX_STATES
            IF( ln_trc_sbc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indsbc(jn)
               slf_i(jl)    = sn_trcsbc(jn)
               rf_trsfac(jl) = rn_trsfac(jn)
                                            ALLOCATE( sf_trcsbc(jl)%fnow(jpi,jpj,1)   , STAT=ierr2 )
               IF( sn_trcsbc(jn)%ln_tint )  ALLOCATE( sf_trcsbc(jl)%fdta(jpi,jpj,1,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc_init : unable to allocate passive tracer SBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         ENDDO
         !                         ! fill sf_trcsbc with slf_i and control print
         CALL fld_fill( sf_trcsbc, slf_i, cn_dir, 'trc_sbc', 'Passive tracer SBC data', 'namtrc_bc' )
         !
      ENDIF
      !
      IF( nb_trccbc > 0 ) THEN       !  allocate only if the number of tracer to initialise is greater than zero
         ALLOCATE( sf_trccbc(nb_trccbc), rf_trcfac(nb_trccbc), STAT=ierr1 )
         IF( ierr1 > 0 ) THEN
            CALL ctl_stop( 'trc_dta_ini: unable to allocate  sf_trccbc structure' )   ;   RETURN
         ENDIF
         !
         DO jn = 1, NO_D3_BOX_STATES
            IF( ln_trc_cbc(jn) ) THEN      ! update passive tracers arrays with input data read from file
               jl = n_trc_indcbc(jn)
               slf_i(jl)    = sn_trccbc(jn)
               rf_trcfac(jl) = rn_trcfac(jn)
                                            ALLOCATE( sf_trccbc(jl)%fnow(jpi,jpj,1)   , STAT=ierr2 )
               IF( sn_trccbc(jn)%ln_tint )  ALLOCATE( sf_trccbc(jl)%fdta(jpi,jpj,1,2) , STAT=ierr3 )
               IF( ierr2 + ierr3 > 0 ) THEN
                 CALL ctl_stop( 'trc_bc : unable to allocate passive tracer CBC data arrays' )   ;   RETURN
               ENDIF
            ENDIF
            !   
         ENDDO
         !                         ! fill sf_trccbc with slf_i and control print
         CALL fld_fill( sf_trccbc, slf_i, cn_dir, 'trccbc', 'Passive tracer CBC data', 'namtrc_bc' )
         !
      ENDIF

   END SUBROUTINE trc_bc_init
#endif

   !!======================================================================
END MODULE trcbc
