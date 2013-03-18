#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bfm
!
! !INTERFACE:
   module api_bfm
!
! !DESCRIPTION: 
! API for the BFM. 
! Storage of variables and diagnostics
! To be used in all the coupled applications except
! GOTM, where it actually originated from.
! The GOTM module netcdfout is needed
! Appropriate functions are already available in the GOTM structure

!
! !USE:
   use global_mem, only:RLEN,ZERO,bfm_lwp,LOGUNIT
   use mem,        only:NO_D3_BOX_STATES
   implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bfm, find
!
! !PUBLIC DATA MEMBERS:
   logical                            :: bio_calc,bioshade_feedback,bfm_rstctl
   integer                            :: bio_setup  =1
   integer                            :: bfm_init  =0
   integer                            :: surface_flux_method=-1
   integer                            :: bottom_flux_method=-1
   integer                            :: n_surface_fluxes=-1
   integer                            :: calc_init_bennut_states
   character(len=PATH_MAX)            :: out_dir,out_fname,out_title
   integer                            :: out_units
   integer                            :: out_delta,out_secs
   character(len=PATH_MAX)            :: rst_fname

   !---------------------------------------------
   ! parameters for massive parallel computation
   ! the following are the default values if the 
   ! macro BFM_PARALLEL is not defined
   !---------------------------------------------
   logical                            :: parallel = .FALSE.
   logical                            :: parallel_log = .FALSE.
   integer                            :: parallel_rank = 0
   character(LEN=4)                   :: str

   !---------------------------------------------
   ! Dimension lengths for output
   !---------------------------------------------
   integer        :: lon_len
   integer        :: lat_len
   integer        :: depth_len
   integer        :: ocepoint_len,surfpoint_len,botpoint_len

   !---------------------------------------------
   ! BFM variable information for output
   !---------------------------------------------
   integer, dimension(:), allocatable      :: var_ids
   logical, dimension(:), allocatable      :: var_ave
   real(RLEN)                              :: ave_count
   logical                                 :: ave_ctl = .false.
   real(RLEN),allocatable,dimension(:,:)   :: D3ave
   real(RLEN),allocatable,dimension(:,:)   :: D2ave
   character(len=64), dimension(:), allocatable :: var_names
   character(len=64), dimension(:), allocatable :: var_units
   character(len=64), dimension(:), allocatable :: var_long
   !---------------------------------------------
   ! Indices of the various output variables
   !---------------------------------------------
   integer,public                            :: stPelStateS=0
   integer,public                            :: stPelDiagS=0
   integer,public                            :: stPelFluxS=0
   integer,public                            :: stBenStateS=0
   integer,public                            :: stBenDiagS=0
   integer,public                            :: stBenFluxS=0
   integer,public                            :: stPelStateE=0
   integer,public                            :: stPelDiagE=0
   integer,public                            :: stPelFluxE=0
   integer,public                            :: stBenStateE=0
   integer,public                            :: stBenDiagE=0
   integer,public                            :: stBenFluxE=0

   !---------------------------------------------
   ! Additional output variables
   !---------------------------------------------
   real(RLEN), dimension(:), allocatable   :: c1dim

   !---------------------------------------------
   ! BFM variable information for data input
   ! integer init: select the initialization
   !               0 = homogeneous
   !               1 = analytical
   !               2 = from file
   ! options for init==1
   ! real anv1: value in the surface layer
   ! real anz1: depth of the surface layer
   ! real anv2: value in the bottom layer
   ! real anz2: depth of the bottom layer
   ! options for init==2
   ! char filename: name of the input file
   ! char  varname: name of the var in input file
   ! Options currently used when coupled with NEMO
   ! logical obc: variable has open boundary data
   ! logical sbc: variable has surface boundary data
   ! logical cbc: variable has coastal boundary data
   !---------------------------------------------
   type InputInfo
      integer           :: init
      character(LEN=40) :: filename
      character(LEN=40) :: varname
      real(RLEN)        :: anz1
      real(RLEN)        :: anv1
      real(RLEN)        :: anz2
      real(RLEN)        :: anv2
      logical           :: obc
      logical           :: sbc
      logical           :: cbc
   end type InputInfo
   type(InputInfo),dimension(NO_D3_BOX_STATES) :: InitVar

   !---------------------------------------------
   ! Additional 1D arrays
   !---------------------------------------------
   ! indices of bottom and surface points
   integer,allocatable,dimension(:),public     :: BOTindices,SRFindices
   ! real mask of river points at surface
   real(RLEN),allocatable,dimension(:),public  :: RIVmask
   ! Total amount for each variable
   real(RLEN),allocatable,dimension(:),public  :: D3STATE_tot,D2STATE_tot

#ifdef BFM_NEMO
   !---------------------------------------------
   ! Additional 3D arrays
   !---------------------------------------------
   real(RLEN),allocatable,dimension(:,:,:),public  :: ZEROS
   ! 3D boolean Land-sea mask
   logical,allocatable,dimension(:,:,:),public     :: SEAmask
   ! 3D boolean sea-bottom mask
   logical,allocatable,dimension(:,:,:),public     :: BOTmask
   ! 3D boolean mask of the surface points
   logical,allocatable,dimension(:,:,:),public     :: SRFmask

   !---------------------------------------------
   ! 3D Indices of the wet points
   !---------------------------------------------
   integer,allocatable,dimension(:),public         :: iwet,jwet,kwet

   !---------------------------------------------
   ! Additional integration arrays
   ! for leapfrog scheme
   !---------------------------------------------
   real(RLEN),allocatable,dimension(:,:),public  :: D3STATEB
   real(RLEN),allocatable,dimension(:,:),public  :: D2STATEB

   !---------------------------------------------
   ! Additional allocatable temporary arrays
   !---------------------------------------------
   logical,allocatable,dimension(:),public        :: btmp1D
   logical,allocatable,dimension(:,:),public      :: btmp2D
   logical,allocatable,dimension(:,:,:),public    :: btmp3D
   integer,allocatable,dimension(:),public        :: itmp1D
   integer,allocatable,dimension(:,:),public      :: itmp2D
   integer,allocatable,dimension(:,:,:),public    :: itmp3D
   real(RLEN),allocatable,dimension(:),public     :: rtmp1D
   real(RLEN),allocatable,dimension(:,:),public   :: rtmp2D
   real(RLEN),allocatable,dimension(:,:,:),public :: rtmp3Da
   real(RLEN),allocatable,dimension(:,:,:),public :: rtmp3Db
#endif

#ifdef BFM_POM
   !---------------------------------------------
   ! Additional 3D arrays
   !---------------------------------------------
   real(RLEN),allocatable,dimension(:,:,:),public  :: ZEROS
   ! 3D boolean Land-sea mask
   logical,allocatable,dimension(:,:,:),public     :: SEAmask
   ! 3D boolean sea-bottom mask
   logical,allocatable,dimension(:,:,:),public     :: BOTmask
   ! 3D boolean mask of the surface points
   logical,allocatable,dimension(:,:,:),public     :: SRFmask

   !---------------------------------------------
   ! 3D Indices of the wet points
   !---------------------------------------------
   integer,allocatable,dimension(:),public         :: iwet,jwet,kwet

   !---------------------------------------------
   ! Additional integration arrays
   ! for leapfrog scheme
   !---------------------------------------------
   real(RLEN),allocatable,dimension(:,:),public  :: D3STATEB
   real(RLEN),allocatable,dimension(:,:),public  :: D2STATEB

   !---------------------------------------------
   ! Additional allocatable temporary arrays
   !---------------------------------------------
   logical,allocatable,dimension(:),public        :: btmp1D
   logical,allocatable,dimension(:,:),public      :: btmp2D
   logical,allocatable,dimension(:,:,:),public    :: btmp3D
   integer,allocatable,dimension(:),public        :: itmp1D
   integer,allocatable,dimension(:,:),public      :: itmp2D
   integer,allocatable,dimension(:,:,:),public    :: itmp3D
   real(RLEN),allocatable,dimension(:),public     :: rtmp1D
   real(RLEN),allocatable,dimension(:,:),public   :: rtmp2D
   real(RLEN),allocatable,dimension(:,:,:),public :: rtmp3Da
   real(RLEN),allocatable,dimension(:,:,:),public :: rtmp3Db
#endif

!
! !REVISION HISTORY:
!  Author(s): Marcello Vichi and Piet Ruardij
!  Uses functions and portions of code from GOTM
!  by Hans Burchard and Karsten Bolding
!
! !LOCAL VARIABLES:
!
! !BUGS
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bfm module
!
! !INTERFACE:
   subroutine init_bfm(namlst)
!
! !DESCRIPTION:
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,            &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,    &
                  NO_D2_BOX_STATES, NO_BOXES_XY,         &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_STATES, Depth, NO_D3_BOX_FLUX,      &
                  NO_D2_BOX_FLUX
   use global_mem, only: LOGUNIT
   use time, only: bfmtime
#if defined key_obcbfm
   use global_mem, only: LOGUNITOBC
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  Adapted from GOTM code by Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: rc,n
   character(len=PATH_MAX)   :: logfname
#if defined key_obcbfm
   character(len=PATH_MAX)   :: logfnameobc
#endif
   namelist /bfm_nml/ bio_calc,bio_setup,bfm_init,bfm_rstctl,   &
                      out_fname,out_dir,out_units,out_title,    &
                      out_delta,out_secs,bioshade_feedback,     &
                      parallel_log
!EOP
!-----------------------------------------------------------------------
!BOC


   !---------------------------------------------
   ! Provide sensible values for namelist parameters
   !---------------------------------------------
   bio_calc   = .TRUE.
   bio_setup  = 1
   bfm_init   = 0
   bfm_rstctl = .FALSE.
   out_fname  = 'bfm'
   rst_fname  = 'bfm_restart'
   out_dir    = '.'
   out_title  = 'Another great BFM simulation!'
   out_units  = 0
   out_delta  = 100
   out_secs   = 100
   bioshade_feedback=.FALSE.

   !---------------------------------------------
   !  Open and read the namelist
   !---------------------------------------------
   open(namlst,file='BFM_General.nml',action='read',status='old',err=99)
   read(namlst,nml=bfm_nml,err=98)
   close(namlst)

#ifdef BFM_PARALLEL
   LEVEL2 "BFM is running in Parallel"
   parallel = .TRUE.
   ! variable parallel_rank must have been assigned previously
   ! in the coupling with the ocean model 
   ! check if logs have to be produced for each process
   ! and provide a different log file name 
   LOGUNIT = 1069 + parallel_rank
   write(str,'(I4.4)') parallel_rank
    if (parallel_log) then
       if (parallel_rank == 0) then
          bfm_lwp = .TRUE.
          logfname = 'bfm.log'
       else 
          bfm_lwp = .FALSE.
          logfname = '/dev/null'
       end if
    else 
       ! logs are produced for every process
       bfm_lwp = .TRUE.
       logfname = 'bfm_'//str//'.log'
    end if
    open(LOGUNIT,file=logfname,action='write',  &
        form='formatted',err=100)

   ! provide different file names for each process domain
   out_fname = trim(out_fname)//'_'//str
   rst_fname = trim(rst_fname)//'_'//str

   LEVEL1 'init_bfm'
   LEVEL3 "Producing log for process rank:",parallel_rank

#else
   LEVEL1 'init_bfm'
   LOGUNIT = 1069
   logfname = 'bfm.log'
   bfm_lwp = .TRUE.
   open(LOGUNIT,file=logfname,action='write',  &
        form='formatted',err=100)

#endif
   !-------------------------------------------------------
   ! Write to log bfmtime setting
   !-------------------------------------------------------
   LEVEL2 'BFM time informations:'
   WRITE(LOGUNIT,*) 'Start Date, Julianday0, JuliandayEnd, step0, timestep, stepnow, stepEnd'
   WRITE(LOGUNIT,*) bfmtime   
   !
   LEVEL2 "Writing NetCDF output to file: ",trim(out_fname)
   LEVEL3 "Output frequency every ",out_delta,"time-steps"

#ifndef INCLUDE_BEN
   ! force bio_setup = 1 when benthic memory is disabled with macro
   bio_setup=1
#endif

   select case (bio_setup)
      case (0)
      case (1) ! Pelagic only
        LEVEL2 "Using a Pelagic setup (bio_setup=1)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
      case (2) ! Benthic only
        LEVEL2 "Using a Benthic-only setup (bio_setup=2)"
        LEVEL3 'benthic variables =',NO_D2_BOX_STATES
        LEVEL3 'benthic diagnostic variables=', NO_D2_BOX_DIAGNOSS
      case (3) ! Pelagic-Benthic coupling
        LEVEL2 "Using a Pelagic-Benthic coupled setup (bio_setup=3)"
        LEVEL3 'pelagic variables =',NO_D3_BOX_STATES
        LEVEL3 'pelagic transported variables ='
        LEVEL3 'pelagic diagnostic variables =', NO_D3_BOX_DIAGNOSS
        LEVEL3 'benthic variables =',NO_D2_BOX_STATES
        LEVEL3 'benthic diagnostic variables=', NO_D2_BOX_DIAGNOSS
   end select

   LEVEL2 'Dimensional informations:'
   LEVEL3 'NO_BOXES_X=',NO_BOXES_X
   LEVEL3 'NO_BOXES_Y=',NO_BOXES_Y
   LEVEL3 'NO_BOXES_Z=',NO_BOXES_Z
   LEVEL3 'NO_BOXES=',NO_BOXES
   LEVEL3 'NO_BOXES_XY=',NO_BOXES_XY
   LEVEL3 'NO_STATES=',NO_STATES
   LEVEL3 'Step 1 of BFM initialisation done ...'
   ! dimension lengths used in the netcdf output
   lon_len = NO_BOXES_X
   lat_len = NO_BOXES_Y
   depth_len = NO_BOXES_Z
   ocepoint_len = NO_BOXES
   surfpoint_len = NO_BOXES_XY
   botpoint_len = NO_BOXES_XY

   !---------------------------------------------
   ! Allocate arrays with attributes of state
   ! variables
   !---------------------------------------------
   ! total number of output states
   n = NO_D3_BOX_STATES+NO_D3_BOX_FLUX+NO_D3_BOX_DIAGNOSS+ &
      NO_D2_BOX_STATES+NO_D2_BOX_FLUX+NO_D2_BOX_DIAGNOSS
   allocate(var_ids(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_ids'
   var_ids=0;
   allocate(var_ave(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_ave'
   var_ave=.false.;
   allocate(var_names(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_names'
   allocate(var_units(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_units'
   allocate(var_long(1:n),stat=rc)
   if (rc /= 0) stop 'init_bfm(): Error allocating var_long'
   ! temporary diagnostic variable
   allocate(c1dim(1:NO_BOXES),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating c1dim'

   return
98 FATAL 'I could not read BFM_General.nml'
   stop 'init_bfm'
99 LEVEL2 'I could not open BFM_General.nml'
   LEVEL2 'Simulation starting without the BFM'
   bio_calc = .false.
100 FATAL 'Cannot create log file: ',trim(logfname)
   stop 'init_bfm'

  end subroutine init_bfm
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
   function find(vector,nt)
!
! !DESCRIPTION:
!  Finds the location of true elements in logical arrays
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   logical,intent(IN) :: vector(:)
   integer,intent(IN) :: nt   ! number of true elements in vector
                              ! nt = count(vector)
                              ! enter as an argument for check
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer            :: find(nt)
!
! !REVISION HISTORY:
!  Author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer            :: l,m,n
!
!EOP
!-----------------------------------------------------------------------
!BOC

    if (nt /= count(vector)) stop '#### Error in find: check the input array ####'
    m = size(vector,1)
    n = 1
    do l = 1,m
      if ( vector(l) ) then
        find(n) = l
        n = n + 1
      end if
    end do

   return
   end function find

!EOC


!-----------------------------------------------------------------------

   end module api_bfm

!-----------------------------------------------------------------------
!Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!Copyright (C) 2006 - Marcello Vichi

