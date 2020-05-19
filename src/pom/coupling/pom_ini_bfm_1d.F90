#include "INCLUDE.h"
#include "cppdefs.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: pom_ini_bfm
!
!
!DESCRIPTION:
!  Initialise the BFM in POM
!  Main communication of array dimensions between
!  BFM and POM
!  Initialisation of variables and netcdf output
!  THIS IS THE 1D VERSION
!
! !INTERFACE:
!
  subroutine pom_ini_bfm_1d
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   use mem
!
   use POM, ONLY: DTI, IEND, KB, H, DZ, ihotst, ALAT, ALON, ZZ, DZZ
!
   use global_mem, only:ZERO
!
   use api_bfm
!
   use Service, ONLY: savef, deltat, nitend
!
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   use netcdf_bfm, only: init_netcdf_rst_bfm,read_rst_bfm

!------------------------------------------------------------------------!
!
!BOC
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Implicit typing is never allowed
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   IMPLICIT NONE
!
!  ----COUNTER-----
!
   integer    :: k,i
!
!  ----READING UNITS-----
!
   integer,parameter    :: namlst=10,unit=11
!
!  ----OCEAN WET POINTS (IN 1D = 1)-----
!
   integer,allocatable  :: ocepoint(:), &
!
!  -----OCEAN SURFACE POINTS (IN 1D=1)-----
!
                          surfpoint(:), &
!
!  ----OCEAN BOTTOM POINTS (IN 1 D=1)-----
!
                          botpoint(:)
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Model timestep
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
     deltat=dti
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Iterations needed for an idays simulations
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
     nitend=iend
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !out_delta = saving frequency
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
     out_delta = savef
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Set the BFM dimensions
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!    *************************************************************************
!    *************************************************************************
!    **                                                                     **
!    ** SINCE THIS IS A 1D IMPLEMENTATION THE SIZE OF THE ARRAYS IS SHRUNK  **
!    ** TO A 1D VECTOR  ALONG the "VERTICAL" DIMENSION (1, 1, KB-1)         **
!    **                                                                     **
!    *************************************************************************
!    *************************************************************************
!
     NO_BOXES=KB-1
     NO_BOXES_X=1
     NO_BOXES_Y=1
     NO_BOXES_Z=NO_BOXES
     NO_BOXES_XY=1
!
!   -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!    Allocate masks for array packing
!   -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   allocate(SEAmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
   allocate(BOTmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
   allocate(SRFmask(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
!
   SEAmask = .TRUE.
   BOTmask = .TRUE.
   SRFmask = .TRUE.
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !allocate ancillary pack mask
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   allocate(ZEROS(NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z))
   ZEROS = ZERO
!
#ifdef INCLUDE_BEN
!
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES_BEN
   NO_BOXES_Z_BEN  = 0
   NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
   NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN
!
#else
!
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES
!
#endif
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Compressed coordinates for netcdf output
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   lon_len = NO_BOXES_X
   lat_len = NO_BOXES_Y
!
   allocate(ocepoint(NO_BOXES))
   ocepoint = 1
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Prepares the array containing the indices of the
   ! elements in pelagic BFM 1D arrays that have a benthic layer
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   allocate(BOTindices(NO_BOXES_XY)); BOTindices = KB-1
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Prepares the array containing the indices of the
   ! elements in pelagic BFM 1D arrays that have a surface layer
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   allocate(SRFindices(NO_BOXES_XY)); SRFindices = 1
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Initialise ancillary arrays for output
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   call init_bfm(namlst)
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Initialise state variable names and diagnostics
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   call set_var_info_bfm
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Allocate memory and give initial values
   ! the argument list is mandatory with BFM
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   call init_var_bfm(namlst,'bfm.nml',unit,bio_setup)
!
!       -----THE LEADING RESTART FLAG IS IHOTST, SO THE VALUE IN bfm_init-------
!                          -----IS OVERWRITTEN-----
!
        bfm_init=ihotst
!
!       -----SET THE THICKNESS OF THE KB-1 LAYERS-----
!
        do k = 1 , NO_BOXES_Z
               Depth(k) = dz(k)*h
        end do
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! initialise netcdf output
   ! MAV: CHECK THE MAPPING OF OCEPOINT
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   call calcmean_bfm(INIT)

   call init_netcdf_bfm(title=out_title,start_time='01-01-0000', &
             time_unit=0,lat=alat,lon=alon,z=zz,dz=dzz,      &
             oceanpoint=ocepoint,                  &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/))
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Allocate and initialise additional
   ! integration arrays
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
   allocate(D3STATEB(NO_D3_BOX_STATES,NO_BOXES))
   D3STATEB = ZERO
!
#ifdef INCLUDE_BEN
!
   allocate(D2STATEB_BEN(NO_D2_BOX_STATES_BEN,NO_BOXES_XY))
   D2STATEB_BEN = ZERO
!
#endif

!  -----DEFINE INITIAL CONDITIONS-----

#ifdef BFM17

   call set_initial_conditions_bfm17

#else

   call set_initial_conditions_bfm56

#endif
!
   call init_save_bfm
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Initialise prior time step for leap-frog
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
     D3STATEB = D3STATE
!
#ifdef INCLUDE_BEN
     D2STATEB_BEN = D2STATE_BEN
#endif
!
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Read restart (Bfm_init = 1 in bfm.nml)
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!======================GIULIA==============================================
   call init_netcdf_rst_bfm(rst_fname,start_time='01-01-0000',  &
             time_unit=0,lat=alat,lon=alon,z=zz,dz=dzz,   &
             oceanpoint=ocepoint,      &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/))
!==========================================================================
!
   !--------------------------------------------------
   ! Read restart file (if flag)
   ! Overwrite previous initialization
   !--------------------------------------------------
   write(6,*) 'before read rst'
   if (bfm_init == 1) call read_rst_bfm(rst_fname)
!
   return
!
 end subroutine pom_ini_bfm_1d
!EOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
