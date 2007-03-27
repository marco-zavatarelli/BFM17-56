#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_ini_bfm.F90
!
! !INTERFACE:
   subroutine trc_ini_bfm()
!
! !DESCRIPTION:
!  Initialise the BFM in NEMO
!  Main communication of array dimensions between
!  BFM and NEMO
!  Initialization of variables and netcdf output
!
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_STATES,Depth,Depth_ben, D3STATE, D2STATE
   use global_mem, only:RLEN,ZERO,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
   use api_bfm
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   ! NEMO modules
   USE trctrp_lec, only: l_trczdf_exp,ln_trcadv_cen2,ln_trcadv_tvd    
   use oce_trc

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   !! * Substitutions
#include "domzgr_substitute.h90"

   integer    :: i,j,k
   integer    :: status
   integer,parameter    :: namlst=10,unit=11
   integer,allocatable  :: ocepoint(:),surfpoint(:),botpoint(:)
   logical,allocatable  :: mask1d(:)
!EOP
!-----------------------------------------------------------------------
!BOC
   allocate(SEAmask(jpi,jpj,jpk))  ! allocate  masks for
   allocate(BOTmask(jpi,jpj,jpk))  ! array packing
   allocate(SRFmask(jpi,jpj,jpk))
   SEAmask = .FALSE.
   BOTmask = .FALSE.
   SRFmask = .FALSE.

   allocate(ZEROS(jpi,jpj,jpk))   ! allocate ancillary pack mask
   ZEROS = ZERO
   allocate(rtmp3Da(jpi,jpj,jpk)) ! allocate temporary 3D array

   where (tmask > ZERO)         ! assign logical Land-Sea mask
     SEAmask = .TRUE.
   elsewhere
     SEAmask = .FALSE.
   end where


   !-------------------------------------------------------
   ! Prepares the masks for the bottom and surface
   ! grid points.
   ! 3D boolean arrays with .T. at the first
   ! and last ocean layers only
   !-------------------------------------------------------

   do j = 1,jpj
      do i = 1,jpi
!if mbathy is defined uncomment the next lines
! minimum of mbathy is now 2. Check if still works
!          if (mbathy(i,j) > 0) then
!            SRFmask(i,j,1) =  .TRUE.
!            BOTmask(i,j,mbathy(i,j)) =  .TRUE.
!          end if
         rtmp3Da(i,j,:) = fse3t(i,j,:)
         if (SEAmask(i,j,1)) then
            SRFmask(i,j,1) = .TRUE.
            do k = 1,jpk
               if (.not.SEAmask(i,j,k)) then
                  BOTmask(i,j,k-1) =  .TRUE.
                  exit
               end if
            end do
         end if
      end do
   end do

   !---------------------------------------------
   ! Set the dimensions
   !---------------------------------------------
   NO_BOXES_X  = jpi
   NO_BOXES_Y  = jpj
   NO_BOXES_Z  = jpk
   NO_BOXES    = count(SEAmask)
   NO_BOXES_XY = count(SRFmask)
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES * NO_BOXES_XY

   !-------------------------------------------------------
   ! Compressed coordinates for netcdf output
   !-------------------------------------------------------
   lon_len = NO_BOXES_X
   lat_len = NO_BOXES_Y
   allocate(ocepoint(NO_BOXES))
   allocate(mask1d(1:NO_BOXES_X*NO_BOXES_Y*NO_BOXES_Z))
   mask1d = reshape(SEAmask,(/NO_BOXES_X*NO_BOXES_Y*NO_BOXES_Z/))
   ocepoint = find(mask1d,NO_BOXES)
   deallocate(mask1d)

   allocate(surfpoint(NO_BOXES_XY))
   allocate(botpoint(NO_BOXES_XY))


!MAV:WIND?
   !   allocate(EWind(NO_BOXES)) ! allocate wind as 1-D variable

   !-------------------------------------------------------
   ! allocate the arrays for surface and bottom boundary conditions
   ! At the moment only for inorganic components (nutrients and gases)
   !-------------------------------------------------------
!      allocate(bc2D_bio_surf(NO_D3_BOX_STATES,NO_BOXES),stat=status)
!      if (status /= 0)
!     &  stop "# FATAL ERROR in ini_trc_bfm: Error allocating bc2D_bio_surf"
!      bc2D_bio_surf = ZERO
!      allocate(bc2D_bio_bot(NO_D3_BOX_STATES,NO_BOXES),stat=status)
!      if (status /= 0)
!     &  stop "# FATAL ERROR in ini_trc_bfm: Error allocating bc2D_bio_bot"
!      bc2D_bio_surf = ZERO

   !-------------------------------------------------------
   ! Prepares the array containing the indices of the
   ! elements in pelagic BFM 1D arrays that have a benthic layer
   !-------------------------------------------------------
   allocate(BOTindices(NO_BOXES)); BOTindices = 0
   allocate(btmp1D(NO_BOXES))
   btmp1D = pack(BOTmask,SEAmask)
   BOTindices = find(btmp1D,NO_BOXES_XY)

   !-------------------------------------------------------
   ! Prepares the array containing the indices of the
   ! elements in pelagic BFM 1D arrays that have a surface layer
   !-------------------------------------------------------
   allocate(SRFindices(NO_BOXES)); SRFindices = 0
   btmp1D = pack(SRFmask,SEAmask)
   SRFindices = find(btmp1D,NO_BOXES_XY)
   deallocate(btmp1D)

   !---------------------------------------------
   ! Initialise ancillary arrays for output
   !---------------------------------------------
   call init_bfm(namlst)

   !---------------------------------------------
   ! Initialise state variable names and diagnostics
   !---------------------------------------------
   call set_var_info_bfm

   !---------------------------------------------
   ! Allocate memory and give initial values
   !---------------------------------------------
   ! the argument list is mandatory with BFM
   call init_var_bfm(namlst,'bfm.nml',unit,bio_setup)

   !-------------------------------------------------------
   ! Prepares the BFM 1D array containing the thickness
   ! at each sea gridpoint (Depth(NO_BOXES))
   !---------------------------------------------
   Depth = pack(rtmp3Da,SEAmask)
   deallocate(rtmp3Da)

   !---------------------------------------------
   ! initialise netcdf output
   ! MAV: MISSING THE MAPPING OF OCEPOINT
   !---------------------------------------------
   call calcmean_bfm(INIT)
   call init_netcdf_bfm(out_title,'01-01-0000',0,  &
             lat2d=gphit,lon2d=glamt,z=gdept_0,      &
             oceanpoint=ocepoint,                  &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/),  &
             mask3d=tmask)
   call init_save_bfm

   if ( l_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) then
      !---------------------------------------------
      ! Allocate and initialise additional 
      ! integration arrays
      !---------------------------------------------
      allocate(D3STATEB(NO_D3_BOX_STATES,NO_BOXES))
      allocate(D2STATEB(NO_D2_BOX_STATES,NO_BOXES))

      !---------------------------------------------
      ! Initialise prior time step for leap-frog 
      !---------------------------------------------
      D3STATEB = D3STATE
      D2STATEB = D2STATE
   end if

   return

   end subroutine trc_ini_bfm
!EOC

