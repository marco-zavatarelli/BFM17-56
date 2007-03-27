#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Light and other environmental forcing used in the BFM
!
! !INTERFACE
   subroutine envforcing_bfm()
!
! !DESCRIPTION
!
! !USES
! BFM modules
   use constants,  only: E2W
   use global_mem, only:RLEN,ZERO
   use mem_param,  only: p_PAR, p_small
   use mem,        only: xEPS, ESS, ETW, ESW, EWIND,    &
                        Depth, EIR
   use api_bfm
   use oce_trc

IMPLICIT NONE
!
! !INPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:

! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer             :: k

!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Assign temperature and salinity
   !---------------------------------------------
      ETW = pack(tn,SEAmask)
      ESW = pack(sn,SEAmask)

#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
   !---------------------------------------------
   ! Assign wind speed
   !---------------------------------------------
      EWIND = pack(vatm,SRFmask(:,:,1) )
#else
      !write(numout,*) 'BFM WARNING: wind speed is not defined in NEMO!'
#endif

   !---------------------------------------------
   ! Temporary 3D array for the storage of the 
   ! light environment.
   ! Assign surface irradiance to the first layer.
   ! (converted to PAR and uE,
   ! add parametric zero for nighttime)
   !---------------------------------------------
      allocate(rtmp3Da(jpi,jpj,jpk)); rtmp3Da = ZERO
      rtmp3Da(:,:,1) = p_PAR*(qsr(:,:)+p_small)/E2W 

   !---------------------------------------------
   ! Compute extinction coefficient
   !---------------------------------------------
      call CalcVerticalExtinction( )

   !---------------------------------------------
   ! temporarely unpack it from the 1D array 
   ! to 3D grid (apply land-sea mask)
   ! Then compute the light climate and repack
   !---------------------------------------------
      allocate(rtmp3Db(jpi,jpj,jpk))
      rtmp3Db = unpack(xEPS,SEAmask,ZEROS)
      do k = 1,jpk-1
         rtmp3Da(:,:,k+1) = rtmp3Da(:,:,k)*exp(-rtmp3Db(:,:,k)*e3t_0(k))
      end do 
      EIR = pack(rtmp3Da,SEAmask)
      deallocate(rtmp3Da)
      deallocate(rtmp3Db)

   !---------------------------------------------
   ! bioshade is instead derived in the
   ! middle of the layer and it's non-dimensional
   !---------------------------------------------
   !if (bioshade_feedback) &
   !  bioshade(1:nlev) =  EIR(:)*exp(-xEPS(:)*Depth(:)*0.5)/ EIR(nlev)

   end subroutine envforcing_bfm
!EOC
!-----------------------------------------------------------------------
