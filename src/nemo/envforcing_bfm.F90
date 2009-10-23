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
   use global_mem, only:RLEN,ZERO,LOGUNIT
   use mem_param,  only: p_PAR, p_small
   use mem,        only: xEPS, ESS, ETW, ESW, EWIND,    &
                         Depth, EIR, ERHO, EICE,        &
                         NO_BOXES, NO_BOXES_XY
   use api_bfm
#ifdef INCLUDE_PELCO2
   use mem,        only: EPCO2air
   use mem_CO2,    only: pco2air
#endif
! OPA modules
   use oce_trc
   use trc_oce, only: etot3
IMPLICIT NONE
! OPA domain substitutions
#include "domzgr_substitute.h90"
!
! !INPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:

! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer             :: i,j,k,n

!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Assign temperature, salinity and density
   !---------------------------------------------
#ifdef USEPACK
      ETW = pack(tn_io,SEAmask)
      ESW = pack(sn_io,SEAmask)
      ERHO = pack(rhop_io,SEAmask)
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily || defined key_flx_ecmwf_mfs
   !---------------------------------------------
   ! Assign wind speed
   !---------------------------------------------
      EWIND = pack(vatm_io,SRFmask(:,:,1) )
#else
      !MAV: this must be temporary! 
      EWIND = 5.0_RLEN
#endif
   !---------------------------------------------
   ! Assign Sea-ice cover
   !---------------------------------------------
      EICE = pack(freeze_io,SRFmask(:,:,1) )
#else
      DO n = 1,NO_BOXES
         ETW(n) = tn_io(iwet(n),jwet(n),kwet(n))
         ESW(n) = sn_io(iwet(n),jwet(n),kwet(n))
         ERHO(n) = rhop_io(iwet(n),jwet(n),kwet(n))
      END DO

      DO n = 1,NO_BOXES_XY
         !---------------------------------------------
         ! Assign wind speed
         !---------------------------------------------
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily || defined key_flx_ecmwf_mfs
         EWIND(n) = vatm_io(iwet(n),jwet(n))

#else
         !MAV: this must be temporary! 
         EWIND(n) = 5.0_RLEN
#endif
         !---------------------------------------------
         ! Assign Sea-ice cover
         !---------------------------------------------
         EICE(n) = freeze_io(iwet(n),jwet(n))
      END DO

#endif



#ifdef INCLUDE_PELCO2
   !---------------------------------------------
   ! Assign atmospheric pCO2
   !---------------------------------------------
      EPCO2air = pco2air !costant
      !EPCO2air = pack(,SRFmask(:,:,1) )
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
#ifdef USEPACK
      rtmp3Db = unpack(xEPS,SEAmask,ZEROS)
#else
      DO n = 1,NO_BOXES
         rtmp3Db(iwet(n),jwet(n),kwet(n)) = xEPS(n)
      END DO
#endif
      do k = 1,jpk-1
         do j = 1,jpj
            do i = 1,jpi
               rtmp3Da(i,j,k+1) = rtmp3Da(i,j,k)* &
                       exp(-rtmp3Db(i,j,k)*fse3t(i,j,k))
            end do 
         end do 
      end do 
#ifdef USEPACK
      EIR = pack(rtmp3Da,SEAmask)
#else
      DO n = 1,NO_BOXES
         EIR(n) = rtmp3Da(iwet(n),jwet(n),kwet(n))
      END DO
#endif

   !---------------------------------------------
   ! bioshading is stored to be passed to OPA
   ! (converted back to W m-2)
   ! The dynamics of active tracers is indeed
   ! computed after the BFM.
   ! It already includes the abiotic part, so that
   ! the BFM extinction coefficients are used
   ! and not the OPA ones
   !---------------------------------------------
      etot3(:,:,:) = rtmp3Da(:,:,:)*E2W/p_PAR

   !---------------------------------------------
   ! Deallocate temporary arrays
   !---------------------------------------------
      deallocate(rtmp3Da)
      deallocate(rtmp3Db)

   end subroutine envforcing_bfm
!EOC
!-----------------------------------------------------------------------
