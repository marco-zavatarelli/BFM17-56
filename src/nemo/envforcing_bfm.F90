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
                         Depth, EIR, ERHO, EICE, EPR,   &
                         NO_BOXES, NO_BOXES_XY
   use api_bfm
   use SystemForcing, only : FieldRead
#ifdef INCLUDE_PELCO2
   use mem_CO2,    only: AtmCO20, AtmCO2, AtmSLP, AtmTDP
#endif
! OPA modules
   use oce_trc
   use trc_oce, only: etot3
! MFS bulk 
#if defined  INCLUDE_PELCO2 && defined USE_MFSBULK
   use sbcblk_mfs, ONLY: sf 
#endif

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
      ETW = pack(tsn(:,:,:,jp_tem),SEAmask)
      ESW = pack(tsn(:,:,:,jp_sal),SEAmask)
      ERHO = pack(rhop(:,:,:),SEAmask)
   !---------------------------------------------
   ! Assign wind speed
   !---------------------------------------------
      EWIND = pack(wndm(:,:),SRFmask(:,:,1) )
   !---------------------------------------------
   ! Assign Sea-ice cover
   !---------------------------------------------
      EICE = pack(fr_i(:,:),SRFmask(:,:,1) )
#else
      DO n = 1,NO_BOXES
         ETW(n) = tsn(iwet(n),jwet(n),kwet(n),jp_tem)
         ESW(n) = tsn(iwet(n),jwet(n),kwet(n),jp_sal)
         ERHO(n) = rhop(iwet(n),jwet(n),kwet(n))
      END DO

      DO n = 1,NO_BOXES_XY
         !---------------------------------------------
         ! Assign wind speed
         !---------------------------------------------
         EWIND(n) = wndm(iwet(n),jwet(n))
         !---------------------------------------------
         ! Assign Sea-ice cover
         !---------------------------------------------
         EICE(n) = fr_i(iwet(n),jwet(n))
      END DO

#endif

#ifdef INCLUDE_PELCO2
   !---------------------------------------------
   ! Assign atmospheric pCO2
   !---------------------------------------------
   if (AtmCO2%init .ne. 0) then
      call FieldRead(AtmCO2)
   endif
   ! Water column pressure 
   ! (need better approximation to convert from m to dbar)
   EPR = pack(fsdept(:,:,:),SEAmask)
#ifdef USE_MFSBULK
   !
   ! Atmospheric sea level pressure (MFS index jp_msl  = 4)
   if ( allocated(AtmSLP%fnow))  then 
      if (AtmSLP%init .eq.4) AtmSLP%fnow = pack(sf(4)%fnow(:,:,1),SRFmask(:,:,1) )
   endif
   !
   ! Atmospheric Dew Point Temperature (MFS index jp_rhm  = 6)
   if ( allocated(AtmTDP%fnow)) then 
      if (AtmTDP%init .eq.4) AtmTDP%fnow = pack(sf(6)%fnow(:,:,1),SRFmask(:,:,1) )
   endif
#endif
#endif

   !---------------------------------------------
   ! Temporary 3D array for the storage of the 
   ! light environment.
   ! Assign surface irradiance to the first layer.
   ! (converted to PAR and uE,
   ! add parametric zero for nighttime)
   ! Initialise the bioshading array if ln_qsr_bio
   !---------------------------------------------
      allocate(rtmp3Da(jpi,jpj,jpk)); rtmp3Da = ZERO
      rtmp3Da(:,:,1) = p_PAR*(qsr(:,:)+p_small)/E2W 
      if (ln_qsr_bio) etot3(:,:,1) = qsr(:,:)

   !---------------------------------------------
   ! Compute extinction coefficient
   !---------------------------------------------
      call CalcVerticalExtinction( )

   !---------------------------------------------
   ! temporarely unpack it from the 1D array 
   ! to 3D grid (apply land-sea mask)
   !---------------------------------------------
      allocate(rtmp3Db(jpi,jpj,jpk))
#ifdef USEPACK
      rtmp3Db = unpack(xEPS,SEAmask,ZEROS)
#else
      DO n = 1,NO_BOXES
         rtmp3Db(iwet(n),jwet(n),kwet(n)) = xEPS(n)
      END DO
#endif

   !---------------------------------------------
   ! Compute the light climate and repack
   ! Note that in BFM light is defined at the
   ! top of each level (W grid in OPA)
   !---------------------------------------------
   ! Bioshading is also stored if ln_qsr_bio 
   ! is true in namelist and passed to OPA
   ! (converted back to W m-2 and on the T grid)
   ! The dynamics of active tracers is 
   ! computed after the BFM call.
   ! It already includes the abiotic part, so that
   ! the BFM extinction coefficients are used
   ! and not the OPA ones
   !---------------------------------------------
      do k = 1,jpkm1
         do j = 1,jpj
            do i = 1,jpi
               rtmp3Da(i,j,k+1) = rtmp3Da(i,j,k)*               &
                                  exp(-rtmp3Db(i,j,k)*fse3w(i,j,k))
               if (ln_qsr_bio) & 
			       etot3(i,j,k+1) = etot3(i,j,k)*    &
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
   ! Deallocate temporary arrays
   !---------------------------------------------
      deallocate(rtmp3Da)
      deallocate(rtmp3Db)

   end subroutine envforcing_bfm
!EOC
!-----------------------------------------------------------------------
