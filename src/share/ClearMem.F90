#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ClearMem.F90
!
! !INTERFACE
subroutine ClearMem
!
! !DESCRIPTION
!  Deallocate all the arrays
!
! !USES
   use mem
   use api_bfm
   implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV)
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
!

   integer :: i,j,origin,destination

#ifndef D1SOURCE
    ! from AllocateMem
    deallocate(flx_calc_nr)
    deallocate(flx_CalcIn)
    deallocate(flx_option)
    deallocate(flx_t)
    deallocate(flx_SS)
    deallocate(flx_states)
    deallocate(flx_ostates)
#else
    ! free willy, free the fluxes
    origin=0
    do i=stPelStateS,stPelStateE
       origin=origin+1
       destination=0
       do j=stPelStateS,stPelStateE
          destination=destination+1
          if( allocated(D3FLUX_MATRIX(origin,destination)%p) ) deallocate(D3FLUX_MATRIX(origin,destination)%p)
       end do
    end do
    deallocate(D3FLUX_MATRIX)
    deallocate(D3FLUX_FUNC)
#endif

    ! from api_bfm 
    deallocate(var_ids)
    deallocate(var_ave)
    deallocate(var_names)
    deallocate(var_units)
    deallocate(var_long)
    deallocate(c1dim)
    if (allocated(D3ave)) deallocate(D3ave)
    if (allocated(D2ave)) deallocate(D2ave)

#ifndef NOT_STANDALONE
     deallocate(D3STATE)
     deallocate(D3SOURCE)
#ifndef ONESOURCE
     deallocate(D3SINK)
#endif
     deallocate(D3STATETYPE)
     deallocate(D3DIAGNOS)
     deallocate(D2DIAGNOS)
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
     deallocate(D2STATE)
     deallocate(D2SOURCE)
     deallocate(D2SINK)
     deallocate(D2STATETYPE)
#endif

#endif

end subroutine ClearMem
!EOC
!-----------------------------------------------------------------------
