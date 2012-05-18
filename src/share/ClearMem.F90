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
    ! from AllocateMem
    deallocate(flx_calc_nr)
    deallocate(flx_CalcIn)
    deallocate(flx_option)
    deallocate(flx_t)
    deallocate(flx_SS)
    deallocate(flx_states)
    deallocate(flx_ostates)
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
