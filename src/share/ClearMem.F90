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
#ifndef NOT_STANDALONE
     deallocate(D3STATE)
     deallocate(D3SOURCE)
     deallocate(D3SINK)
     deallocate(D3STATETYPE)
     deallocate(D2STATE)
     deallocate(D2SOURCE)
     deallocate(D2SINK)
     deallocate(D2STATETYPE)
     deallocate(D3DIAGNOS)
     deallocate(D2DIAGNOS)

#endif

end subroutine ClearMem
!EOC
!-----------------------------------------------------------------------
