#include "cppdefs.h"
#include "INCLUDE.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ResetFluxes.F90
!
! !INTERFACE
subroutine ResetFluxes
!
! !DESCRIPTION
!  Reset the arrays for the next integration
!
! !USES
   use global_mem, only:ZERO
   use mem, ONLY: NO_D3_BOX_STATES,D3SOURCE, D3SINK, &
         PELBOTTOM, PELSURFACE
#ifdef INCLUDE_BEN
   use mem, ONLY: NO_D2_BOX_STATES, D2SINK, D2SOURCE
#endif
   implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
!
   ! Reset source term arrays 
   ! only the diagonal
   do i=1,NO_D3_BOX_STATES
      D3SOURCE(i,i,:) = ZERO
#ifndef ONESOURCE
      D3SINK(i,i,:) = ZERO
#endif
   end do
#ifdef INCLUDE_BEN
   do i=1,NO_D2_BOX_STATES
      D2SOURCE(i,i,:) = ZERO
#ifndef ONESOURCE
      D2SINK(i,i,:) = ZERO
#endif
   end do
#endif

   PELBOTTOM(:,:) = ZERO
   PELSURFACE(:,:) = ZERO

end subroutine ResetFluxes
!EOC
!-----------------------------------------------------------------------
