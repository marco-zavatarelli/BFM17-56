#include "cppdefs.h"
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
   use mem, ONLY: NO_D2_BOX_STATES,NO_BOXES_XY,D2SOURCE, &
         NO_D3_BOX_STATES,NO_BOXES,D3SOURCE, D3SINK, D2SINK, &
         PELBOTTOM, PELSURFACE
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
#ifdef DEBUG
   D3SOURCE = ZERO
   D3SINK = ZERO
   D2SOURCE = ZERO
   D2SINK = ZERO
#else
   ! only the diagonal
   do i=1,NO_D3_BOX_STATES
      D3SOURCE(i,i,:) = 0.0
      D3SINK(i,i,:) = 0.0
   end do
   do i=1,NO_D2_BOX_STATES
      D2SOURCE(i,i,:) = 0.0
      D2SINK(i,i,:) = 0.0
   end do
   PELBOTTOM = ZERO
   PELSURFACE = ZERO
#endif

end subroutine ResetFluxes
!EOC
!-----------------------------------------------------------------------
