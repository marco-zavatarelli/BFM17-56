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
#ifdef NOPOINTERS
  use mem
#else
   use mem, ONLY: NO_D3_BOX_STATES,D3SOURCE, D3SINK, &
         PELBOTTOM, PELSURFACE
#ifdef INCLUDE_BEN
   use mem, ONLY: NO_D2_BOX_STATES, D2SINK, D2SOURCE
#endif
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
#ifdef D1SOURCE
   ! Reset the 1-dimensional source term arrays 
   D3SOURCE(:,:) = ZERO
#  ifdef INCLUDE_BEN
   D2SOURCE(:,:) = ZERO
#  endif
#else
#  ifdef  ONESOURCE
   ! Reset the whole source term array 
   D3SOURCE(:,:,:) = ZERO
#     ifdef INCLUDE_BEN
   D2SOURCE(:,:,:) = ZERO
#     endif
#  else
   ! Reset source and sink term arrays 
   ! only the diagonal
   do i=1,NO_D3_BOX_STATES
      D3SOURCE(i,i,:) = ZERO
      D3SINK(i,i,:) = ZERO
   end do
#    ifdef INCLUDE_BEN
   do i=1,NO_D2_BOX_STATES
      D2SOURCE(i,i,:) = ZERO
      D2SINK(i,i,:) = ZERO
   end do
#    endif
#  endif
#endif
   ! reset surface and bottom fluxes
   do i=1,NO_D3_BOX_STATES
       PELBOTTOM(i,:) = ZERO
       PELSURFACE(i,:) = ZERO 
   end do

end subroutine ResetFluxes
!EOC
!-----------------------------------------------------------------------
