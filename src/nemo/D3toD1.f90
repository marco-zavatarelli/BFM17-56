!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Specialized routine to couple BFM with NEMO
!
! !INTERFACE:
   function D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ) result(BoxNumber)
!
! !DESCRIPTION:
!  This routine resolve the mapping between BFM 1D structure 
!  and 3D OGCMs (volume fields)
!
! !USES:
   use mem, ONLY: NO_BOXES_X,NO_BOXES_Y
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   integer :: BoxNumberX,BoxNumberY,BoxNumberZ
   integer :: BoxNumber

   BoxNumber = (BoxNumberZ-1) * (NO_BOXES_Y*NO_BOXES_X) +  &
               (BoxNumberY+NO_BOXES_Y*(BoxNumberX-1)) 

   return
   end function D3toD1
!EOC
!-----------------------------------------------------------------------


