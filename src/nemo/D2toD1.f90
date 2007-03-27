!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Specialized routine to couple BFM with NEMO
!
! !INTERFACE:
   function D2toD1(BoxNumberX,BoxNumberY) result(BoxNumber)
!
! !DESCRIPTION:
!  This routine resolve the mapping between BFM 1D structure 
!  and 3D OGCMs (surface fields)
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   integer :: BoxNumberX,BoxNumberY
   integer :: BoxNumber

   BoxNumber = 0

   return
   end function D2toD1
!EOC
!-----------------------------------------------------------------------


