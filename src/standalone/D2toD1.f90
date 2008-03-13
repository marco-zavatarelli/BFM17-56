function D2toD1(x,y)
!
! Mapping 2D variables on 1D
! Dummy routine for the standalone model
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Momme Butenschoen, March 2006 !
! Dipartimento di Fisica        !
! Universita' di Bologna        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
 use mem, only: NO_BOXES,NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z
 IMPLICIT NONE
 integer,intent(IN) :: x,y
 integer            :: D2toD1
!
 D2toD1=1
!
 end function D2toD1
