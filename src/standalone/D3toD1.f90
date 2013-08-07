 function D3toD1(x,y,z)
!
! Mapping 3D variables on 1D
! Simply returns the Z index
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Momme Butenschoen, March 2006 !
! Dipartimento di Fisica        !
! Universita' di Bologna        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
 IMPLICIT NONE
 integer,intent(IN) :: x,y,z
 integer            :: D3toD1
!
 D3toD1=z
!
end function D3toD1
