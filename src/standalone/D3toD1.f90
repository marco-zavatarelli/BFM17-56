 function D3toD1(x,y,z)
!
! !ROUTINE: Mapping 3D variables on 1D
!           Simply returns the Z index
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO)
!
! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!
 IMPLICIT NONE
 integer,intent(IN) :: x,y,z
 integer            :: D3toD1
!
 D3toD1=z
!
end function D3toD1
