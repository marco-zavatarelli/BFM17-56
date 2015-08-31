function D2toD1(x,y)
!
! ROUTINE: Mapping 2D variables on 1D
!          Dummy routine for the standalone model
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
 IMPLICIT NONE
 integer,intent(IN) :: x,y
 integer            :: D2toD1
!
 D2toD1=1
!
 end function D2toD1
