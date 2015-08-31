!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Get the numeric timestep
!
! !INTERFACE:
   function GetDelta() result(Delta)
!
! !DESCRIPTION:
!  Transfer the integration time step to the BFM
!  Unit conversion from seconds to days
!
! !USES:
   use constants, only: RLEN, SEC_PER_DAY
   use time,      only: timestep
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
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
!EOP
!-----------------------------------------------------------------------
!BOC
   real(RLEN) :: Delta

   Delta = timestep/SEC_PER_DAY

   return
   end function GetDelta
!EOC
!-----------------------------------------------------------------------


