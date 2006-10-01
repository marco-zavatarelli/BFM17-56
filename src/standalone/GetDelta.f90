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
   use time, only: timestep
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
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


