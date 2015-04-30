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
   use oce_trc,   only: rdt   !: time step for the dynamics
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
!EOP
!-----------------------------------------------------------------------
!BOC
   real(RLEN) :: Delta

   Delta = rdt/SEC_PER_DAY

   return
   end function GetDelta
!EOC
!-----------------------------------------------------------------------


