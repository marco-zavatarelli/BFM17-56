!$Id: analytical_profile.F90,v 1.6 2007-01-06 11:49:15 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: analytical_profile
!
! !INTERFACE:
   subroutine analytical_profile(nlev,z,z1,v1,z2,v2,prof)
!
! !DESCRIPTION:
! This routine creates a vertical profile {\tt prof} with value
! {\tt v1} in a surface layer down to depth {\tt z1} and a bottom
! layer of value {\tt v2} reaching from depth {\tt z2} down to the bottom.
! Both layers are connected by an intermediate layer reaching from {\tt z1}
! to {\tt z2} with values linearly varying from {\tt v1} to {\tt v2}.
!
! !USES:
   use global_mem, only:RLEN,ZERO,bfm_lwp,LOGUNIT
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   real(RLEN), intent(in)                :: z(nlev)
   real(RLEN), intent(in)                :: z1,v1,z2,v2
!
! !OUTPUT PARAMETERS:
   real(RLEN), intent(out)               :: prof(nlev)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   real(RLEN)                  :: alpha
!
!-----------------------------------------------------------------------
!BOC
   if (z2-z1 .gt. -1.e-15) then
         alpha = (v2-v1)/(z2-z1+2.e-15)
   else
      STDERR '**********************************************'
      STDERR '* Error detected by analytical_profile.F90:  *'
      STDERR '*  anz2 should be larger than anz1.          *'
      STDERR '*   Please edit BFM_General.nml or bio_bfm.nml.      *'
      STDERR '**********************************************'
      stop
   end if

   do i=nlev,1,-1
!PROVA      if(-1.*z(i) .le. z1) then
      if(z(i) .le. z1) then
         prof(i) = v1
      end if
      if (alpha.le.1.e15) then
!PROVA         if(-1.*z(i) .gt. z1 .and. -1.*z(i) .le. z2) then
         if(z(i) .gt. z1 .and. z(i) .le. z2) then
            prof(i) = v1 + alpha*(z(i)-z1)
         end if
      end if
!PROVA      if(-1.*z(i) .gt. z2) then
      if(z(i) .gt. z2) then
         prof(i) = v2
      end if
   end do

   return
   end subroutine analytical_profile
!EOC

!-----------------------------------------------------------------------
! Copyright 2013 BFM System Team (bfm_st@lists.cmcc.it)
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
