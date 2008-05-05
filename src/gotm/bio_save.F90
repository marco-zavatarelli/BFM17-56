!$Id: bio_save.F90,v 1.8 2007-03-14 12:46:07 kbk Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Storing the results
!
! !INTERFACE:
   subroutine bio_save(nlev,totn)
!
! !DESCRIPTION:
! This is a modified version of the standard bio_save routine.
! It just acts as a wrapper for the BFM saving routines.
!
! !USES:
   use bio_var
   use output, only: out_fmt
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: totn
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  Adaptation to the BFM: M. Vichi
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
!EOP
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      init_saved_vars=.false.
      first=.true.
   end if
   select case (out_fmt)
      case (ASCII)
         FATAL 'ASCII output format is not valid with BFM'
         stop 'bio_save'
      case (NETCDF)
!         call prepare_bio_output(1,nlev,_ZERO_)
         call bio_save_bfm(nlev)
      case default
         FATAL 'A non valid output format has been chosen'
         stop 'bio_save'
   end select

   return
   end subroutine bio_save
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
