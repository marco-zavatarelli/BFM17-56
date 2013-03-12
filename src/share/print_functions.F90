#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: print_functions --- 
!
! !INTERFACE:
   MODULE print_functions
!
! !DESCRIPTION:
!   Print functions to output arrays
! !USES:
   use global_mem, ONLY: RLEN,ZERO,ONE
   IMPLICIT NONE

   public prxy
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: the famous PRXY from glorious POM
!
! !INTERFACE:
     subroutine prxy(theunit,label,array, &
                     im,iskp,jm,jskp,scala)
!
! !DESCRIPTION:
!
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Writes a horizontal 2-D field.                      *
! *                                                                    *
! *                label ....... label for output                      *
! *                array(im,jm). array to be printed                   *
! *                iskp ........ skipping interval for i               *
! *                jskp ........ skipping interval for j               *
! *                scala ....... < 0 for floating point numbers output *
! *                              0 for integer output, divisor for a   *
! *                                based on magnitudes of |a| values   *
! *                              > 0 for integer output, divisor for a *
! *                                given by scala                      *
! *                                                                    *
! **********************************************************************
!
! !INPUT PARAMETERS:
     integer, intent(in)                      :: theunit, im, iskp, &
                                                 jm, jskp
     character(len=*), intent(in)             :: label
     real(RLEN), dimension(im,jm), intent(in) :: &
         array                       ! integer 2D array to be print
     real(RLEN), intent(in)                   :: scala
!
! !REVISION HISTORY:
!  Original author(s): G. Mellor
!
!EOP
!
! !LOCAL VARIABLES:
    real(RLEN)                             :: amx,scale
    integer                                :: i,ib,ie,j,jwr,cols
!-------------------------------------------------------------------------
!BOC
      if(scala.ge.ZERO) then
        cols=24
      else
        cols=12
      endif

      if (scala.lt.ZERO) scale = ONE
      if (scala.eq.ZERO) then
        amx=1.e-12_RLEN
        do j=1,jm,jskp
          do i=1,im,iskp
            amx=max(abs(array(i,j)),amx)
          end do
        end do
          scale=10.e0_RLEN**(int(log10(amx)+100.e0_RLEN)-103)
        endif
      if(scala.gt.ZERO) scale=scala

      write(theunit,1) label
    1 format(1x,a40/)
      write(theunit,2) scale
    2 format('  multiply all values by ',1f8.2)

      do ib=1,im,cols*iskp

        ie=ib+(cols-1)*iskp
        if(ie.gt.im) ie=im

        if(scala.ge.ZERO) then
          write(theunit,3) (i,i=ib,ie,iskp)
    3     format(/,2x,24i5,/)
        else
          write(theunit,4) (i,i=ib,ie,iskp)
    4     format(/,12i10,/)
        endif

        do j=1,jm,jskp
          jwr=jm+1-j
          if(scala.ge.ZERO) then
            write(theunit,5) jwr,(nint(array(i,jwr)/scale),i=ib,ie,iskp)
    5       format(1x,i3,24i5)
          else
            write(theunit,6) jwr,(array(i,jwr),i=ib,ie,iskp)
    6       format(1x,i2,12(e10.2))
          endif
        end do

        write(theunit,7)
    7   format(//)

      end do

      end subroutine prxy
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------

end module print_functions

!-----------------------------------------------------------------------
