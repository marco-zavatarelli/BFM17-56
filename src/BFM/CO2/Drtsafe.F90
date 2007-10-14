#include "DEBUG.h"
#include "INCLUDE.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Drtsafe
!
! !INTERFACE:
function drtsafe(x1,x2,xacc)
!
! !DESCRIPTION:
!   find roots of the Total Alkalinity function ta_iter_1 
!   by Newton-Raphson and bisection
!   Adapted and optimized from Numerical Recipes rtsafe.f90 
!   (error checking have been removed)
!
! !USES:
   use mem_co2
   use global_mem
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(IN)  :: x1,x2,xacc
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   real(RLEN)             :: drtsafe
!
! !REVISION HISTORY:
!  Author(s):
!   WH Press, SA Teukolsky, WT Vetterling and BP Flannery
!   Adapted from OCMIP standard files by M. Vichi
!
! !LOCAL VARIABLES:
   real(RLEN)             :: swap  
   real(RLEN)             :: df,dx,dxold,f,fh,fl,temp,xh,xl
   integer                :: j
   integer,parameter      :: MAXIT=100
!
!EOP
!-----------------------------------------------------------------------
!BOC


   call Ta_Iter_1(x1,fl,df)
   call Ta_Iter_1(x2,fh,df)

   if (fl == ZERO) then
     drtsafe=x1
     return
   else if (fh == ZERO) then
     drtsafe=x2
     return
   else if (fl < ZERO) then
     xl=x1
     xh=x2
   else
     xh=x1
     xl=x2
     swap=fl
     fl=fh
     fh=swap
   end if

   drtsafe=0.5_RLEN*(x1+x2)
   dxold=abs(x2-x1)
   dx=dxold
   call Ta_Iter_1(drtsafe,f,df)

do j=1,MAXIT
  if (((drtsafe-xh)*df-f)*((drtsafe-xl)*df-f) >= ZERO .OR. &
       abs(2.0_RLEN*f) > abs(dxold*df) ) then
     dxold=dx
     dx=0.5_RLEN*(xh-xl)
     drtsafe=xl+dx
     if (xl == drtsafe) return
  else
     dxold=dx
     dx=f/df
     temp=drtsafe
     drtsafe=drtsafe-dx
     if (temp == drtsafe) return
  end if
  if (abs(dx) < xacc) return
  call Ta_Iter_1(drtsafe,f,df)
  if (f < ZERO) then
     xl=drtsafe
     fl=f
  else
     xh=drtsafe
     fh=f
  end if
end do

!-----------------------------------------------------------------------
end function drtsafe
!-----------------------------------------------------------------------
!EOC
