!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: calc_sigma_depth
!
! !INTERFACE:
   subroutine calc_sigma_depth(nlev,ddu,maxdepth,arr)
!
! !DESCRIPTION:
!  This routine defines the depth distribution
!  of vertical levels for the output of nutrient profiles 
!  in the benthic nutrient model.
!  This routine is a simplification of
!  the calculation used in gotm/getm
!
! !USES:
   use global_mem, ONLY: RLEN,ONE,ZERO
   IMPLICIT NONE
! !INPUT PARAMETERS:
     integer,intent(IN)           :: nlev
     real(RLEN),intent(IN)        :: ddu
     real(RLEN),intent(IN)        :: maxdepth
! !OUTPUT PARAMETERS:
     real(RLEN),intent(OUT)       :: arr(1:nlev)
! !LOCAL PARAMETERS:
     real(RLEN)                   :: r,s
     integer                      :: i
!
!
! !REVISION HISTORY:
!  Original by Piet Ruardij
!
!EOP
!-----------------------------------------------------------------------
!BOC

   r =ZERO
   do i=1,nlev
      s= maxdepth*(ONE - tanh(ddu*float(nlev-i)/float(nlev))/tanh(ddu))
      arr(i)=(s+r) * 0.5_RLEN
      r=s
   end do
   return
   end subroutine calc_sigma_depth
!EOC
!-----------------------------------------------------------------------
