#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read seaice data, interpolate in time
!
! !INTERFACE:
   subroutine external_seaice
!
! !DESCRIPTION:

!
! !USES:
   use global_mem, only: RLEN,ZERO
   use mem_Param,  only: p_PAR
   use constants,  only: E2W
   ! seaice forcings
   use mem,        only: EVB,ETB,ESB,EIB,EHB,ESI
   use time, only: julianday, secondsofday, time_diff, &
                   julian_day,calendar_date
   use envforcing, only: init_forcing_vars, daylength, density, &
                         unit_forcing, read_obs
   use standalone, only: latitude
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!  modified for BFM by: Marcello Vichi
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   real(RLEN)                :: t,alpha
   real(RLEN), save          :: dt
   integer, save             :: data_jul1,data_secs1
   integer, save             :: data_jul2=0,data_secs2=0
   real(RLEN), save          :: obs1(4),obs2(4)=0.
   integer                   :: rc
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'external_forcing (jul,sec): ',julianday,secondsofday
   call  calendar_date(julianday,yy,mm,dd)
   LEVEL2 'Calendar day:',yy,mm,dd
#endif
   ! constant sea-ice fraction (changed below if Sea-ice model is used)
   EICE = ZERO

   if (init_forcing_vars) then
     data_jul2=0
     data_secs2=0
     obs2(:)=ZERO
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .lt. 0) then
      do
         data_jul1 = data_jul2
         data_secs1 = data_secs2
         obs1 = obs2
         call read_obs(unit_forcing,yy,mm,dd,hh,min,ss,4,obs2,rc)
         call julian_day(yy,mm,dd,data_jul2)
         data_secs2 = hh*3600 + min*60 + ss
         if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .gt. 0) EXIT
      end do
      dt = time_diff(data_jul2,data_secs2,data_jul1,data_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(julianday,secondsofday,data_jul1,data_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   ETW = obs1(1) + t*alpha
   alpha = (obs2(2)-obs1(2))/dt
   ESW = obs1(2) + t*alpha
   alpha = (obs2(3)-obs1(3))/dt
   ! convert from irradiance (W/m2) to PAR in uE/m2/s
   EIR = (obs1(3) + t*alpha)*p_PAR/E2W
   alpha = (obs2(4)-obs1(4))/dt
   EWIND = obs1(4) + t*alpha

! sea-ice environmental forcings
  EICE(:) = obs(1)
  EVB(:) = obs(2)
  ETB(:) = obs(3
  ESB(:) = obs(4)
  ! convert from irradiance to PAR in uE/m2/s
  EIB(:) = obs(5)/E2W
  EHB(:) = obs(6)
  ESI(:) = obs(7)

#ifdef DEBUG
   LEVEL2 'EICE=',EICE
#endif
  return
   end subroutine external_forcing
!EOC

