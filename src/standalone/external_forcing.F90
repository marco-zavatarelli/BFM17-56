#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read forcing data, interpolate in time
!
! !INTERFACE:
   subroutine external_forcing
!
! !DESCRIPTION:

!
! !USES:
   use global_mem, only: RLEN,ZERO
   use mem,        only: ETW,ESW,EIR,SUNQ,EWIND,  &
                         EICE,ERHO,Depth
   use mem_Param,  only: p_PAR
   use constants,  only: E2W
#ifdef INCLUDE_PELCO2
   use mem,        only: EPCO2air
   use mem_CO2,    only: pco2air
#endif
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

   ! leap years not considered (small error)
   SUNQ=daylength(real(julianday,RLEN),latitude,ylength=365.0_RLEN)
   ERHO = density(ETW,ESW,Depth/2.0_RLEN)
#ifdef DEBUG
   LEVEL2 'ETW=',ETW
   LEVEL2 'ESW=',ESW
   LEVEL2 'EIR=',EIR
   LEVEL2 'ERHO=',ERHO
   LEVEL2 'EWIND=',EWIND
   LEVEL2 'EICE=',EICE
   LEVEL2 'SUNQ=',ETW
#endif
  return
   end subroutine external_forcing
!EOC
