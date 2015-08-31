#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read seaice data, interpolate in time
!
! !INTERFACE:
   subroutine external_seaice
#ifdef INCLUDE_SEAICE
!
! !DESCRIPTION:
! Coupling of sea ice processes for Kobbefjord site
!
! !USES:
   use global_mem, only: RLEN,ZERO,ONE
   use constants,  only: E2W, SEC_PER_DAY
   ! seaice forcings
   use mem,        only: EICE,EVB,ETB,ESB,EIB,EHB,ESI,EDH,EDS
   use mem,        only: ESW,EIR
   use time,       only: julianday, secondsofday, time_diff, &
                         julian_day,calendar_date
   use envforcing, only: init_forcing_vars, daylength, density, &
                         unit_seaice, read_obs, END_OF_FILE, READ_ERROR
   use standalone, only: latitude
   use bfm_error_msg
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!  modified for BFM by: Marcello Vichi
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
!
!EOP
!
! !LOCAL VARIABLES:
   integer,parameter                  :: NSI=9
   integer                            :: yy,mm,dd,hh,minutes,ss
   real(RLEN)                         :: t,alpha,jday
   real(RLEN), save                   :: dt
   integer, save                      :: data_jul1,data_secs1
   integer, save                      :: data_jul2=0,data_secs2=0
   real(RLEN), save                   :: obs1(NSI),obs2(NSI)=0.
   integer                            :: ierr,jh,jn
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'external_forcing (jul,sec): ',julianday,secondsofday
   call  calendar_date(real(julianday,RLEN),yy,mm,dd,jh,jn)
   LEVEL2 'Calendar day:',yy,mm,dd
#endif

   ! constant sea-ice fraction (changed below if Sea-ice model is used)
   EICE = ZERO

   if (init_forcing_vars) then
     data_jul2=0
     data_secs2=0
     obs2(:)=ZERO
     ! check consistency of initial date
     call read_obs(unit_seaice,yy,mm,dd,hh,minutes,ss,NSI,obs2,ierr)
     call julian_day(yy,mm,dd,0,0,jday)
     if (jday > julianday )  &
        call bfm_error('external_seaice','Model start date is earlier than start date of sea ice data')
     rewind(unit_seaice)
   end if
!  This part initialise and read in new values if necessary.
   if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .lt. 0) then
      do
         data_jul1 = data_jul2
         data_secs1 = data_secs2
         obs1 = obs2
         call read_obs(unit_seaice,yy,mm,dd,hh,minutes,ss,NSI,obs2,ierr)
         select case (ierr)
           case (READ_ERROR)
              call bfm_error('external_forcing','Error reading forcing data')
           case (END_OF_FILE)
              call bfm_error('external_forcing','Model end date is beyond forcing data')
         end select
         call julian_day(yy,mm,dd,0,0,jday)
         data_jul2 = int(jday)
         data_secs2 = hh*3600 + minutes*60 + ss
         if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .gt. 0) EXIT
      end do
      dt = time_diff(data_jul2,data_secs2,data_jul1,data_secs1)
   end if

!  Do the time interpolation
   t  = time_diff(julianday,secondsofday,data_jul1,data_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   EICE(:) = obs1(1) + t*alpha
   alpha = (obs2(2)-obs1(2))/dt
   EHB(:) = obs1(2) + t*alpha
   alpha = (obs2(3)-obs1(3))/dt
   EVB(:) = obs1(3) + t*alpha
   alpha = (obs2(4)-obs1(4))/dt
   ETB(:) = obs1(4) + t*alpha
   alpha = (obs2(5)-obs1(5))/dt
   ESB(:) = obs1(5) + t*alpha
   ! convert from irradiance (W/m2) to PAR in uE/m2/s for BFM to run 
   alpha = (obs2(6)-obs1(6))/dt
  ! the PAR/IRRAD is 0.5 for ocean, but for sea ice it is computed by 
  ! the physical model so it is not used
   EIB(:) = ((obs1(6) + t*alpha)/E2W)
   alpha = (obs2(7)-obs1(7))/dt
   ESI(:) = obs1(7) + t*alpha
   alpha = (obs2(8)-obs1(8))/dt
   EDH(:)=(obs1(8) + t*alpha)*SEC_PER_DAY
   alpha = (obs2(9)-obs1(9))/dt
   EDS(:)=(obs1(9) + t*alpha)*SEC_PER_DAY

#ifdef DEBUG
   LEVEL2 'EICE=',EICE
   LEVEL2 'EHB=',EHB
   LEVEL2 'EVB=',EVB
   LEVEL2 'ETB=',ETB
   LEVEL2 'ESB=',ESB
   LEVEL2 'EIB=',EIB
   LEVEL2 'ESI=',ESI
   LEVEL2 'EDH=',EDH
   LEVEL2 'EDS=',EDS
#endif
  return
#endif
   end subroutine external_seaice
!EOC

