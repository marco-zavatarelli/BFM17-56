#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read other external data, interpolate in time
!
! !INTERFACE:
   subroutine external_data
!
! !DESCRIPTION:

!
! !USES:
   use global_mem, only: RLEN,ZERO
   use constants,  only: SEC_PER_DAY
   use mem,        only: jbotR6c,jbotR6n,jbotR6p,jbotR6s, &
                         N1p,N3n,N4n,N5s,                 &
                         R6c,R6n,R6p,R6s,O2o
   use mem,        only: iiC,iiN,iiP,iiS
   use time,       only: julianday, secondsofday, time_diff, &
                         julian_day,calendar_date,dayofyear
   use envforcing
   use api_bfm
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!  modified for BFM by: Marcello Vichi
!
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
   integer,parameter         :: NOBS=4
   integer                   :: yy,mm,dd,hh,min,sec
   real(RLEN)                :: t,alpha
   real(RLEN), save          :: dt
   integer, save             :: data_jul1,data_secs1
   integer, save             :: data_jul2=0,data_secs2=0
   real(RLEN), save          :: obs1(NOBS),obs2(NOBS)=0.
   integer                   :: rc,dyear,jh,jn
   real(RLEN)                :: dfrac,jday
!-----------------------------------------------------------------------
!BOC

   if (use_external_data) then
#ifdef DEBUG
      LEVEL1 'external_data (jul,sec): ',julianday,secondsofday
      call  calendar_date(real(julianday,RLEN),yy,mm,dd,jh,jn)
      LEVEL2 'Calendar day:',yy,mm,dd
#endif
      if (init_forcing_vars) then
        data_jul2=0
        data_secs2=0
        obs2(:)=ZERO
      end if
      !  This part initialise and read in new values if necesecary.
      if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .lt. 0) then
         do
            data_jul1 = data_jul2
            data_secs1 = data_secs2
            obs1 = obs2
            call read_obs(unit_data,yy,mm,dd,hh,min,sec,NOBS,obs2,rc)
            call julian_day(yy,mm,dd,0,0,jday)
            data_jul2 = int(jday)
            data_secs2 = hh*3600 + min*60 + sec
            if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .gt. 0) EXIT
         end do
         dt = time_diff(data_jul2,data_secs2,data_jul1,data_secs1)
      end if

      !  Do the time interpolation
      t  = time_diff(julianday,secondsofday,data_jul1,data_secs1)
      alpha = (obs2(1)-obs1(1))/dt
      N1p = obs1(1) + t*alpha
      alpha = (obs2(2)-obs1(2))/dt
      N3n = obs1(2) + t*alpha
      alpha = (obs2(3)-obs1(3))/dt
      N4n = obs1(3) + t*alpha
      alpha = (obs2(4)-obs1(4))/dt
      N5s = obs1(4) + t*alpha
   end if

   ! Bottom deposition (must be negative)
   ! (mg C m^-2 d^-1 or mmol NUT m^-2 d^-1)
   ! parameters read from namelist
   ! (set to zero for no deposition)
   call dayofyear(julianday,dyear)
   dfrac = secondsofday/SEC_PER_DAY
   if (botdep_c>ZERO) jbotR6c(:) = -deposition(dyear,dfrac,botdep_c,iiC)
   if (botdep_n>ZERO) jbotR6n(:) = -deposition(dyear,dfrac,botdep_n,iiN)
   if (botdep_p>ZERO) jbotR6p(:) = -deposition(dyear,dfrac,botdep_p,iiP)
   if (botdep_si>ZERO) jbotR6s(:) = -deposition(dyear,dfrac,botdep_si,iiS)
#ifdef DEBUG
   LEVEL2 'day of year: ',dyear,dfrac
   LEVEL2 'botdep_c',botdep_c
#endif

  return
   end subroutine external_data
!EOC

