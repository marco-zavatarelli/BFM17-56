#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read other external event data (instantaneous events)
!
! !INTERFACE:
   subroutine event_data
!
! !DESCRIPTION:
!  This routine reads values of model state variables at a certain 
!  moment in time and assign them to the STATE array
!  (thus overriding any value computed during the integration)
!  This routine is primarily meant for inoculation experiments. 
!  It now reads phytoplankton concentrations but the user needs 
!  to adapt the reading and assignements to its own need.
!
! !USES:
   use global_mem, only: RLEN,ZERO
   use constants,  only: SEC_PER_DAY
   use mem
   use time,       only: julianday, secondsofday, time_diff, &
                         julian_day,calendar_date,dayofyear
   use envforcing
   use api_bfm
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!                      using functions from Karsten Bolding
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
!EOP
!
! !LOCAL VARIABLES:
   ! The number of event variables is now equal to the number 
   ! of Phyto groups. Additional variables can be added
   integer,parameter         :: NOBS=iiPhytoPlankton

   integer                   :: yy,mm,dd,hh,min,sec
   real(RLEN)                :: t,alpha,jday
   real(RLEN), save          :: dt
   integer, save             :: data_jul1,data_secs1
   integer, save             :: data_jul2=0,data_secs2=0
   real(RLEN), save          :: obs1(NOBS),obs2(NOBS)=0.
   integer                   :: rc,i,jh,jn
   ! Redfield and other constant ratios
   real(RLEN),parameter      :: nc_ratio=0.0126_RLEN 
   real(RLEN),parameter      :: pc_ratio=0.7862e-3_RLEN
   real(RLEN),parameter      :: sc_ratio=0.0145_RLEN
   real(RLEN),parameter      :: lc_ratio=0.03_RLEN
!-----------------------------------------------------------------------
!BOC

   if (use_event_data) then
#ifdef DEBUG
      LEVEL1 'event_data (jul,sec): ',julianday,secondsofday
      call  calendar_date(real(julianday,RLEN),yy,mm,dd,jh,jn)
      LEVEL2 'Calendar day:',yy,mm,dd
#endif
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
            call read_obs(unit_event,yy,mm,dd,hh,min,sec,NOBS,obs2,rc)
            call julian_day(yy,mm,dd,0,0,jday)
            data_jul2 = int(jday)
            data_secs2 = hh*3600 + min*60 + sec
            if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .gt. 0) EXIT
         end do
         dt = time_diff(data_jul2,data_secs2,data_jul1,data_secs1)
         !  Assign the instantaneous value (ex. phytoplankton)
         do i=1,iiPhytoPlankton
            if (obs2(i) > ZERO) then
               D3STATE(ppPhytoPlankton(i,iiC),:) = obs2(i) 
               D3STATE(ppPhytoPlankton(i,iiN),:) = obs2(i) * nc_ratio
               D3STATE(ppPhytoPlankton(i,iiP),:) = obs2(i) * pc_ratio
               D3STATE(ppPhytoPlankton(i,iiS),:) = obs2(i) * sc_ratio
               D3STATE(ppPhytoPlankton(i,iiL),:) = obs2(i) * lc_ratio
            end if
         end do
      end if

   end if

  return
   end subroutine event_data
!EOC

