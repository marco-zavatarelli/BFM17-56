#include "cppdefs.h"
#include "INCLUDE.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read forcing data, interpolate in time
! 
! NOTE: this is a template file and must be edited by the user 
!       depending on the type of input forcings
!
! !INTERFACE:
   subroutine external_forcing
!
! !DESCRIPTION:

!
! !USES:
   use global_mem, only: RLEN,ZERO
   use mem_Param,  only: p_PAR
   use constants,  only: E2W
#ifdef NOPOINTERS
   use mem
#else
   use mem,        only: ETW,ESW,EIR,SUNQ,EWIND,  &
                         EICE,ERHO,Depth
#ifdef INCLUDE_PELCO2
   use mem,        only: EPCO2air
#endif
#endif
   use time,       only: julianday, secondsofday, time_diff, &
                         julian_day,calendar_date,dayofyear
   use envforcing, only: init_forcing_vars, daylength, density, &
                         unit_forcing, read_obs, READ_ERROR, END_OF_FILE
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
!
!EOP
!
! !LOCAL VARIABLES:
! --------- to be edited by user --------------------------------------
   integer,parameter         :: NOBS=4 ! no of coulumns
   integer,parameter         :: iETW=1 ! column for temperature [degC]
   integer,parameter         :: iESW=2 ! column for salinity [-]
   integer,parameter         :: iWND=3 ! column for wind speed [m/s]
   integer,parameter         :: iEIR=4 ! column for irradiance [W/m2]
! --------- end part edited by user -----------------------------------
   integer                   :: yy,mm,dd,hh,min,ss,dyear
   real(RLEN)                :: t,alpha,jday
   real(RLEN), save          :: dt
   integer, save             :: data_jul1,data_secs1
   integer, save             :: data_jul2=0,data_secs2=0
   real(RLEN), save          :: obs1(NOBS),obs2(NOBS)=0.
   integer                   :: jh,jn
   integer                   :: ierr
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'external_forcing (jul,sec): ',julianday,secondsofday
   call  calendar_date(real(julianday,RLEN),yy,mm,dd,jh,jn)
   LEVEL2 'Calendar day:',yy,mm,dd
#endif
   ! constant sea-ice fraction (changed below if Sea-ice model is used)
   EICE(:) = ZERO

   if (init_forcing_vars) then
     data_jul2=0
     data_secs2=0
     obs2(:)=ZERO
     ! check consistency of initial date
     call read_obs(unit_forcing,yy,mm,dd,hh,min,ss,NOBS,obs2,ierr)
     call julian_day(yy,mm,dd,0,0,jday)
     if (jday > julianday )  &
        call bfm_error('external_forcing','Model start date is earlier than start date of forcing data')
     rewind(unit_forcing)
   end if
   !  This part initialise and read in new values if necessary.
   if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .lt. 0) then
      do
         data_jul1 = data_jul2
         data_secs1 = data_secs2
         obs1 = obs2
         call read_obs(unit_forcing,yy,mm,dd,hh,min,ss,NOBS,obs2,ierr)
         select case (ierr)
           case (READ_ERROR)  
              call bfm_error('external_forcing','Error reading forcing data')
           case (END_OF_FILE)  
              call bfm_error('external_forcing','Model end date is beyond forcing data')
         end select
         call julian_day(yy,mm,dd,0,0,jday)
         data_jul2 = int(jday)
         data_secs2 = hh*3600 + min*60 + ss
         if(time_diff(data_jul2,data_secs2,julianday,secondsofday) .gt. 0) EXIT
      end do
      dt = time_diff(data_jul2,data_secs2,data_jul1,data_secs1)
   end if

   !  Do the time interpolation
   t  = time_diff(julianday,secondsofday,data_jul1,data_secs1)
!-----------------------------------------------------------------------
! This is an example with required forcing variables
! Additional ones can be added updating the NOBS parameter and 
! adding lines like:
!   alpha = (obs2(iOBS)-obs1(iOBS))/dt
!   OBS(:) = obs1(OBS) + t*alpha
!-----------------------------------------------------------------------
   alpha = (obs2(iETW)-obs1(iETW))/dt
   ETW(:) = obs1(iETW) + t*alpha
   alpha = (obs2(iESW)-obs1(iESW))/dt
   ESW(:) = obs1(iESW) + t*alpha
   alpha = (obs2(iWND)-obs1(iWND))/dt
   EWIND(:) = obs1(iWND) + t*alpha
   ! Irradiance is assumed to be at the top of the box
   ! convert from irradiance (W/m2) to PAR in uE/m2/s
   alpha = (obs2(iEIR)-obs1(iEIR))/dt
   EIR(:) = (obs1(iEIR) + t*alpha)*p_PAR/E2W 

   ! Compute day length: leap years not considered (small error)
   call dayofyear(julianday,dyear)
   SUNQ=daylength(real(dyear,RLEN),latitude,ylength=365.0_RLEN)
   ! Compute density at the middle of the layer with simplified equation of state
   ERHO(:) = density(ETW(:),ESW(:),Depth(:)/2.0_RLEN)
#ifdef DEBUG
   LEVEL2 'ETW=',ETW(:)
   LEVEL2 'ESW=',ESW(:)
   LEVEL2 'EIR=',EIR(:)
   LEVEL2 'ERHO=',ERHO(:)
   LEVEL2 'EWIND=',EWIND(:)
   LEVEL2 'EICE=',EICE(:)
   LEVEL2 'SUNQ=',SUNQ
#endif
  return
   end subroutine external_forcing
!EOC

