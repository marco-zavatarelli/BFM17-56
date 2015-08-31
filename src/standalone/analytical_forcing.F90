#include "cppdefs.h"
#include "INCLUDE.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: analytical_forcing
!
! !INTERFACE
subroutine analytical_forcing
!
! !DESCRIPTION
!   Define analytical forcings used in the BFM
! !USES
   use api_bfm
   use global_mem, only: RLEN,ONE
#ifdef NOPOINTERS
   use mem
#else
   use mem,        only: ETW,ESW,EIR,SUNQ,ThereIsLight,EWIND,  &
                         EICE,jbotR6c,jbotR6n,jbotR6p,jbotR6s, &
                         R6c,R6n,R6p,R6s,O2o,ERHO,Depth,xEPS
   use mem,        only: iiC,iiN,iiP,iiS
   use mem_param,  only: p_small
#ifdef INCLUDE_SEAICE
   ! seaice forcings
   use mem,        only: EVB,ETB,ESB,EIB,EHB,ESI
#endif
#endif
   use mem_PAR,    only: ChlAttenFlag, P_PARRGB, P_PAR, &
                         R_EPS, B_EPS, G_EPS,           &
                         EIRR, EIRB, EIRG
   use constants,  only: E2W, SEC_PER_DAY
   use standalone, only: timesec,latitude
   use envforcing
#ifdef INCLUDE_PELCO2
   use mem_CO2,    only: AtmCO20, AtmCO2, AtmSLP, AtmTDP
#endif
   use time,       only: julianday, secondsofday, timefmt, &
                         julian_day,calendar_date,dayofyear
   use SystemForcing, ONLY :FieldRead
   use bfm_error_msg

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO), M. Vichi (CMCC)
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
! !LOCAL VARIABLES:
   real(RLEN)          :: dfrac,wlight,dtime
   integer             :: dyear
   real(RLEN),external :: GetDelta
   real(RLEN)          :: biodelta
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'envforcing_bfm: analytical'
   LEVEL2 'time=',timesec
#endif
   !---------------------------------------------
   ! Computes all the forcings
   !---------------------------------------------
   dtime = timesec/SEC_PER_DAY
   if (timefmt==2) then
      call dayofyear(julianday,dyear)
      dfrac = secondsofday/SEC_PER_DAY
   else
      dfrac=(dtime-floor(dtime)) ! fraction of the day
      dyear=mod(dtime,360._RLEN) ! Day of the year
   end if
   SUNQ(:)=daylength(REAL(dyear,RLEN),latitude)
   wlight=light(dyear,dfrac)
   select case(ltype)
    case (3) ! light on/off distribution for daylight average
      ThereIsLight(:)=lightAtTime(dfrac,SUNQ(1))
      wlight=wlight*ThereIsLight(1)
    case (4) ! light scaled up by photoperiod
      wlight=wlight*24.0_RLEN/SUNQ(1)
    case (1) ! instantaneous light distribution
      wlight=instLight(wlight,SUNQ(1),dfrac)
    case default ! light constant during the day
   end select
   ETW(:) = temperature(dyear,dfrac)
   ESW(:) = salinity(dyear,dfrac)
   ! compute density at the middle of the layer
   ERHO(:) = density(ETW(:),ESW(:),Depth(:)/2.0_RLEN)

   ! convert from irradiance to PAR in uE/m2/s
   select case (ChlAttenFlag)
      case (1) ! Broadband
         EIR(:) = wlight*p_PAR/E2W
      case (2) ! RGB
         EIR(:) = p_PARRGB*(wlight+p_small)/E2W
         EIRR(:) = EIR(:) * exp ( -R_eps(:) )
         EIRG(:) = EIR(:) * exp ( -G_eps(:) )
         EIRB(:) = EIR(:) * exp ( -B_eps(:) )
         EIR(:) = EIRB(:) + EIRG(:) + EIRR(:)
         ! weighted broadband diffuse attenuation coefficient for diagnostics
         xEPS(:) = (EIRB(:)*B_eps(:) + EIRG(:)*G_eps(:) + EIRR(:)*R_eps(:))/EIR(:)
      case default 
         call BFM_ERROR("analytical_forcing","Bad value for ChlAttenFlag.")
   end select
      
   ! analytical wind velocity 
   EWIND(:) = wind(dyear,dfrac)
   ! constant sea-ice fraction 
   EICE(:) = 0.0_RLEN
#ifdef INCLUDE_PELCO2
   if (AtmCO2%init == 0) then
      ! increase of initial CO2 concentration in the air by CO2inc [%]
      AtmCO2%fnow = AtmCO20*(ONE + (CO2inc / 100_RLEN) / 365._RLEN * dtime)
   else
      call FieldRead(AtmCO2)
   endif
!  Update sea level pressure
   if ( allocated(AtmSLP%fnow)) CALL FieldRead(AtmSLP)

!  Update Air Dew Point temperature
   if ( allocated(AtmTDP%fnow)) CALL FieldRead(AtmTDP)
#endif
#ifdef DEBUG
   LEVEL2 'ETW=',ETW(:)
   LEVEL2 'ESW=',ESW(:)
   LEVEL2 'EIR=',EIR(:)
   LEVEL2 'ERHO=',ERHO(:)
   LEVEL2 'EWIND=',EWIND(:)
   LEVEL2 'EICE=',EICE(:)
#endif

   ! Bottom deposition (must be negative)
   ! (mg C m^-2 d^-1 or mmol NUT m^-2 d^-1)
   ! parameters read from namelist
   ! (set to zero for no deposition)
   if (botdep_c>ZERO) jbotR6c(:) = -deposition(dyear,dfrac,botdep_c,iiC)
   if (botdep_n>ZERO) jbotR6n(:) = -deposition(dyear,dfrac,botdep_n,iiN)
   if (botdep_p>ZERO) jbotR6p(:) = -deposition(dyear,dfrac,botdep_p,iiP)
   if (botdep_si>ZERO) jbotR6s(:) = -deposition(dyear,dfrac,botdep_si,iiS)

#ifdef INCLUDE_SEAICE
! sea-ice environmental forcings
! overwritten by the seaice data file, if present
  EVB(:) = ONE
  ETB(:) = ONE
  ESB(:) = ONE
  ! convert from irradiance to PAR in uE/m2/s
  EIB(:) = ONE/E2W
  EHB(:) = ONE
  ESI(:) = ONE
#endif

end subroutine analytical_forcing
!EOC
!-----------------------------------------------------------------------

