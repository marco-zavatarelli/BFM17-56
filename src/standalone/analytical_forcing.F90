#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcings used in the BFM
!
! !INTERFACE
   subroutine analytical_forcing
!
! !DESCRIPTION
!
! !USES
   use api_bfm
   use global_mem, only: RLEN
   use mem,        only: ETW,ESW,EIR,SUNQ,ThereIsLight,EWIND,  &
                         EICE,jbotR6c,jbotR6n,jbotR6p,jbotR6s, &
                         R6c,R6n,R6p,R6s,O2o,ERHO,Depth
   use mem,        only: iiC,iiN,iiP,iiS
   use mem_Param,  only: LightForcingFlag,p_PAR
   use constants,  only: E2W, SEC_PER_DAY
   use standalone, only: timesec,latitude
   use envforcing
#ifdef INCLUDE_PELCO2
   use mem,        only: EPCO2air
   use mem_CO2,    only: pco2air
#endif
#ifdef INCLUDE_SEAICE
   ! seaice forcings
   use mem,        only: EVB,ETB,ESB,EIB,EHB,ESI
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO)
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
   SUNQ=daylength(dtime,latitude)
   dfrac=(dtime-floor(dtime)) ! fraction of the day
   dyear=mod(dtime,360._RLEN) ! Day of the year
   wlight=light(dyear,dfrac)
   select case(LightForcingFlag)
    case (3) ! light on/off distribution for daylight average
      ThereIsLight=lightAtTime(dfrac,sunq)
      wlight=wlight*ThereIsLight
    case (1) ! instantaneous light distribution
      wlight=instLight(wlight,sunq,dfrac)
    case default ! light constant during the day
   end select
   ETW = temperature(dyear,dfrac)
   ESW = salinity(dyear,dfrac)
   ERHO = density(ETW,ESW,Depth/2.0_RLEN)
   ! convert from irradiance to PAR in uE/m2/s
   EIR = wlight*p_PAR/E2W
   ! analytical wind velocity 
   EWIND = wind(dyear,dfrac)
   ! constant sea-ice fraction 
   EICE = 0.0_RLEN
#ifdef INCLUDE_PELCO2
   ! 1% increase of pCO2 in the air 
   EPCO2air = pco2air*(ONE+0.01_RLEN/360._RLEN*dtime)
#endif
#ifdef DEBUG
   LEVEL2 'ETW=',ETW
   LEVEL2 'ESW=',ESW
   LEVEL2 'EIR=',EIR
   LEVEL2 'ERHO=',ERHO
   LEVEL2 'EWIND=',EWIND
   LEVEL2 'EICE=',EICE
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
! Reading from file to be added
  EVB = ONE
  ETB = ONE
  ESB = ONE
  ! convert from irradiance to PAR in uE/m2/s
  EIB = ONE/E2W
  EHB = ONE
  ESI = ONE
#endif

   end subroutine analytical_forcing
!EOC
!-----------------------------------------------------------------------

