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
! Coupling of sea ice processes for Kobbefjord site
!
! !USES:
   use global_mem, only: RLEN,ZERO,ONE
   use mem_Param,  only: p_PAR,CalcSeaiceAlgae,CalcPhytoPlankton
   use mem_Param,  only: CalcSeaiceZoo,CalcSeaiceBacteria 
   use mem_Param,  only: p_eps0,p_epsR6,p_epsChla, p_small
   use constants,  only: E2W, SEC_PER_DAY
   ! seaice forcings
   use mem,        only: EICE,EVB,ETB,ESB,EIB,EHB,ESI,EDH,EDS,F3c, &
                         F2o,I1p,I3n,I4n,I5s,S1l,S2l,U6c
   use mem,        only: iiBen,ppI1p,N1p,N3n,N4n,N5s,O2o,O3c,P1l,P2l,NO_BOXES_XY, &
                         flux_vector, ppI3n,ppI4n,ppI5s,ppF2o,ppF3c,&
                         ppS1l,ppS2l
   use mem,        only: SeaiceAlgae,ppSeaiceAlgae,PhytoPlankton,ppPhytoPlankton,PELSURFACE, &
                         iiS1,iiS2,iiP1,iiP2,iiSeaiceAlgae,iiPhytoPlankton,iiPel
   use mem,        only: SeaiceDetritus,ppSeaiceDetritus,PelDetritus,ppPelDetritus, &
                         iiU1,iiU6,iiR1,iiR6
   use mem,        only: SeaiceBacteria,ppSeaiceBacteria,PelBacteria,ppPelBacteria, &
                         iiT1,iiB1
   use mem,        only: SeaiceZoo,ppSeaiceZoo,MicroZooPlankton,ppMicroZooPlankton, &
                         iiX1,iiZ5
   use mem,        only: iiC,iiN,iiP,iiS,iiL
   use mem,        only: jsurN1p,jsurN3n,jsurN4n,jsurN5s,jsurO3c, &
                         jsurO2o,jsurP1l,jsurP2l,Depth,NO_BOXES,ppN1p,ppN3n,ppN4n,ppN5s,ppO2o
   use mem,        only: ESW,EIR
   use time,       only: julianday, secondsofday, time_diff, &
                         julian_day,calendar_date
   use envforcing, only: init_forcing_vars, daylength, density, &
                         unit_seaice, read_obs
   use standalone, only: latitude
   use api_bfm,    only: SRFindices
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
   integer,parameter                  :: NSI=9
   integer                            :: yy,mm,dd,hh,minutes,ss
   real(RLEN)                         :: t,alpha
   real(RLEN), save                   :: dt
   integer, save                      :: data_jul1,data_secs1
   integer, save                      :: data_jul2=0,data_secs2=0
   real(RLEN), save                   :: obs1(NSI),obs2(NSI)=0.
   integer                            :: rc
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_N1,flux_pel_ice_N3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_N4,flux_pel_ice_N5
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_O2,flux_pel_ice_O3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice, flux_atm_N1, flux_atm_N3
   real(RLEN), dimension(NO_BOXES_XY) :: I1p_tilde,I3n_tilde,I4n_tilde,I5s_tilde
   real(RLEN), dimension(NO_BOXES_XY) :: F2o_tilde,F3c_tilde
   real(RLEN), dimension(NO_BOXES_XY) :: I1p_star,I3n_star,I4n_star,I5s_star
   real(RLEN), dimension(NO_BOXES_XY) :: F2o_star,F3c_star

   integer                            :: i,j,p
   real(RLEN), dimension(:), pointer  :: lcl_PelagicVar,lcl_SeaiceVar
   real(RLEN)                         :: tmpflux(NO_BOXES)
   real(RLEN)                         :: localdelta
   real(RLEN), external               :: GetDelta
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
         call read_obs(unit_seaice,yy,mm,dd,hh,minutes,ss,NSI,obs2,rc)
         call julian_day(yy,mm,dd,data_jul2)
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
  ! the PAR/IRRAD is 0.5 for ocean, but for sea ice I do compute it from the physical model so it is not needed!
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

   end subroutine external_seaice
!EOC

