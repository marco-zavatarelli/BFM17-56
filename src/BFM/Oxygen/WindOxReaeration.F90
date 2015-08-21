#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: WindOxReaeration_3
!
! DESCRIPTION
!   Model describes reaeration between air and water column.
!       as forced by temperature and wind.
!
!       The equation and correlation used in this routine
!       are found in the 
!		R. Wanninkhof (1992), Relationship between windspeed and gas
!		exchange over the oecean
!               J. GeoPhys. Res. 97, 7373-7382
!
! !INTERFACE
  subroutine OxygenReaerationDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use constants, ONLY: HOURS_PER_DAY,ZERO_KELVIN,MW_C
  use global_mem, ONLY: RLEN,ZERO,ONE,LOGUNIT
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: O2o, cxoO2, EICE, EWIND, ETW
  use mem, ONLY: ppO2o, jsurO2o, NO_BOXES_XY, NO_BOXES,  &
     Depth, iiBen, iiPel, flux_vector 
#endif
  use mem_Param,  ONLY: AssignAirPelFluxesInBFMFlag
!  use mem_WindOxReaeration_3
#ifdef BFM_GOTM
  use bio_var, ONLY: SRFindices
#else
  use api_bfm, ONLY: SRFindices
#endif
!  use global_interface,   ONLY: CalcSchmidtNumberOx

!
! !AUTHORS
!   11 March 1998 Original version by P. Ruardij
!	              JWB 1999/03/25 Corrected k 
!
! !REVISION_HISTORY
!   2012 T. Lovato : Rearraged code and added Chemical enhancement
!   
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij and M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
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
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN) :: temp(NO_BOXES_XY)  ! deg C
  real(RLEN) :: wind(NO_BOXES_XY)  ! m/s
  real(RLEN) :: ice(NO_BOXES_XY)   ! fraction
  real(RLEN) :: tmpflux(NO_BOXES)

  real(RLEN),parameter  :: C1=1953.4_RLEN
  real(RLEN),parameter  :: C2=128.0_RLEN
  real(RLEN),parameter  :: C3=3.9918_RLEN
  real(RLEN),parameter  :: C4=0.050091_RLEN
  real(RLEN),parameter  :: O2SCHMIDT=660.0_RLEN
  real(RLEN),parameter  :: d=0.31_RLEN
  real(RLEN),parameter  :: CM2M=0.01_RLEN
  integer, save :: first=0
  real(RLEN),allocatable,save,dimension(:) :: pschmidt,temp2,bt,  &
                                         kun,O2AIRFlux,ScRatio
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BEGIN compute
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer :: AllocStatus, DeallocStatus
   if (first==0) then
      first=1
      allocate(pschmidt(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating pschmidt"
      allocate(temp2(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating temp2"
      allocate(bt(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating bt"
      allocate(kun(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating kun"
      allocate(O2AIRFlux(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating O2AIRFlux"
      allocate(ScRatio(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating ScRatio"
   end if

    temp = ETW(SRFindices)
    wind = EWIND(:)
    ice = EICE(:)
    tmpflux(:) = ZERO

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate Schmidt number,
    ! ratio between the kinematic viscosity and the molecular 
    ! diffusivity of carbon dioxide.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    temp2 = temp*temp
    pschmidt = (C1 - C2*temp + C3*temp2 - C4*temp2*temp)

    ScRatio = O2SCHMIDT / pschmidt
    !
    ! ScRatio is limited to 0 when T > 40 °C 
    WHERE(ScRatio .le. 0.0_RLEN); ScRatio=0.0_RLEN ; END WHERE
    !
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Calculate wind dependency, including conversion cm/hr => m/day :
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    kun = (d*wind*wind)*sqrt(ScRatio)* CM2M*HOURS_PER_DAY
    !
    ! This is the old formulation used before 2012 modifications
    !kun_old = (0.074E00_RLEN*wind*wind)*sqrt(O2SCHMIDT/pschmidt)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! flux o2 in mmol/m3/day   
    ! cxoO2 [mMol/m3] is the O2 concentration at saturation 
    ! computed using ETW
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !
    O2AIRFlux(:) = (ONE-ice(:)) * kun * ( cxoO2(SRFindices)- O2o(SRFindices)) 
    ! Update flux
    jsurO2o(:)  = jsurO2o(:) + O2AIRFlux(:)
    ! Convert to mmol/m2/day
    tmpflux(SRFindices) = jsurO2o(:) / Depth(SRFindices)
    if ( AssignAirPelFluxesInBFMFlag) then
        call flux_vector( iiPel, ppO2o, ppO2o, tmpflux ) 
    end if
#ifdef DEBUG
    write(LOGUNIT,*) ' Oxygen Reareation'
    write(LOGUNIT,'(A,I3,A,F12.6,A,F12.6)') 'Idx: ', SRFindices(1),' DOSat ',cxoO2(1),' DOwater ',O2o(1)
    write(LOGUNIT,'(A,F12.6,A,F10.2,A,F12.6)') 'O2 flux  ', jsurO2o(1),' Depth ',Depth(1),' kun ', kun(1)
    write(LOGUNIT,'(A,F12.6,A,F12.6,A,F12.6)') 'New Flux ', O2AIRFlux(1), ' wind ', wind(1),' temp ',temp(1)
    write(LOGUNIT,*)
#endif
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  return

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end subroutine OxygenReaerationDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
