#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenOxygen
!
! DESCRIPTION
!   Description of first order oxic processes in the sediment
!       and computation of oxygen penetration depth
!
! !INTERFACE
  subroutine BenOxygenDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: D1m, G2o, D2m
  use mem, ONLY: ppD1m, ppG2o, ppD2m, InitializeModel, LocalDelta, shiftD1m, &
    jbotO2o, ETW_Ben, irrenh, rrBTo, jG2K3o, jG2K7o, O2o_Ben, NO_BOXES_XY, iiBen, &
    iiPel, flux_vector,KNO3,M3n, BoxNumberXY_ben
#endif
  use constants,  ONLY: SEC_PER_DAY, ONE_PER_DAY, BENTHIC_BIO,STANDARD,EQUATION
  use mem_Param,  ONLY: p_poro, p_small, p_d_tot, CalcBenthicFlag
  use mem_BenOxygen
  use bennut_interface, ONLY: CalculateFromSet
  use mem_Param,  ONLY: p_d_tot, p_clD1D2m

!  
!
! !AUTHORS
!   P. Ruardij
!
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij & M.Vichi
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
  real(RLEN),dimension(NO_BOXES_XY)  :: r
  real(RLEN),dimension(NO_BOXES_XY)  :: diff
  real(RLEN),dimension(NO_BOXES_XY)  :: zmG2o
  real(RLEN),dimension(NO_BOXES_XY)  :: D1mNew
  real(RLEN),dimension(NO_BOXES_XY)  :: G2oNew
  real(RLEN),dimension(NO_BOXES_XY)  :: jG2O2o
  real(RLEN),dimension(NO_BOXES_XY)  :: unc_shiftD1m
  real(RLEN)                             :: dummy


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Emperical equation derived from Broecker and Peng (1973)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  diff = SEC_PER_DAY* 1.0e-9_RLEN* (10.0_RLEN)**((- 984.26_RLEN/( 273.0_RLEN+ &
    ETW_Ben(:))+ 3.672_RLEN))* p_poro* irrenh(:)* p_exsaf

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Recalculate total consumption from /m2 to /m3 pw:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  zmG2o  =  ( rrBTo(:)+ jG2K3o(:)+ jG2K7o(:))/( D1m(:)* p_poro)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Determine new thickness:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  D1mNew  =   sqrt(  2.0_RLEN* diff* O2o_Ben(:)/( p_small+ zmG2o))
  D1mNew  = min(D1mnew ,p_d_tot-2.0_RLEN*p_clD1D2m);

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate rate of change of thickness of the aerobic layer:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  shiftD1m(:)  =  ( max(  p_mD1m,  D1mNew)- D1m(:))/ ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Damping the change of D1m in case of large changes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( InitializeModel== 0) then
     shiftD1m(:) = shiftD1m(:)* (D1m(:)/( D1m(:)+ &
      abs(shiftD1m(:))))**(p_xdampingD1m)*( p_chD1m/( p_chD1m+ D1m(:)))
     if ( CalcBenthicFlag > BENTHIC_BIO) then
        do BoxNumberXY_ben = 1,NO_BOXES_XY
           r(1) = CalculateFromSet( KNO3(BoxNumberXY_ben), EQUATION, &
                  STANDARD, D1m(BoxNumberXY_ben), dummy)/M3n(BoxNumberXY_ben)
           if ( r(1) .lt.ZERO) then
             write(LOGUNIT,*) "BFM Warning: BenOxygen proportion M3n(D1m)/M3n(0..D2m)=",r(1)
           endif
        end do
     end if
     shiftD1m(:)= shiftD1m(:) * max(ZERO,min(ONE,r*2.0_RLEN));
  end if


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Damping the change of D1m in case of too thick D1m
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  r= min( shiftD1m(:),max(ZERO,p_d_tot-p_chD1m-D1m(:)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! recalculate the new D1mNew at the actual time step:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  D1mNew  =   D1m(:)+ r* LocalDelta

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! calculate the consumption which belongs to the corrected D1mNew
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  zmG2o  =  ( 2.0_RLEN* diff* O2o_Ben(:))/( D1mNew* D1mNew)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! New oxygen conc. in the sediment:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  G2oNew = D1mNew*( O2o_Ben(:)- 0.66667_RLEN* zmG2o* D1mNew* D1mNew/( 2.0_RLEN* &
    diff))* p_poro

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! flux to pelagic: correct flux for rate of change of G2o
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  jG2O2o = -( rrBTo(:)+ jG2K3o(:)+ jG2K7o(:))-( G2oNew- G2o(:))/ &
    ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Assign fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( InitializeModel== 0) then
    shiftD1m(:)=r
    call flux_vector( iiBen, ppD1m,ppD1m, r )
    call flux_vector( iiBen, ppG2o,ppG2o,-jG2O2o )
    jbotO2o(:)=jbotO2o(:)+jG2O2o
    if ( CalcBenthicFlag== BENTHIC_BIO) then
      ! Compute shifting of the denitrification layer here in case of running only
      ! the benthic submodel and NOT the benthic nutrient model.
      r  =   p_d_tot- D2m(:)
      call flux_vector( iiBen, ppD2m,ppD2m, shiftD1m(:)* r/( r+ 0.01_RLEN) )
    end if
  else
      G2o(:)=G2oNew
  endif

  end subroutine BenOxygenDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
