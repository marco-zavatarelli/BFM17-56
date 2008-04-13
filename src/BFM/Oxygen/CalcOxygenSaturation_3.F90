#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcOxygenSaturation_3
!
! DESCRIPTION
!   !       Function for calculation of the Oxygen Saturation
!
!	based on:
!
!	WEISS 1970 DEEP SEA RES 17, 721-735.
!	units of ln(ml(STP)/l)
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE CalcOxygenSaturation()
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): O2o
  ! The following Pelagic 1-d global boxvars are modified : cxoO2
  ! The following Pelagic 1-d global boxvars got a value: eO2mO2
  ! The following Pelagic 1-d global boxvars  are used: ETW, ESW
  ! The following global constants are used: RLEN,ZERO_KELVIN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use constants, ONLY:ZERO_KELVIN
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: O2o
#endif
  use mem, ONLY: ppO2o, cxoO2, eO2mO2, ETW, ESW, NO_BOXES, iiBen, iiPel, &
    flux_vector
  use mem_Param,  ONLY: p_small
  IMPLICIT NONE
!  
!
! !AUTHORS
!   Piet Ruardij
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
  real(RLEN),dimension(NO_BOXES)  :: h
  real(RLEN),dimension(NO_BOXES)  :: abt
!EOP
!-------------------------------------------------------------------------!
!BOC

  ! calc absolute temperature divided by 100.0;
  ! input for the next empirical equation.
  abt  =  ( ETW(:)- ZERO_KELVIN)/ 100.0_RLEN

  !   calc theoretical oxygen saturation for temp + salinity
  !   From WEISS 1970 DEEP SEA RES 17, 721-735.
  !   units of ln(ml(STP)/l)
  h = - 173.4292_RLEN+ 249.6339_RLEN/ abt+ 143.3483_RLEN* dlog( &
    abt)- 21.8492_RLEN* abt+ ESW(:)*(- 0.033096_RLEN+ 0.014259_RLEN* abt- &
    0.0017_RLEN* (abt)**(2.0_RLEN))

  ! convert units to ml(STP)/l
  h  =   dexp(h)

  ! convert to mMol/m3
  !   calc volume of an ideal gas at standard temp (25C) and
  !   pressure (1.e-3 atm)
  !   p_videal = (8.3145 * 298.15 / 101325.0) = 24.4665e-3;
  cxoO2(:)  =   h/ 24.4665E-3_RLEN
  eO2mO2(:)  =   max(p_small,O2o(:))/ cxoO2(:)

  end subroutine CalcOxygenSaturation
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
