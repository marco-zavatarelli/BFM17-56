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
!
! !INTERFACE
  SUBROUTINE CalcOxygenSaturation()
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use constants, ONLY:ZERO_KELVIN
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: O2o
  use mem, ONLY: ppO2o, cxoO2, eO2mO2, ETW, ESW, NO_BOXES, iiBen, iiPel, &
    flux_vector
#endif
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
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij  (rua@nioz.nl)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
     integer,save ::first=0
     real(RLEN),allocatable,save,dimension(:) :: h,abt
     integer :: AllocStatus, DeallocStatus
!EOP
!-------------------------------------------------------------------------!
!BOC
  if (first==0) then
     first=1
     allocate(h(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating h"
     allocate(abt(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating abt"
  endif

  ! calc absolute temperature divided by 100.0;
  ! input for the next empirical equation.
  abt  =  ( ETW(:)- ZERO_KELVIN)/ 100.0_RLEN

  !   calc theoretical oxygen saturation for temp + salinity
  !   From WEISS 1970 DEEP SEA RES 17, 721-735.
  !   units of ln(ml(STP)/l)
  h = - 173.4292_RLEN+ 249.6339_RLEN/ abt+ 143.3483_RLEN* log(abt) - &
    21.8492_RLEN* abt+ ESW(:)*(- 0.033096_RLEN+ 0.014259_RLEN* abt- &
    0.0017_RLEN* (abt)**(2.0_RLEN))

  ! convert units to ml(STP)/l
  h  =   exp(h)

  ! convert ml/l to mMol/m3
  ! Use the volume of a mole of DO at STP : 22.391 l (ICES conversions)
  ! Conversion:  1 ml/l = 10^3/22.391 = 44.661 uMol/L [ = mMol/m3 ]
  cxoO2(:)  =   h * 44.661_RLEN
  eO2mO2(:)  =   max(p_small,O2o(:))/ cxoO2(:)

  end subroutine CalcOxygenSaturation
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
