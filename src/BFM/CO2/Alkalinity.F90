#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Alkalinity
!
! DESCRIPTION
!
! !INTERFACE
  subroutine AlkalinityDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifdef INCLUDE_PELCO2
use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: ppO3h, ppN3n, ppN6r, Source_D3_vector, NO_BOXES, iiBen, iiPel, &
    flux_vector, ppN4n
#endif
  use mem_param,  ONLY: p_qro
!  
!
! !AUTHORS
!   M. Vichi from Wolf-Gladrow et al. (2007)
!
!
!
! !REVISION_HISTORY
!   
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
  integer,save ::first=0 
  real(RLEN),allocatable,save,dimension(:) :: rateN
  integer :: AllocStatus
  if (first==0) then
     first=1
     allocate(rateN(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rateN in Alkalinity.F90"
  end if
 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Changes in alkalinity due to N uptake, (de)nitrification
  ! From Wolf-Gladrow et al. (2007)
  ! - the assimilation of 1 mole of nitrogen (atoms) leads to 
  ! (i) an increase of alkalinity by 1 mole when nitrate or 
  !     nitrite is the N source, 
  ! (ii) to a decrease of alkalinity by 1 mole when ammonia is used
  ! (iii) to no change of alkalinity when molecular N is the source
  ! - Nitrification leads to a decrease of TA by 2 moles per 
  !   mole of NO3-  formed 
  ! - Denitrification leads to an increase of TA by 1 mole per 
  !   mole of nitrate converted
  !
  ! It is computed this way
  ! net_uptakeNO3=dNO3/dt+denit-nit
  ! net_uptakeNH4=dNH4/dt+nit
  ! dTA/dt = -net_uptakeNO3 + net_uptakeNH4 - 2*nit + denit
  !        = -dNO3/dt+dNH4/dt 
  ! Sulfur reactions associated to reduction equivalents are not
  ! considered as included in the operational TA definition
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  rateN(:) = - Source_D3_vector(ppN3n) +  &
             Source_D3_vector(ppN4n) 
  call flux_vector( iiPel, ppO3h,ppO3h, rateN)

#endif

  end subroutine AlkalinityDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
