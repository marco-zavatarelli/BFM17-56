#include "DEBUG.h"
#include "INCLUDE.h"
#ifdef INCLUDE_BEN
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicSystem
!
! DESCRIPTION
!   !   This is the top level of the benthic submodel.
!       All the biological processes affecting the benthic dynamics are called
!       in a specified sequence according to the calculation flags.
!
!
! !INTERFACE
  subroutine BenthicSystemDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: ppY1c, ppY1n, ppY1p, ppY2c, ppY2n, ppY2p, ppY4c, ppY4n, &
    ppY4p, ppY5c, ppY5n, ppY5p, ppH1c, ppH1n, ppH1p, ppH2c, ppH2n, ppH2p, iiY1, &
    iiY2, iiY3, iiY4, iiY5, iiH1, iiH2, iiBen, iiPel, flux
#endif
  use mem_Param, ONLY: CalcBenOrganisms, CalcBenBacteria

!
!
! !AUTHORS
!   P. Ruardij and M. Vichi
!       
! !REVISION_HISTORY
!   !
!
! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi
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


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of turenh and irrenh in Bioturbation...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call BioturbationDynamics

  call BenGlobalDynamics

  if ( CalcBenOrganisms(iiY1)) then
    call BenOrganismDynamics( iiY1, ppY1c, ppY1n, ppY1p)
  end if

  if ( CalcBenOrganisms(iiY2)) then
    call BenOrganismDynamics( iiY2, ppY2c, ppY2n, ppY2p)
  end if

  if ( CalcBenOrganisms(iiY4)) then
    call BenOrganismDynamics( iiY4, ppY4c, ppY4n, ppY4p)
  end if

  if ( CalcBenOrganisms(iiY5)) then
    call BenOrganismDynamics( iiY5, ppY5c, ppY5n, ppY5p)
  end if

  if ( CalcBenOrganisms(iiY3)) then
    call FilterFeederDynamics
  end if

  if ( CalcBenBacteria(iiH1)) then
    call BenBacDynamics( iiH1, ppH1c, ppH1n, ppH1p)
  end if

  if ( CalcBenBacteria(iiH2)) then
    call BenBacDynamics( iiH2, ppH2c, ppH2n, ppH2p)
  end if

  end subroutine BenthicSystemDynamics
#endif
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
