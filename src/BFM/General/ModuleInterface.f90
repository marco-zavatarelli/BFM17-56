!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleInterface
!
! DESCRIPTION
!   Definition of Explicit Interfaces
!
! !INTERFACE
  MODULE global_interface

! !AUTHORS
!   mfstep/ERSEM team, especially  Momme Butenschoen:
!
! !REVISION_HISTORY
!   ---
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
  INTERFACE

  subroutine MesoZooDynamics(zoo)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: zoo
  end subroutine MesoZooDynamics
  end INTERFACE

  INTERFACE

  subroutine MicroZooDynamics(zoo)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: zoo
  end subroutine MicroZooDynamics
  end INTERFACE

  INTERFACE

  subroutine PhytoDynamics(phyto)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: phyto
  end subroutine PhytoDynamics
  end INTERFACE

  INTERFACE

  subroutine PhotoAvailableRadiation(phyto)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: phyto
  end subroutine PhotoAvailableRadiation
  end INTERFACE

  INTERFACE

  subroutine LightAdaptationDynamics(phyto)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: phyto
  end subroutine LightAdaptationDynamics
  end INTERFACE

  INTERFACE

  subroutine BenOrganismDynamics(y, ppyc, ppyn, ppyp)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: y
  integer,intent(IN) :: ppyc
  integer,intent(IN) :: ppyn
  integer,intent(IN) :: ppyp
  end subroutine BenOrganismDynamics
  end INTERFACE

  INTERFACE

  subroutine BenBacDynamics(hx, pphxc, pphxn, pphxp)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: hx
  integer,intent(IN) :: pphxc
  integer,intent(IN) :: pphxn
  integer,intent(IN) :: pphxp
  end subroutine BenBacDynamics
  end INTERFACE

  INTERFACE

  FUNCTION eTq(temp,p_q10)
  use global_mem, ONLY:RLEN
  implicit none
  real(RLEN),intent(IN) :: temp
  real(RLEN),intent(IN) :: p_q10
  real(RLEN) :: eTq
  end FUNCTION eTq
  end INTERFACE

  INTERFACE

  SUBROUTINE CalcChlorophylla()
  use global_mem, ONLY:RLEN
  implicit none
  end SUBROUTINE CalcChlorophylla
  end INTERFACE

  INTERFACE

  SUBROUTINE CalcVerticalExtinction()
  use global_mem, ONLY:RLEN
  implicit none
  end SUBROUTINE CalcVerticalExtinction
  end INTERFACE

  INTERFACE

  SUBROUTINE CalcOxygenSaturation()
  use global_mem, ONLY:RLEN
  implicit none
  end SUBROUTINE CalcOxygenSaturation
  end INTERFACE

  INTERFACE

  FUNCTION BoxAbove(box_no)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: box_no
  integer :: BoxAbove
  end FUNCTION BoxAbove
  end INTERFACE

  INTERFACE

  FUNCTION BoxBeneath(box_no)
  use global_mem, ONLY:RLEN
  implicit none
  integer,intent(IN) :: box_no
  integer :: BoxBeneath
  end FUNCTION BoxBeneath
  end INTERFACE

  INTERFACE

  FUNCTION CalcSchmidtNumberOx(Temp)
  use global_mem, ONLY:RLEN
  implicit none
  real(RLEN),intent(IN) :: Temp
  real(RLEN) :: CalcSchmidtNumberOx
  end FUNCTION CalcSchmidtNumberOx
  end INTERFACE

  INTERFACE

  FUNCTION CalcSchmidtNumberCO2(Temp)
  use global_mem, ONLY:RLEN
  implicit none
  real(RLEN),intent(IN) :: Temp
  real(RLEN) :: CalcSchmidtNumberCO2
  end FUNCTION CalcSchmidtNumberCO2
  end INTERFACE

  end MODULE global_interface
