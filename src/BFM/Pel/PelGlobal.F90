#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelGlobal
!
! DESCRIPTION
!   Compute global pelagic diagnostic variables that are needed in other
!   subroutines
!
! !INTERFACE
  subroutine PelGlobalDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ZERO
  use mem_Param,  ONLY: p_small
  use mem_PelGlobal
  use mem
  use mem_Settling
  use mem_Phyto,  ONLY: p_rPIm
#ifdef BFM_GOTM
  use bio_var, ONLY: BOTindices
#else
  use api_bfm, ONLY: BOTindices
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
!
! !AUTHORS
!   Piet Ruardij
!
! !REVISION_HISTORY
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
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Reset total flux variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  flPTN6r(:)  =   ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in pelagic organic matter
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1, iiPelDetritus
    if ( ppPelDetritus(i,iiP) > 0 ) &
      qpcOMT(i,:)  =  PelDetritus(i,iiP)/( p_small+ PelDetritus(i,iiC))
    if ( ppPelDetritus(i,iiN) > 0 ) &
      qncOMT(i,:)  =  PelDetritus(i,iiN)/( p_small+ PelDetritus(i,iiC))
    if ( ppPelDetritus(i,iiS) > 0 ) &
      qscOMT(i,:)  =   PelDetritus(i,iiS)/( p_small+ PelDetritus(i,iiC))
#ifdef INCLUDE_PELFE
    if ( ppPelDetritus(i,iiF) > 0 ) &
      qfcOMT(i,:)  =   PelDetritus(i,iiF)/( p_small+ PelDetritus(i,iiC))
#endif
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in microzooplankton
  ! in case of fixed quota the values are constant and assigned 
  ! in Initialize.F90
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1, iiMicroZooPlankton
    if ( ppMicroZooPlankton(i,iiP) > 0 ) &
      qpcMIZ(i,:)  =   MicroZooPlankton(i,iiP)/( p_small+ MicroZooPlankton(i,iiC))
    if ( ppMicroZooPlankton(i,iiN) > 0 ) &
      qncMIZ(i,:)  =   MicroZooPlankton(i,iiN)/( p_small+ MicroZooPlankton(i,iiC))
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in mesozooplankton
  ! in case of fixed quota the values are constant and assigned 
  ! in Initialize.F90
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1, iiMesoZooPlankton
    if ( ppMesoZooPlankton(i,iiP) > 0 ) &
      qpcMEZ(i,:)  =   MesoZooPlankton(i,iiP)/( p_small+ MesoZooPlankton(i,iiC))
    if ( ppMesoZooPlankton(i,iiN) > 0 ) &
      qncMEZ(i,:)  =   MesoZooPlankton(i,iiN)/( p_small+ MesoZooPlankton(i,iiC))
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute constituents quota in phytoplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1, iiPhytoPlankton
    if ( ppPhytoPlankton(i,iiP) > 0 ) &
      qpcPPY(i,:)  =   PhytoPlankton(i,iiP)/( p_small+ PhytoPlankton(i,iiC))
    if ( ppPhytoPlankton(i,iiN) > 0 ) &
      qncPPY(i,:)  =   PhytoPlankton(i,iiN)/( p_small+ PhytoPlankton(i,iiC))
    if ( ppPhytoPlankton(i,iiS) > 0 ) &
      qscPPY(i,:)  =   PhytoPlankton(i,iiS)/( p_small+ PhytoPlankton(i,iiC))
    if ( ppPhytoPlankton(i,iiL) > 0 ) &
      qlcPPY(i,:)  =   PhytoPlankton(i,iiL)/( p_small+ PhytoPlankton(i,iiC))
#ifdef INCLUDE_PELFE
    if ( ppPhytoPlankton(i,iiF) > 0 ) &
      qfcPPY(i,:)  =   PhytoPlankton(i,iiF)/( p_small+ PhytoPlankton(i,iiC))
#endif
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in Pelagic Bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1, (iiPelBacteria)
    qpcPBA(i,:)  =   PelBacteria(i,iiP)/( p_small+ PelBacteria(i,iiC))
    qncPBA(i,:)  =   PelBacteria(i,iiN)/( p_small+ PelBacteria(i,iiC))
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Prescribe background costant sedimentation velocities
  ! The velocity at the bottom interface is set equal to the burial
  ! velocity. This needs to be removed if there is no bottom level.
  ! Check also Settling.F90
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sediR2(:)     = ZERO
  sediR2(BOTindices) = ZERO
  sediR6(:)  =   p_rR6m
  sediR6(BOTindices) = p_burvel_R6
  do i = 1 , ( iiPhytoPlankton)
    sediPPY(i,:)  =   p_rPIm( i)
    sediPPY(i,BOTindices)  =   p_burvel_PI
  end do

  end subroutine PelGlobalDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
