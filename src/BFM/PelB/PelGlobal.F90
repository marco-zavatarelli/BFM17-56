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
!   !
!
! !INTERFACE
  subroutine PelGlobalDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
  use mem_Param,  ONLY: p_small
  use mem_PelGlobal
  use mem
!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Tue Apr 20 09:11:59 AM CEST 2004
!
!
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
  integer  :: i

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Reset var in which silica fluxes are collected:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  flP1R6s(:)  =   0.0D+00
  flPTN6r(:)  =   0.0D+00
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in pelagic detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  qpR6c(:)  =   R6p(:)/( p_small+ R6c(:))
  qnR6c(:)  =   R6n(:)/( p_small+ R6c(:))
  qsR6c(:)  =   R6s(:)/( p_small+ R6c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in microzooplankton and HNAN
  ! in case of fixed quota qp_mz and qn_mz are one time calculated in the Initialize.F90
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiMicroZooPlankton)
    if ( ppMicroZooPlankton(i,iiP) > 0 ) &
      qp_mz(i,:)  =   MicroZooPlankton(i,iiP)/( p_small+ MicroZooPlankton(i,iiC))
    if ( ppMicroZooPlankton(i,iiN) > 0 ) &
      qn_mz(i,:)  =   MicroZooPlankton(i,iiN)/( p_small+ MicroZooPlankton(i,iiC))
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in omnivorous and herbivorous mesozooplankton
  ! in case of fixed quota qp_mz and qn_mz are one time calculated in the Initialize.F90
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiMesoZooPlankton)
    if ( ppMesoZooPlankton(i,iiP) > 0 ) &
      qpZc(i,:)  =   MesoZooPlankton(i,iiP)/( p_small+ MesoZooPlankton(i,iiC))
    if ( ppMesoZooPlankton(i,iiN) > 0 ) &
      qnZc(i,:)  =   MesoZooPlankton(i,iiN)/( p_small+ MesoZooPlankton(i,iiC))
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in phytoplankton
  ! Compute light prop.or chl. quota in phytoplankton (dep. on ChlLightFlag)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiPhytoPlankton)

    qpPc(i,:)  =   PhytoPlankton(i,iiP)/( p_small+ PhytoPlankton(i,iiC))
    qnPc(i,:)  =   PhytoPlankton(i,iiN)/( p_small+ PhytoPlankton(i,iiC))
    qlPc(i,:)  =   PhytoPlankton(i,iiL)/( p_small+ PhytoPlankton(i,iiC))
  end do

  qsPc(iiP1,:)  =   P1s(:)/( p_small+ P1c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in Pelagic Bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  qpB1c(:)  =   B1p(:)/( p_small+ B1c(:))
  qnB1c(:)  =   B1n(:)/( p_small+ B1c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute sedimentation velocities
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sediR2(:)     = ZERO
  sediR6(:)  =   p_rR6m
  do i = 1 , ( iiPhytoPlankton)
    sediPI(i,:)  =   p_rPIm( i)
  end do

  end subroutine PelGlobalDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
