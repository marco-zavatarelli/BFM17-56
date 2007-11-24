#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: SeaiceGlobal
!
! DESCRIPTION
! Computation of the global diagnostics in the sea-ice (nutrient quota)
!
!
! !INTERFACE
  subroutine SeaiceGlobalDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: 23STATE
#else
  use mem, ONLY: U6p, U6c, U6n, U6s, S1s, S1c, T1p, &
    T1c, T1n, SeaiceZoo, SeaiceAlgae, SeaiceDetritus, &
    iiSeaiceZoo, iiSeaiceAlgae, iiSeaiceDetritus, ppSeaiceDetritus,&
    ppSeaiceZoo, ppSeaiceAlgae, iiC, iiP, iiN, iiS, iiL, iiU6, iiS1
#endif
  use mem, ONLY: qnSc, qpSc, qsSc, qlSc, qpXc, qnXc, &
                 qpUc, qnUc, qsUc, qpT1c, qnT1c
  use mem_Param,  ONLY: p_small
  use mem_PelGlobal

!  
!
! !AUTHORS
!   Marcello Vichi and Letizia Tedesco
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2007 the BFM team
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
  ! Compute nutrient quota in seaice detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiSeaiceDetritus)
    qpUc(i,:)  =   SeaiceDetritus(i,iiP)/( p_small+ SeaiceDetritus(i,iiC))
    qnUc(i,:)  =   SeaiceDetritus(i,iiN)/( p_small+ SeaiceDetritus(i,iiC))
  end do

  qsUc(iiU6,:)  =   U6s(:)/( p_small+ U6c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in seaicezoo
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiSeaiceZoo)
     if ( ppSeaiceZoo(i,iiP) > 0 ) &
        qpXc(i,:)  =   SeaiceZoo(i,iiP)/( p_small+ SeaiceZoo(i,iiC))
     if ( ppSeaiceZoo(i,iiN) > 0 ) &
        qnXc(i,:)  =   SeaiceZoo(i,iiN)/( p_small+ SeaiceZoo(i,iiC))
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in Seaicealgae
  ! Compute light prop.or chl. quota in Seaicealgae 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiSeaiceAlgae)
    qpSc(i,:)  =   SeaiceAlgae(i,iiP)/( p_small+ SeaiceAlgae(i,iiC))
    qnSc(i,:)  =   SeaiceAlgae(i,iiN)/( p_small+ SeaiceAlgae(i,iiC))
    qlSc(i,:)  =   SeaiceAlgae(i,iiL)/( p_small+ SeaiceAlgae(i,iiC))
  end do

  qsSc(iiS1,:)  =   S1s(:)/( p_small+ S1c(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in seaice Bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  qpT1c(:) = T1p(:)/( p_small+ T1c(:))
  qnT1c(:) = T1n(:)/( p_small+ T1c(:))

  end
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
