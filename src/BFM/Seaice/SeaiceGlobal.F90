#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
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
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: U6p, U6c, U6n, U6s, T1p, X1n, X1p, X1c, &
    U1p, U1n, U1c, &
    T1c, T1n, T1c, SeaiceZoo, SeaiceAlgae, SeaiceDetritus, ppSeaiceDetritus, &
    iiSeaiceZoo, iiSeaiceAlgae, iiSeaiceBacteria, iiSeaiceDetritus, &
    ppSeaiceZoo, ppSeaiceAlgae, ppSeaiceBacteria, SeaiceBacteria, &
    iiC, iiP, iiN, iiS, iiL, iiU6, iiU1
#endif
  use mem, ONLY: qnSc, qpSc, qsSc, qlSc, qpXc, qnXc, &
                 qpUc, qnUc, qsUc, qpTc, qnTc
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
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
  ! All sea ice algae may have a silica component (switched off in the 
  ! GlobalDefs file)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiSeaiceAlgae)
    qpSc(i,:)  =   SeaiceAlgae(i,iiP)/( p_small+ SeaiceAlgae(i,iiC))
    qnSc(i,:)  =   SeaiceAlgae(i,iiN)/( p_small+ SeaiceAlgae(i,iiC))
    qlSc(i,:)  =   SeaiceAlgae(i,iiL)/( p_small+ SeaiceAlgae(i,iiC))
    qsSc(i,:)  =   SeaiceAlgae(i,iiS)/( p_small+ SeaiceAlgae(i,iiC))
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute nutrient quota in seaice Bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiSeaiceBacteria)
     if ( ppSeaiceBacteria(i,iiP) > 0 ) &
        qpTc(i,:)  =   SeaiceBacteria(i,iiP)/( p_small+ SeaiceBacteria(i,iiC))
     if ( ppSeaiceZoo(i,iiN) > 0 ) &
        qnTc(i,:)  =   SeaiceBacteria(i,iiN)/( p_small+ SeaiceBacteria(i,iiC))
  end do

  end
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
