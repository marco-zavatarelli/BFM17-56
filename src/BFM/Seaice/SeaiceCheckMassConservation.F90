#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CheckMassConservation
!
! DESCRIPTION
!   !
!
! !INTERFACE
  subroutine SeaiceCheckMassConservation
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
use global_mem, ONLY:RLEN,ZERO
  use constants, ONLY: MW_P, MW_N, MW_SI
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: T1c, T1p, T1n, I1p, I3n, I4n, I5s
  use mem, ONLY: ppT1c,ppT1p,ppT1n, &
     ppI1p, ppI3n, ppI4n, ppI5s, ppF3c, &
     qpcSZO,qncSZO
  use mem, ONLY: flux_vector,ppSeaiceZoo,SeaiceZoo, &
    totseaicec,totseaicep,totseaicen,totseaices, &
    iiSeaiceZoo,NO_BOXES_ICE,iiC,iiN,iiP,iiS,&
    SeaiceAlgae,iiSeaiceAlgae,ppSeaiceAlgae, &
    SeaiceDetritus,iiSeaiceDetritus,ppSeaiceDetritus
#endif  
!
! !AUTHORS
!   Letizia Tedesco and Marcello Vichi
!
!
!
! !REVISION_HISTORY
!   Created at Mon Nov 21 09:44:23 CET 2005
!
!
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2007 the BFM team
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
  real(RLEN),dimension(NO_BOXES_ICE) :: s
  integer                            :: i,j
  
  totseaicec(:)=ZERO
  totseaicep(:)=ZERO
  totseaicen(:)=ZERO
  totseaices(:)=ZERO
  do i=1, iiSeaiceAlgae
     s=SeaiceAlgae(i,iiC)
     totseaicec(:)=totseaicec(:) + s
     s=SeaiceAlgae(i,iiN)
     totseaicen(:)=totseaicen(:) + s
     s=SeaiceAlgae(i,iiP)
     totseaicep(:)=totseaicep(:) + s
     j=ppSeaiceAlgae(i,iiS)
     if ( j/=0) then
        s=SeaiceAlgae(i,iiS)
        totseaices(:)=totseaices(:) + s
     end if
  end do
  do i=1, iiSeaiceZoo
     j=max(1,ppSeaiceZoo(i,iiN))
     s=SeaiceZoo(i,j)
     if ( j==1) s=s*qncSZO(i, :)
     totseaicen(:)=totseaicen(:) + s
     j=max(1,ppSeaiceZoo(i,iiP))
     s=SeaiceZoo(i,j)
     if ( j==1) s=s*qpcSZO(i, :)
     totseaicep(:)=totseaicep(:) + s
  end do
  do i=1, iiSeaiceDetritus
     if ( ppSeaiceDetritus(i,iiC)/=0) then
        s=SeaiceDetritus(i,iiC)
        totseaicec(:)=totseaicec(:) + s
     end if
     if ( ppSeaiceDetritus(i,iiN)/=0) then
        s=SeaiceDetritus(i,iiN)
        totseaicen(:)=totseaicen(:) + s
     end if
     if ( ppSeaiceDetritus(i,iiP)/=0) then
        s=SeaiceDetritus(i,iiP)
        totseaicep(:)=totseaicep(:) + s
     end if
     if ( ppSeaiceDetritus(i,iiS)/=0) then
        s=SeaiceDetritus(i,iiS)
        totseaices(:)=totseaices(:) + s
     end if
  end do

  totseaicec(:) = totseaicec(:)+ T1c(:)

  ! Convert from default units to g 
  totseaicec(:) = totseaicec(:)/1000.0_RLEN
  totseaicen(:) = totseaicen(:)*MW_N/1000.0_RLEN
  totseaicep(:) = totseaicep(:)*MW_P/1000.0_RLEN
  totseaices(:) = totseaices(:)*MW_Si/1000.0_RLEN


  end subroutine SeaiceCheckMassConservation
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
