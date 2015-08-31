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
  use mem, ONLY: I1p, I3n, I4n, I5s
  use mem, ONLY: qpcSZO,qncSZO
  use mem, ONLY: flux_vector,ppSeaiceZoo,SeaiceZoo, &
    totseaicec,totseaicep,totseaicen,totseaices, &
    iiSeaiceZoo,iiC,iiN,iiP,iiS,&
    SeaiceAlgae,iiSeaiceAlgae,ppSeaiceAlgae, &
    SeaiceDetritus,iiSeaiceDetritus,ppSeaiceDetritus, &
    SeaiceBacteria,iiSeaiceBacteria
#endif  
!
! !AUTHORS
!   Letizia Tedesco and Marcello Vichi
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
  integer                            :: i,j
  
  totseaicec(:) = ZERO
  totseaicep(:) = I1p
  totseaicen(:) = I3n+I4n
  totseaices(:) = I5s
  do i=1, iiSeaiceAlgae
     totseaicec(:)=totseaicec(:) + SeaiceAlgae(i,iiC)
     totseaicen(:)=totseaicen(:) + SeaiceAlgae(i,iiN)
     totseaicep(:)=totseaicep(:) + SeaiceAlgae(i,iiP)
     j=ppSeaiceAlgae(i,iiS)
     if ( j/=0) then
        totseaices(:)=totseaices(:) + SeaiceAlgae(i,iiS)
     end if
  end do
  do i=1, iiSeaiceZoo
     j=max(1,ppSeaiceZoo(i,iiN))
     if ( j==1) &
        totseaicen(:)=totseaicen(:) + SeaiceZoo(i,j)*qncSZO(i, :)
     j=max(1,ppSeaiceZoo(i,iiP))
     if ( j==1) &
        totseaicep(:)=totseaicep(:) + SeaiceZoo(i,j)*qpcSZO(i, :)
  end do
  do i=1, iiSeaiceDetritus
     if ( ppSeaiceDetritus(i,iiC)/=0) then
        totseaicec(:)=totseaicec(:) + SeaiceDetritus(i,iiC)
     end if
     if ( ppSeaiceDetritus(i,iiN)/=0) then
        totseaicen(:)=totseaicen(:) + SeaiceDetritus(i,iiN)
     end if
     if ( ppSeaiceDetritus(i,iiP)/=0) then
        totseaicep(:)=totseaicep(:) + SeaiceDetritus(i,iiP)
     end if
     if ( ppSeaiceDetritus(i,iiS)/=0) then
        totseaices(:)=totseaices(:) + SeaiceDetritus(i,iiS)
     end if
  end do

  do i=1, iiSeaiceBacteria
     totseaicec(:)=totseaicec(:) + SeaiceBacteria(i,iiC)
     totseaicen(:)=totseaicen(:) + SeaiceBacteria(i,iiN)
     totseaicep(:)=totseaicep(:) + SeaiceBacteria(i,iiP)
  end do

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
