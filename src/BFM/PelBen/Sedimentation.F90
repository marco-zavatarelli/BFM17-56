#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Sedimentation
!
! DESCRIPTION
!   Define all fluxes of material which enters the benthic system:
!	a, fluxes of detritus (slow degradable and labile organic detritus)
!       b. Changes in organic matter distribution depths 
!
! !INTERFACE
  subroutine SedimentationDynamics
!
#ifdef INCLUDE_BEN
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: Q6c, Q6n, Q6p, Q6s, Q1c, Q1n, Q1p, D1m, D6m, D7m, D8m, D9m
  use mem, ONLY: NO_BOXES_XY, ppQ6c, ppQ6n, ppQ6p, ppQ6s, ppQ1c, ppQ1n, ppQ1p,&
    ppD6m, ppD7m, ppD8m, ppD9m, jbotR6c, jbotR6n, jbotR6p, jbotR6s, jbotR1c,  &
    jbotR1n, jbotR1p, iiBen, iiPel, flux_vector
#endif
!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Mon Nov 21 09:11:50 CET 2005
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, the BFM team
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
  real(RLEN)            :: newDm(NO_BOXES_XY)
  REAL(RLEN)            :: Delta
  real(RLEN), external  :: GetDelta

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculated Fluxes from Pelagic to Benthic
  ! These fluxes are the sum of sedimentation flux + the flux of
  ! material excreted by filterfeeders and originating from the pelagic.
  ! NOTE: ALL DETRITUS FLUXES FROM THE PELAGIC TO THE SEDIMENT ARE 
  !       DIRECTED TO Q6 VIA R6 IN SETTLING.F90
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector( iiBen, ppQ6c,ppQ6c, -jbotR6c(:) )
  call flux_vector( iiBen, ppQ6n,ppQ6n, -jbotR6n(:) )
  call flux_vector( iiBen, ppQ6p,ppQ6p, -jbotR6p(:) )
  call flux_vector( iiBen, ppQ6s,ppQ6s, -jbotR6s(:) )
  call flux_vector( iiBen, ppQ1c,ppQ1c, -jbotR1c(:) )
  call flux_vector( iiBen, ppQ1n,ppQ1n, -jbotR1n(:) )
  call flux_vector( iiBen, ppQ1p,ppQ1p, -jbotR1p(:) )

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes in the depth distribution state variables 
  ! due to sedimentation of detritus (burial is thus also included)
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   Delta=GetDelta( )
   call RecalcPenetrationDepth( D1m(:), D6m(:), &
        -jbotR6c(:)*Delta, Q6c(:),newDm(:) )
   call flux_vector(iiBen, ppD6m,ppD6m,(newDM(:)- D6m(:))/Delta)
   call RecalcPenetrationDepth( D1m(:), D7m(:), &
        -jbotR6n(:)*Delta, Q6n(:),newDm(:) )
   call flux_vector(iiBen, ppD7m,ppD7m,(newDM(:)- D7m(:))/Delta)
   call RecalcPenetrationDepth( D1m(:), D8m(:), &
        -jbotR6p(:)*Delta, Q6p(:),newDm(:) )
   call flux_vector(iiBen, ppD8m,ppD8m,(newDM(:)- D8m(:))/Delta)
   call RecalcPenetrationDepth(D1m(:), D9m(:), &
        -jbotR6s(:)*Delta, Q6s(:),newDm(:) )
   call flux_vector(iiBen, ppD9m,ppD9m,(newDM(:)- D9m(:))/Delta)

#endif
  end subroutine SedimentationDynamics

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
