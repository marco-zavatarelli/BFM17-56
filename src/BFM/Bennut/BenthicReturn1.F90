#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicReturn1
!
! DESCRIPTION
!   This process is a very simple parameterisation of benthic 
!   remineralisation.
!   Benthic organism sub-model cannot be used with this model.
!   A constant portion of the organic matter in the sediments is
!   released to the water column as inorganic nutrients.
!   Oxygen consumption is stoichiometrically associated to carbon 
!   remineralisation rates and nitrogen remineralisation is partitioned 
!   into ammonium and nitrate flux with a constant value. 
!   Fluxes are then used as boundary conditions for pelagic variables
!   in BentoPelCoup.F90
!
!
! !INTERFACE
  subroutine BenthicReturn1Dynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE_BEN, D2DIAGNOS
#else
  use mem,  ONLY: Q6c, Q1c, Q6p, Q1p, K1p, Q6n, Q1n, K3n, K4n, Q6s, K5s
#endif
  use mem, ONLY: ppQ6c, ppQ1c, ppQ6p, ppQ1p, ppK1p, ppQ6n, ppQ1n, &
    ppK3n, ppK4n, ppQ6s, ppK5s, jbotO2o, ppO2o, jbotN1p, ppN1p, jbotN3n, ppN3n, jbotN4n, ppN4n, jbotN5s, ppN5s, &
    NO_BOXES_XY, iiBen, iiPel, flux_vector
  use mem_BenthicReturn1
  use mem_Param, ONLY: CalcConservationFlag

!  
!
! !AUTHORS
!   Piet Ruardij 
!
!
!
! !REVISION_HISTORY
!   Created at Fri Apr 30 21:30:27 CEST 2004
!
!
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
  real(RLEN),dimension(NO_BOXES_XY)  :: rate
#ifndef IMPFLUX
  rate  =   p_reminQ6c* Q6c(:)
  call flux_vector( iiBen, ppQ6c,ppQ6c,-( rate) )
  jbotO2o(:)  =  - rate/ 12.D+00

  rate  =   p_reminQ1c* Q1c(:)
  call flux_vector( iiBen, ppQ1c,ppQ1c,-( rate) )
  jbotO2o(:)  =   jbotO2o(:)- rate/ 12.D+00

  rate  =   p_reminQ6p* Q6p(:)
  call flux_vector( iiBen, ppQ6p,ppQ6p,-( rate) )
  jbotN1p(:)  =   rate

  rate  =   p_reminQ1p* Q1p(:)
  call flux_vector( iiBen, ppQ1p,ppQ1p,-( rate) )
  jbotN1p(:)  =   jbotN1p(:)+ rate

  ! K1.p is not used in this model version
  ! activated only for mass conservation check
  if (CalcConservationFlag) &
  call flux_vector( iiBen, ppK1p,ppK1p,-(- jbotN1p(:)) )

  rate  =   p_reminQ6n* Q6n(:)
  call flux_vector( iiBen, ppQ6n,ppQ6n,-( rate) )
  jbotN3n(:)  =   rate* p_pQIN3
  jbotN4n(:)  =   rate*( 1.0D+00- p_pQIN3)

  rate  =   p_reminQ1n* Q1n(:)
  call flux_vector( iiBen, ppQ1n,ppQ1n,-( rate) )
  jbotN3n(:)  =   jbotN3n(:)+ rate* p_pQIN3
  jbotN4n(:)  =   jbotN4n(:)+ rate*( 1.0D+00- p_pQIN3)

  ! K3.n and K4.n are not used in this model version
  ! activated only for mass conservation check
  if (CalcConservationFlag) then
     call flux_vector( iiBen, ppK3n,ppK3n,-(- jbotN3n(:)) )
     call flux_vector( iiBen, ppK4n,ppK4n,-(- jbotN4n(:)) )
  end if

  rate  =   p_reminQ6s* Q6s(:)
  call flux_vector( iiBen, ppQ6s,ppQ6s,-( rate) )
  jbotN5s(:)  =   rate

  ! K5.s is not used in this model version
  ! activated only for mass conservation check
  if (CalcConservationFlag) &
  call flux_vector( iiBen, ppK5s,ppK5s,-(- jbotN5s(:)) )
#else
  rate = p_reminQ6c_imp
  call flux_vector( iiBen, ppQ6c,ppQ6c,-( rate) )
  jbotO2o(:)  =  - rate/ 12.D+00
  rate = p_reminQ6n_N3n_imp
  call flux_vector( iiBen, ppQ6n,ppQ6n,-( rate) )
  jbotN3n(:)  =   rate
  rate = p_reminQ6n_N4n_imp
  call flux_vector( iiBen, ppQ6n,ppQ6n,-( rate) )
  jbotN4n(:)  =   rate
  rate = p_reminQ6p_imp
  call flux_vector( iiBen, ppQ6p,ppQ6p,-( rate) )
  jbotN1p(:)  =   rate  
  rate = p_reminQ6s_imp
  call flux_vector( iiBen, ppQ6s,ppQ6s,-( rate) )
  jbotN5s(:)  =   rate
#endif
  end subroutine BenthicReturn1Dynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
