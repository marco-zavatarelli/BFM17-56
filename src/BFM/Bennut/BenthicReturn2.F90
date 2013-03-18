#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicReturn2
!
! DESCRIPTION
!   This process is an intermediate parameterisation of benthic 
!   remineralisation processes.
!   The sub-model of benthic organisms is required to use this 
!   parameterisation (at least benthic bacteria).
!   Oxygen demand at the water-sed interface is associated to amount 
!   of reduction equivalents present in the sediments.
!   Nutrients are released to the water column at constant specific
!   rates, according to the pore-water concentration.
!   Nitrogen remineralisation is partitioned 
!   into ammonium and nitrate flux with a constant value. 
!   A constant portion of biogenic silica in the sediments is
!   released to the water column as silicate.
!
! !INTERFACE
  subroutine BenthicReturn2Dynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: G2o, K6r, K1p, K11p, K4n, K14n, Q6s, D1m, D2m
  use mem, ONLY: ppG2o, ppK6r, ppK1p, ppK11p, ppK4n, ppK14n, ppQ6s, &
    ppD1m, ppD2m, jbotN3n, jbotN4n, jG2K7o, jbotN1p, jbotN5s, NO_BOXES_XY, iiBen, &
    iiPel, flux_vector
#endif
  use mem_Param,  ONLY: p_qro, p_d_tot
  use mem_BenthicReturn2
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
  real(RLEN),dimension(NO_BOXES_XY)  :: rate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  ! Oxygen consumption in the sediments
  rate = p_reminO2* max( ZERO, - G2o(:)/ D1m(:)+ K6r(:)/ p_qro/( &
    p_d_tot- D1m(:)))* D1m(:)
  call flux_vector( iiBen, ppG2o,ppG2o,-( rate) )
  call flux_vector( iiBen, ppK6r,ppK6r,-( rate* p_qro) )
  jG2K7o(:)  =   rate

  !----------------------------------------------------------------------
  ! Phosphorus remineralization in the sediments
  !----------------------------------------------------------------------
  rate  =   p_reminN1* K1p(:)
  ! jbotN1p is used in BenPelCoup to define the pelagic flux
  call flux_vector( iiBen, ppK1p,ppK1p,-( rate) )
  jbotN1p(:)  =   rate
  rate  =   p_K11K1p* K11p(:)
  call flux_vector( iiBen, ppK11p,ppK1p, rate )
  ! K21.p is not used in this model version

  !----------------------------------------------------------------------
  ! Nitrogen remineralization in the sediments
  !----------------------------------------------------------------------
  rate  =   p_reminN4* K4n(:)
  ! K3.n is not used in this model version
  jbotN3n(:)  =   rate* p_pQIN3
  jbotN4n(:)  =   rate*( ONE - p_pQIN3)
  ! jbotN3n is used in BenPelCoup to define the pelagic flux
  ! jbotN4n is used in BenPelCoup to define the pelagic flux
  call flux_vector( iiBen, ppK4n,ppK4n,-( jbotN3n(:)+ jbotN4n(:)) )
  rate  =   p_K14K4n* K14n(:)
  call flux_vector( iiBen, ppK14n,ppK4n, rate )

  !----------------------------------------------------------------------
  ! Silicate remineralization in the sediments
  !----------------------------------------------------------------------
  rate  =   p_reminQ6s* Q6s(:)
  call flux_vector( iiBen, ppQ6s,ppQ6s,-( rate) )
  ! K5.s is not used in this model version
  jbotN5s(:)  =   rate

  end subroutine BenthicReturn2Dynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
