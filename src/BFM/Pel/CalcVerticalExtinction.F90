#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcVerticalExtinction
!
! DESCRIPTION
!   Calculates the vertical extinction.
!     
!
! !INTERFACE
  SUBROUTINE CalcVerticalExtinction()
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: R6c, PhytoPlankton
  use mem, ONLY: ppR6c, ppPhytoPlankton, xEPS, ABIO_eps, ESS, iiPhytoPlankton, &
                 iiC, iiL, NO_BOXES, iiBen, iiPel, flux_vector
#endif
  use mem_Param, ONLY: p_eps0, p_epsR6, p_epsESS, ChlDynamicsFlag
  use mem_Phyto, ONLY: p_qlcPPY, p_epsChla
!
! !AUTHORS
!   ERSEM-team
!
!
! !REVISION_HISTORY
!   !
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
  integer  :: i
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton



  select case ( p_eps0 == ZERO)

    case( .TRUE. )
      xEPS(:)  =   ABIO_eps(:) + p_epsR6* R6c(:)

    case( .FALSE. )
      xEPS(:)  =   p_eps0 + p_epsESS*ESS(:)+ p_epsR6* R6c(:)

  end select


  select case ( ChlDynamicsFlag)

    case ( 1 )
      do i = 1 , ( iiPhytoPlankton)
        lcl_PhytoPlankton => PhytoPlankton(i,iiC)
        xEPS(:) = xEPS(:) + p_epsChla(i) * p_qlcPPY(i) * lcl_PhytoPlankton
      end do

    case ( 2 )
      do i = 1 , ( iiPhytoPlankton)
        lcl_PhytoPlankton => PhytoPlankton(i,iiL)
        xEPS(:) = xEPS(:) + p_epsChla(i) * lcl_PhytoPlankton
      end do

  end select

  end subroutine CalcVerticalExtinction
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
