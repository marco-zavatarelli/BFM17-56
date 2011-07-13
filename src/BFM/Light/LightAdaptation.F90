#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: LightAdaptation
!
! DESCRIPTION
!   This routine describes the photoadaptation of the phytoplankton to 
!	the prevailing irradiance level at depth
!
! !INTERFACE
  subroutine LightAdaptationDynamics(phyto, ppphytoc, ppphyton, ppphytop, &
    ppphytos, ppphytol)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,NOTRANSPORT
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: D3STATE, PhytoPlankton
  use mem, ONLY: ppPhytoPlankton, D3STATETYPE, Depth, xEPS, EIR, EPLi, &
    iiL, iiC, Source_D3_vector, NO_BOXES, iiBen, iiPel, flux_vector
#endif
  use mem_LightAdaptation


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol

!  
!
! !AUTHORS
!   Original version by W. Ebenhoeh, Oldenburg University
!                           Hanneke Baretta-Bekker, VKI
!       Translated to OpenSesame by Piet Ruardij
!	Phytoplankton species dependency added by M. Vichi, INGV
!
!
!
! !REVISION_HISTORY
!   File created on 8 feb. 1997
!	Modified by Daji and JWB, 19/6/1998
!	Checked by D.Mills and JWB 030429
!	Horrible error removed in addepth
!
!
! COPYING
!   
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES) :: phytoc
  real(RLEN),dimension(NO_BOXES) :: phyton
  real(RLEN),dimension(NO_BOXES) :: phytop
  real(RLEN),dimension(NO_BOXES) :: phytos
  real(RLEN),dimension(NO_BOXES) :: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: iphytol
  real(RLEN),dimension(NO_BOXES)  :: addepth
  real(RLEN),dimension(NO_BOXES)  :: adfactor
  real(RLEN),dimension(NO_BOXES)  :: rate_PLi
  real(RLEN),dimension(NO_BOXES)  :: rate_EPLi
  real(RLEN),dimension(NO_BOXES)  :: new_EPLi
  real(RLEN),dimension(NO_BOXES)  :: eir_c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#if ! defined GFORTRAN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  phytoc = D3STATE(ppphytoc,:)
  phyton = D3STATE(ppphyton,:)
  phytop = D3STATE(ppphytop,:)
  phytos = D3STATE(ppphytos,:)
  phytol = D3STATE(ppphytol,:)


  ! EPLi[%phyto] has been already calculated 

  ! computation of the change in IRR_OPT
  ! adaptation to the light in the depth p_addepth below surface.
  ! p_addepth is superseded if this function is used in a
  ! high resolution vertical grid. In that case the central depth
  ! of the layer is used.

  addepth  =   min(  Depth(:),  p_addepth(phyto))
  adfactor  =   exp( - xEPS(:)* addepth)

  eir_c  =   EIR(:)* adfactor

  select case ( p_isw(phyto))

    case ( 1 )
      new_EPLi  =   max(  eir_c,  p_clEPLi(phyto))
      new_EPLi  =   min(  new_EPLi,  p_chEPLi(phyto))

    case ( 2 )
      new_EPLi = max( 2.0E+00_RLEN* eir_c* p_chEPLi(phyto)/( &
        eir_c+ p_chEPLi(phyto)), p_clEPLi(phyto))
      new_EPLi  =   min(  new_EPLi,  p_chEPLi(phyto))

  end select


  ! Speed of adaptation is controlled by p_ruPLi ( 1 maximum speed
  !                        0 no adaptation )

  iphytol  =   ppPhytoPlankton(phyto,iiL)
  select case ( D3STATETYPE( iphytol))

    case ( NOTRANSPORT )

      rate_EPLi  =   p_ruEPLi(phyto)*( new_EPLi(:)- EPLi(phyto,:))
      call flux_vector( iiPel, iphytol,iphytol, rate_EPLi )

    case default

      rate_PLi = Source_D3_vector(ppPhytoPlankton(phyto,iiC))* EPLi(phyto,:)+ &
                 p_ruEPLi(phyto)*( new_EPLi- EPLi(phyto,:))* phytoc
      call flux_vector( iiPel, iphytol,iphytol, rate_PLi )

  end select
#endif

  end subroutine LightAdaptationDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
