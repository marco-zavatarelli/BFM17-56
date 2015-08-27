#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PhotoAvailableRadiation
!
! DESCRIPTION
!   This process computes the depth-integrated PAR and expresses
!	this in eiPPY for each phytoplankton. 
!       The change in ELiPPY
!	(optimal irradiance Iopt) due to daily variations is calculated in 
!	LightAdaptation.p.
!	The daily irradiance in each compartment is calculated
!	in the Light/Light.f and Light/VerticalDistribution.f and passed 
!       in EIR. The extinction-coefficient xEPS (/m) is calculated in 
!	in CalcVerticalExtinction.f
!
!   # Switch between diffferent light-production functions ("light type"):
!     iswLtyp = 0 : Steele (old ERSEM)  y*exp(1-y)
!     iswLtyp = 1 : Steele (Simpson)    y*exp(1-y)
!     iswLtyp = 2 : Ebenhoeh            2y/(1+y^2)
!     iswLtyp = 3 : ramp                min(1,y)
!     iswLtyp = 4 : step                1 if y>1 , 0 elsewhere
!     iswLtyp = 5 : Smith_average
!     iswLtyp = 6 : Smith II (actual_Irr)		
!
!     with y = irradiance/optimal
! 
!
! !INTERFACE
  subroutine PhotoAvailableRadiation(phyto)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE,NOTRANSPORT
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D3STATE, PhytoPlankton
  use mem, ONLY: ppPhytoPlankton, D3STATETYPE, EIR, xEPS, Depth, ELiPPY, eiPPY, &
                 iiC, iiL, NO_BOXES, iiBen, iiPel, flux_vector
#endif
  use mem_PAR
  use mem_Phyto, ONLY: p_iswLtyp


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: insw_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: phyto

!  
!
! !AUTHORS
!   Original version by W. Ebenhoeh, Oldenburg University 
!                           Hanneke Baretta-Bekker, VKI
!       Translated to OpenSesame by Piet Ruardij, NIOZ
!	Dependency on phytoplankton species by M. Vichi, INGV
!
!
!
! !REVISION_HISTORY
!   File created on 8 feb. 1997
!	Modified by Daji Huang and JWB, 19/6/1998
!	Scaled the calculation of eiPPY to a max. of 1 (at SUNQ of 16h) JWB040930
!	(Re)introduced the Smith formulations and removed the eiPPY scaling JWB041006
!
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer,dimension(NO_BOXES)  :: rampcontrol
  real(RLEN),dimension(NO_BOXES)  :: pIRRZ
  real(RLEN),dimension(NO_BOXES)  :: pIRR0
  real(RLEN),dimension(NO_BOXES)  :: xd
  real(RLEN),dimension(NO_BOXES)  :: exfac
  real(RLEN),dimension(NO_BOXES)  :: f_0_noon
  real(RLEN),dimension(NO_BOXES)  :: f_z_noon
  real(RLEN),dimension(NO_BOXES)  :: f_0_afternoon
  real(RLEN),dimension(NO_BOXES)  :: f_z_afternoon
  real(RLEN),dimension(NO_BOXES)  :: f_0_mean
  real(RLEN),dimension(NO_BOXES)  :: f_z_mean
  real(RLEN),dimension(NO_BOXES)  :: pirr0_noon
  real(RLEN),dimension(NO_BOXES)  :: pirrz_noon
  real(RLEN),dimension(NO_BOXES)  :: pirr0_afternoon
  real(RLEN),dimension(NO_BOXES)  :: pirrz_afternoon
  real(RLEN),dimension(NO_BOXES)  :: lx0
  real(RLEN),dimension(NO_BOXES)  :: lxz
  real(RLEN),dimension(NO_BOXES)  :: corr_mean
  real(RLEN),dimension(NO_BOXES)  :: corr_irra
  real(RLEN),dimension(NO_BOXES)  :: corr_irrb
  real(RLEN),dimension(NO_BOXES)  :: noon_light
  real(RLEN),dimension(NO_BOXES)  :: afternoon_light
  real(RLEN),dimension(NO_BOXES)  :: sum_state_phyto


  ! Recalculate Optimal light from the transported Px.l
  i  =   ppPhytoPlankton(phyto,iiL)
  select case ( D3STATETYPE( i))

    case ( NOTRANSPORT )
      ELiPPY(phyto,:)  =   PhytoPlankton(phyto,iiL)

    case default
      ELiPPY(phyto,:)  =   PhytoPlankton(phyto,iiL)/ PhytoPlankton(phyto,iiC)

  end select

  noon_light  =   EIR(:)* 1.7596E+00_RLEN  ! magic number = 2 * 2PI/(4+PI)
  afternoon_light  =   EIR(:)* 1.0620E+00_RLEN  ! magic number = (sqrt(2)+1)/2* 2PI/(4+PI)

  xd  =   xEPS(:)* Depth(:)
  exfac  =   exp( - xd)


  pIRR0  =   EIR(:)/ ELiPPY(phyto,:)
  pirr0_noon  =   noon_light/ ELiPPY(phyto,:)
  pirr0_afternoon  =   afternoon_light/ ELiPPY(phyto,:)

  pIRRZ  =   pIRR0* exfac
  pirrz_noon  =   pirr0_noon* exfac
  pirrz_afternoon  =   pirr0_afternoon* exfac

  select case ( p_iswLtyp(phyto))

    case ( 0 )
      ! Steele - Di Toro
      f_0_mean  =   exp(  ONE- pIRR0)
      f_z_mean  =   exp(  ONE- pIRRZ)
      corr_mean  =  ( f_z_mean- f_0_mean)/ xd
      corr_irra  =  ( f_z_mean- f_0_mean)/ xd* 6.0E+00_RLEN
      corr_irrb  =   ZERO

    case ( 1 )
      ! Steele
      f_0_noon  =   exp(  ONE- pirr0_noon)
      f_z_noon  =   exp(  ONE- pirrz_noon)
      f_0_afternoon  =   exp(  ONE- pirr0_afternoon)
      f_z_afternoon  =   exp(  ONE- pirrz_afternoon)
      corr_irra  =  -( f_0_noon- f_z_noon)/ xd
      corr_irrb  =  -( f_0_afternoon- f_z_afternoon)/ xd

    case ( 2 )
      ! Ebenhoeh
      f_0_noon  =   atan(pirr0_noon)
      f_z_noon  =   atan(pirrz_noon)
      f_0_afternoon  =   atan(pirr0_afternoon)
      f_z_afternoon  =   atan(pirrz_afternoon)
      corr_irra  =   2.0E+00_RLEN*( f_0_noon- f_z_noon)/ xd
      corr_irrb  =   2.0E+00_RLEN*( f_0_afternoon- f_z_afternoon)/ xd

    case ( 3 )
      rampcontrol  =   2* int(insw_vector(  pirrz_noon- ONE))
      rampcontrol = min( 2, rampcontrol+ int(insw_vector( pirr0_noon- &
        ONE)))

        WHERE (( rampcontrol)==2)
          corr_irra  =  ( log(  pirr0_noon)- log(  pirrz_noon))/ xd

        ELSEWHERE (( rampcontrol)==1)
          corr_irra  =  ( ONE+ log(  pirr0_noon)- pirrz_noon)/ xd

        ELSEWHERE (( rampcontrol)==0)
          corr_irra  =  ( pirr0_noon- pirrz_noon)/ xd

      END WHERE

      rampcontrol  =   2* int(insw_vector(  pirrz_afternoon- ONE))
      rampcontrol = min( 2, rampcontrol+ int(insw_vector( pirr0_afternoon- ONE)))

      WHERE (( rampcontrol)==2)
          corr_irrb  =   ONE

        ELSEWHERE (( rampcontrol)==1)
          corr_irrb  =  ( ONE+ log(  pirr0_afternoon)- pirrz_afternoon)/ xd

        ELSEWHERE (( rampcontrol)==0)
          corr_irrb  =  ( pirr0_afternoon- pirrz_afternoon)/ xd

      END WHERE


    case ( 4 )
      !Step
      lx0  =   log(  pirr0_noon)
      lxz  =   log(  pirr0_afternoon)
      corr_irra  =  ( max(  ZERO,  lx0)- max(  ZERO,  lx0- xd))/ xd
      corr_irrb  =  ( max(  ZERO,  lxz)- max(  ZERO,  lxz- xd))/ xd

    case ( 5 )
      !Smith:
      f_0_mean = log( ONE+ sqrt( (ONE/( 1.D-80+ pIRR0* &
        exp( ONE)))**(2.0E+00_RLEN)+ ONE))
      f_z_mean = log( exfac+ sqrt( (ONE/( 1.D-80+ pIRR0* &
        exp( ONE)))**(2.0E+00_RLEN)+ exfac* exfac))
      corr_mean  =  ( f_0_mean- f_z_mean)/ xd
      corr_irra  =  ( f_0_mean- f_z_mean)/ xd* 6.0E+00_RLEN
      corr_irrb  =   ZERO

    case ( 6 )
      ! Smith II
      f_0_noon = log( ONE+ sqrt( (pirr0_noon* exp( ONE))**(- &
        2.0E+00_RLEN)+ ONE))
      f_z_noon = log( exfac+ sqrt( (pirr0_noon* exp( ONE))**(- &
        2.0E+00_RLEN)+ exfac* exfac))
      f_0_afternoon = log( ONE+ sqrt( (pirr0_afternoon* exp( &
        ONE))**(- 2.0E+00_RLEN)+ ONE))
      f_z_afternoon = log( exfac+ sqrt( (pirr0_afternoon* exp( ONE))**(- &
        2.0E+00_RLEN)+ exfac* exfac))
      corr_irra  =  ( f_0_noon- f_z_noon)/ xd
      corr_irrb  =  ( f_0_afternoon- f_z_afternoon)/ xd

  end select

  select case ( LightPeriodFlag)

    case ( 1 )
      !rectangular integration:
      eiPPY(phyto,:)  =   min(  ONE,  corr_mean)

    case ( 3 )
      !   Simpson integration is used as default:
      eiPPY(phyto,:)  =  ( corr_irra+ 4.0E+00_RLEN* corr_irrb)/ 6.0E+00_RLEN

  end select


  end subroutine PhotoAvailableRadiation
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
