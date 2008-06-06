#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: InitBenthicNutrient3
!
! DESCRIPTION
!   Initialization of the diagenetic state variables  in the sediment 
!
!
!
! !INTERFACE
  subroutine InitBenthicNutrient3Dynamics
!
! !USES:
  ! The following global scalar vars got a value: LocalDelta
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use bfm_error_msg, ONLY: set_warning_for_getm
  use global_mem, ONLY:RLEN,LOGUNIT
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: LocalDelta, iiBen, iiPel, iiReset,flux
  use mem,  ONLY: K4n,K14n,K24n,K3n,K1p,K11p,K21p,K6r,K16r,K26r,K5s,D1m,D2m,D6m,D7m,G2o;
  use mem,  ONLY: H1c, H2c,H1n,H2n,H1p,H2p,Q1c,Q1n,Q1p,Q11c,Q11n,Q11p,N3n_Ben,N4n_Ben
  use mem,  ONLY: KNH4,KNO3,KRED,KPO4,KSIO3
  use mem,  ONLY: ppD1m,ppG2o, ppD2m,reBTn, reBTp, reATn, reATp,ppG3c,rrATo,rrBTo
  use mem,  ONLY: jbotO2o,rrBTo,jG2K3o,jG2K7o,shiftD1m,shiftD2m,ETW_Ben
  use mem,  ONLY:    NO_BOXES_XY, NO_BOXES_XY, &
                     BoxNumberXY, ERHO_Ben
#ifdef INCLUDE_BENCO2
  use mem,  ONLY:KCO2,G3c,G13c,G23c,KALK,G3h,G13h,G23h,O3h_Ben
#endif
#endif
  use bennut_interface, ONLY: CalculateFromSet
  use constants, ONLY: INTEGRAL,MASS
  use mem_BenSilica, ONLY: p_clD2m,  p_chD2m,  p_chM5s, p_cvM5s,p_q10
  use mem_Param, ONLY: p_d_tot,p_d_tot_2, p_clD1D2m, p_small,p_poro 
  use mem_BenDenitriDepth, ONLY:p_cmD2m

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq
  use mem_globalfun,   ONLY: IntegralExp
  use mem_BenthicNutrient3


!
! !AUTHORS
!   P. Ruardij
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij & Marcllo Vichi
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
   real(RLEN)           :: cD2m, HT_0,alpha,HTc,HTn,HTp,chM5s
   real(RLEN),dimension(NO_BOXES_XY)  :: r
   integer              :: i,nn

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), external  :: GetDelta
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Get actual time step for the calculation of the transient profile of &
  ! the nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  nn=max(10,p_N)
  LocalDelta  =   1.0/real(nn)
  BoxNumberXY=1.0;

  ! Initial caluclation only to gete a reaonable value for D1m and D2m

  G2o(1)=0.0D+00
  jG2K3o(1)=0.0D+00
  jG2K7o(1)=0.0D+00
  

  
  D1m(BoxNumberXY)=0.05;
  D2m(BoxNumberXY)=0.07;
  alpha=1.0/(1.D-80+D6m(BoxNumberXY));
  HTc = (H1c(BoxNumberXY)+H2c(BoxNumberXY)) 
  HTn = (H1n(BoxNumberXY)+H2n(BoxNumberXY)) 
  HTp = (H1p(BoxNumberXY)+H2p(BoxNumberXY)) 
  Q11c(BoxNumberXY)=0.1 * Q1c(BoxNumberXY)
  Q11n(BoxNumberXY)=0.1 * Q1n(BoxNumberXY)
  Q11p(BoxNumberXY)=0.1 * Q1p(BoxNumberXY)
  
  do i=1,p_N
    HT_0 = HTc / IntegralExp( - alpha,p_d_tot )
    H1c(BoxNumberXY)=HT_0*IntegralExp(-alpha,D1m(BoxNumberXY))
    H2c(BoxNumberXY)= HTc-H1c(BoxNUmberXY)
    HT_0 = HTn / IntegralExp( - alpha,p_d_tot )
    H1n(BoxNumberXY)=HT_0*IntegralExp(-alpha,D1m(BoxNumberXY))
    H2n(BoxNumberXY)= HTn-H1n(BoxNUmberXY)
    HT_0 = HTp / IntegralExp( - alpha,p_d_tot )
    H1p(BoxNumberXY)=HT_0*IntegralExp(-alpha,D1m(BoxNumberXY))
    H2p(BoxNumberXY)= HTp-H1p(BoxNUmberXY)

    call BenthicSystemDynamics


    call BenOxygenDynamics

    D1m(1)=D1m(1) + ShiftD1m(1) * LocalDelta
    D1m(1)=min(D1m(1), p_d_tot-2.0 * p_clD1D2m)


    if ( i== 1 )  then
       r  =   p_d_tot- D1m(1)
       D2m(1)=D1m(1)+p_clD1D2m
    else
       call BenDenitriDepthDynamics

       D2m(1)=D2m(1) + ShiftD2m(1)* LocalDelta
       D2m(1)=min( max( D1m(1)+p_clD1D2m, D2m(1) ), &
                     p_d_tot-p_clD1D2m)
    endif

!   G2o(1)= G2o(1)- (jbotO2o(1) +rrBTo(1)+jG2K3o(1)+jG2K7o(1)) * LocalDelta

    call BenAmmoniumDynamics   

    call BenNitrateDynamics

    call BenAnoxicDynamics

    call flux(1,iiReset,1,1,0.0D+00)

    if ( reAtn(1) .le. 0.0 .or. reATp(1).le.0.0) then
       H1c(BoxNumberXY)=0.5*H1c(BoxNumberXY)
       H2c(BoxNumberXY)=0.5*H2c(BoxNumberXY)
       HTc=0.5*HTc
       H1n(BoxNumberXY)=0.5*H1n(BoxNumberXY)
       H2n(BoxNumberXY)=0.5*H2n(BoxNumberXY)
       HTn=0.5*HTn
       H1p(BoxNumberXY)=0.5*H1p(BoxNumberXY)
       H2p(BoxNumberXY)=0.5*H2p(BoxNumberXY)
       HTp=0.5*HTp
    endif
  enddo

  do BoxNumberXY=1,NO_BOXES_XY
        K4n(BoxNumberXY) =max(p_small,CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, MASS, &
                                                0.0D+00, D1m(BoxNumberXY)))
        K14n(BoxNumberXY)=max(p_small,CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY),D2m(BoxNumberXY)))
        K24n(BoxNumberXY)=max(p_small,CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY),p_d_tot_2))
        K3n(BoxNumberXY) =max(p_small,CalculateFromSet( KNO3(BoxNumberXY), INTEGRAL, MASS, &
                                                0.0D+00, D2m(BoxNumberXY)))
        K6r(BoxNumberXY)=max(p_small,CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, MASS, &
                                       0.0D+00,D1m(BoxNumberXY)))
        K16r(BoxNumberXY)=max(0.1D+00,CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY),D2m(BoxNUmberXY)))
        K26r(BoxNumberXY)=max(0.1D+00,CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY),p_d_tot_2))
  enddo

  call BenPhosphateDynamics

  call BenSilicaDynamics

  call BenQ1TransportDynamics


  do BoxNumberXY=1,NO_BOXES_XY
      K1p(BoxNumberXY) =CalculateFromSet( KPO4(BoxNumberXY), INTEGRAL, MASS, &
                                              0.0D+00, D1m(BoxNumberXY))
      K11p(BoxNumberXY)=CalculateFromSet( KPO4(BoxNumberXY), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY),D2m(BoxNumberXY))
      K21p(BoxNumberXY)=CalculateFromSet( KPO4(BoxNumberXY), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY),p_d_tot_2)
  enddo

#ifdef INCLUDE_BENCO2
      call BenCO2TransportDynamics
      call BenAlkalinityDynamics
      do BoxNumberXY=1,NO_BOXES_XY
            G3c(BoxNumberXY) =CalculateFromSet( KCO2(BoxNumberXY), INTEGRAL, MASS, &
                                              0.0D+00, D1m(BoxNumberXY))
            G13c(BoxNumberXY)=CalculateFromSet( KCO2(BoxNumberXY), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY),D2m(boxNUmberXY))
            G23c(BoxNumberXY)=CalculateFromSet( KCO2(BoxNumberXY), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY),p_d_tot)
            ! convert alkalinity from pelagic units (umol/kg) to sediment units (mmol/m2)
            G3h(BoxNumberXY)=O3h_Ben(BoxNumberXY)*D1m(BoxNumberXY) &
                             *p_poro(BoxNumberXY)*ERHO_Ben(BoxNumberXY)/1000._RLEN
            G13h(BoxNumberXY)=O3h_Ben(BoxNumberXY)*(D2m(BoxNumberXY)-D1m(BoxNumberXY)) &
                              *p_poro(BoxNumberXY)*ERHO_Ben(BoxNumberXY)/1000._RLEN
            G23h(BoxNumberXY)=O3h_Ben(BoxNumberXY)*(p_d_tot_2-D2m(BoxNumberXY)) &
                              *p_poro(BoxNumberXY)*ERHO_Ben(BoxNumberXY)/1000._RLEN
      enddo
#endif
  end subroutine InitBenthicNutrient3Dynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
