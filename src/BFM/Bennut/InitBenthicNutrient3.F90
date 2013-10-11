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
  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: LocalDelta, iiBen, iiPel, iiReset,flux
  use mem,  ONLY: K4n,K14n,K24n,K3n,K1p,K11p,K21p,K6r,K16r,K26r,K5s,D1m,D2m,D6m,D7m,G2o;
  use mem,  ONLY: H1c, H2c,H1n,H2n,H1p,H2p,Q1c,Q1n,Q1p,Q11c,Q11n,Q11p,N3n_Ben,N4n_Ben
  use mem,  ONLY: KNH4,KNO3,KRED,KPO4,KSIO3
  use mem,  ONLY: ppD1m,ppG2o, ppD2m,reBTn, reBTp, reATn, reATp,ppG3c,rrATo,rrBTo
  use mem,  ONLY: jbotO2o,rrBTo,jG2K3o,jG2K7o,shiftD1m,shiftD2m,ETW_Ben
  use mem,  ONLY:    NO_BOXES_XY, BoxNumberXY_ben, ERHO_Ben
#ifdef INCLUDE_BENCO2
  use mem,  ONLY:KCO2,G3c,G13c,G23c,KALK,G3h,G13h,G23h,O3h_Ben
#endif
#endif
  use bennut_interface, ONLY: CalculateFromSet
  use constants, ONLY: INTEGRAL,MASS
  use mem_BenSilica, ONLY: p_clD2m,  p_chD2m,  p_chM5s, p_cvM5s,p_q10
  use mem_Param, ONLY: p_d_tot,p_d_tot_2, p_clD1D2m, p_small,p_poro,p_small
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
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
  LocalDelta  =   ONE/real(nn)
  BoxNumberXY_ben=1

  ! Initial calculation only to gete a reasonable value for D1m and D2m

  G2o(BoxNumberXY_ben)=ZERO
  jG2K3o(BoxNumberXY_ben)=ZERO
  jG2K7o(BoxNumberXY_ben)=ZERO
  
  D1m(BoxNumberXY_ben)=0.05_RLEN
  D2m(BoxNumberXY_ben)=0.07_RLEN
  alpha=ONE/(p_small+D6m(BoxNumberXY_ben));
  HTc = (H1c(BoxNumberXY_ben)+H2c(BoxNumberXY_ben)) 
  HTn = (H1n(BoxNumberXY_ben)+H2n(BoxNumberXY_ben)) 
  HTp = (H1p(BoxNumberXY_ben)+H2p(BoxNumberXY_ben)) 
  Q11c(BoxNumberXY_ben)=0.1_RLEN * Q1c(BoxNumberXY_ben)
  Q11n(BoxNumberXY_ben)=0.1_RLEN * Q1n(BoxNumberXY_ben)
  Q11p(BoxNumberXY_ben)=0.1_RLEN * Q1p(BoxNumberXY_ben)
  
  do i=1,p_N
    HT_0 = HTc / IntegralExp( - alpha,p_d_tot )
    H1c(BoxNumberXY_ben)=HT_0*IntegralExp(-alpha,D1m(BoxNumberXY_ben))
    H2c(BoxNumberXY_ben)= HTc-H1c(BoxNumberXY_ben)
    HT_0 = HTn / IntegralExp( - alpha,p_d_tot )
    H1n(BoxNumberXY_ben)=HT_0*IntegralExp(-alpha,D1m(BoxNumberXY_ben))
    H2n(BoxNumberXY_ben)= HTn-H1n(BoxNumberXY_ben)
    HT_0 = HTp / IntegralExp( - alpha,p_d_tot )
    H1p(BoxNumberXY_ben)=HT_0*IntegralExp(-alpha,D1m(BoxNumberXY_ben))
    H2p(BoxNumberXY_ben)= HTp-H1p(BoxNumberXY_ben)

    call BenthicSystemDynamics

    call BenOxygenDynamics

    D1m(BoxNumberXY_ben)=D1m(BoxNumberXY_ben) + ShiftD1m(BoxNumberXY_ben) * LocalDelta
    D1m(BoxNumberXY_ben)=min(D1m(BoxNumberXY_ben), p_d_tot-2.0_RLEN * p_clD1D2m)


    if ( i== 1 )  then
       r  =   p_d_tot- D1m(BoxNumberXY_ben)
       D2m(BoxNumberXY_ben)=D1m(BoxNumberXY_ben)+p_clD1D2m
    else
       call BenDenitriDepthDynamics

       D2m(BoxNumberXY_ben)=D2m(BoxNumberXY_ben) + ShiftD2m(BoxNumberXY_ben)* LocalDelta
       D2m(BoxNumberXY_ben)=min( &
            max( D1m(BoxNumberXY_ben)+p_clD1D2m, D2m(BoxNumberXY_ben) ), &
            p_d_tot-p_clD1D2m)
    endif

!   G2o(1)= G2o(1)- (jbotO2o(1) +rrBTo(1)+jG2K3o(1)+jG2K7o(1)) * LocalDelta

    call BenAmmoniumDynamics   

    call BenNitrateDynamics

    call BenAnoxicDynamics

    call flux(1,iiReset,1,1,0.0)

    if ( reAtn(1) .le. ZERO .or. reATp(BoxNumberXY_ben).le.0.0) then
       H1c(BoxNumberXY_ben)=0.5_RLEN*H1c(BoxNumberXY_ben)
       H2c(BoxNumberXY_ben)=0.5_RLEN*H2c(BoxNumberXY_ben)
       HTc=0.5_RLEN*HTc
       H1n(BoxNumberXY_ben)=0.5_RLEN*H1n(BoxNumberXY_ben)
       H2n(BoxNumberXY_ben)=0.5_RLEN*H2n(BoxNumberXY_ben)
       HTn=0.5_RLEN*HTn
       H1p(BoxNumberXY_ben)=0.5_RLEN*H1p(BoxNumberXY_ben)
       H2p(BoxNumberXY_ben)=0.5_RLEN*H2p(BoxNumberXY_ben)
       HTp=0.5_RLEN*HTp
    endif
  enddo

  do BoxNumberXY_ben=1,NO_BOXES_XY
        K4n(BoxNumberXY_ben) =max(p_small,CalculateFromSet( KNH4(BoxNumberXY_ben), INTEGRAL, MASS, &
                                                ZERO, D1m(BoxNumberXY_ben)))
        K14n(BoxNumberXY_ben)=max(p_small,CalculateFromSet( KNH4(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY_ben),D2m(BoxNumberXY_ben)))
        K24n(BoxNumberXY_ben)=max(p_small,CalculateFromSet( KNH4(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY_ben),p_d_tot_2))
        K3n(BoxNumberXY_ben) =max(p_small,CalculateFromSet( KNO3(BoxNumberXY_ben), INTEGRAL, MASS, &
                                                ZERO, D2m(BoxNumberXY_ben)))
        K6r(BoxNumberXY_ben)=max(p_small,CalculateFromSet( KRED(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       ZERO,D1m(BoxNumberXY_ben)))
        K16r(BoxNumberXY_ben)=max(0.1_RLEN,CalculateFromSet( KRED(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY_ben),D2m(BoxNumberXY_ben)))
        K26r(BoxNumberXY_ben)=max(0.1_RLEN,CalculateFromSet( KRED(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY_ben),p_d_tot_2))
  enddo

  call BenPhosphateDynamics

  call BenSilicaDynamics

  call BenQ1TransportDynamics


  do BoxNumberXY_ben=1,NO_BOXES_XY
      K1p(BoxNumberXY_ben) =CalculateFromSet( KPO4(BoxNumberXY_ben), INTEGRAL, MASS, &
                                              ZERO, D1m(BoxNumberXY_ben))
      K11p(BoxNumberXY_ben)=CalculateFromSet( KPO4(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY_ben),D2m(BoxNumberXY_ben))
      K21p(BoxNumberXY_ben)=CalculateFromSet( KPO4(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY_ben),p_d_tot_2)
  enddo

#ifdef INCLUDE_BENCO2
      call BenCO2TransportDynamics
      call BenAlkalinityDynamics
      do BoxNumberXY_ben=1,NO_BOXES_XY
            G3c(BoxNumberXY_ben) =CalculateFromSet( KCO2(BoxNumberXY_ben), INTEGRAL, MASS, &
                                              ZERO, D1m(BoxNumberXY_ben))
            G13c(BoxNumberXY_ben)=CalculateFromSet( KCO2(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D1m(BoxNumberXY_ben),D2m(BoxNumberXY_ben))
            G23c(BoxNumberXY_ben)=CalculateFromSet( KCO2(BoxNumberXY_ben), INTEGRAL, MASS, &
                                       D2m(BoxNumberXY_ben),p_d_tot)
            ! convert alkalinity from pelagic units (umol/kg) to sediment units (mmol/m2)
            G3h(BoxNumberXY_ben)=O3h_Ben(BoxNumberXY_ben)*D1m(BoxNumberXY_ben) &
                             *p_poro(BoxNumberXY_ben)*ERHO_Ben(BoxNumberXY_ben)/1000._RLEN
            G13h(BoxNumberXY_ben)= &
                 O3h_Ben(BoxNumberXY_ben)* &
                 (D2m(BoxNumberXY_ben)- &
                 D1m(BoxNumberXY_ben))* &
                 p_poro(BoxNumberXY_ben)* &
                 ERHO_Ben(BoxNumberXY_ben)/1000._RLEN
            G23h(BoxNumberXY_ben)= &
                 O3h_Ben(BoxNumberXY_ben)* &
                 (p_d_tot_2-D2m(BoxNumberXY_ben))* &
                 p_poro(BoxNumberXY_ben)*ERHO_Ben(BoxNumberXY_ben)/1000._RLEN
      enddo
#endif
  end subroutine InitBenthicNutrient3Dynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
