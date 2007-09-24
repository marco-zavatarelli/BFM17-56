#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: InitBenthicNutrient3
!
! DESCRIPTION
!   Initialization of the diagenetic state variables  in the sediment 
!   MAV: This is very experimental and needs to be revised
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

USE BFM_ERROR_MSG, ONLY: set_warning_for_getm
use global_mem, ONLY:RLEN,LOGUNIT
use mem,  ONLY: LocalDelta, iiBen, iiPel, iiReset,flux 
use mem,  ONLY: K4n,K14n,K24n,K3n,K1p,K11p,K21p,K6r,K5s,D1m,D2m,D6m,D7m,G2o;
use mem,  ONLY: H1c, H2c,H1n,H2n,H1p,H2p,Q1c,Q1n,Q1p,Q11c,Q11n,Q11p,N3n_Ben,N4n_Ben
use mem,  ONLY: KNH4,KNO3,KRED,KPO4,KSIO3
use mem,  ONLY: ppD1m,ppG2o, ppD2m,reBTn, reBTp, reATn, reATp
use mem,  ONLY: jbotO2o,rrBTo,jG2K3o,jG2K7o,shiftD1m,shiftD2m
use mem,  ONLY: NO_BOXES_XY, BoxNumberXY

use bennut_interface, ONLY: CalculateFromSet

use constants, ONLY: INTEGRAL,MASS
use mem_BenSilica, ONLY: p_clD2m,  p_chD2m
use mem_Param, ONLY: p_d_tot,p_clD1D2m
use mem_BenDenitriDepth, ONLY:p_cmD2m

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! The following sesame functions are used:IntegralExp
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
use mem_globalfun,   ONLY: IntegralExp


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
real(RLEN)           :: cD2m, HT_0,alpha,HTc,HTn,HTp
real(RLEN),dimension(NO_BOXES_XY)  :: r
integer              :: i

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! user defined external functions
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
real(RLEN), external  :: GetDelta
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Get actual time step for the calculation of the transient profile of &
! the nutrients
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

LocalDelta  =   1.0
do BoxNumberXY = 1,NO_BOXES_XY

    ! Initial calculation only to get a reasonable value for D1m and D2m

    G2o(BoxNumberXY)=0.0D+00
    jG2K3o(BoxNumberXY)=0.0D+00
    jG2K7o(BoxNumberXY)=0.0D+00



    D1m(BoxNumberXY)=0.05;
    D2m(BoxNumberXY)=0.07;
    alpha=1.0/(1.D-80+D6m(BoxNumberXY));
    HTc = (H1c(BoxNumberXY)+H2c(BoxNumberXY)) 
    HTn = (H1n(BoxNumberXY)+H2n(BoxNumberXY)) 
    HTp = (H1p(BoxNumberXY)+H2p(BoxNumberXY)) 
    Q11c(BoxNumberXY)=0.1 * Q1c(BoxNumberXY)
    Q11n(BoxNumberXY)=0.1 * Q1n(BoxNumberXY)
    Q11p(BoxNumberXY)=0.1 * Q1p(BoxNumberXY)

    ! 4 iterations are usually sufficient to get reasonable values of D1,D2
    do i=1,4
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

    D1m(BoxNumberXY)=D1m(BoxNumberXY) + ShiftD1m(BoxNumberXY) * LocalDelta
    D1m(BoxNumberXY)=min(D1m(BoxNumberXY), p_d_tot-2.0 * p_clD1D2m)


    if ( i== 1 )  then
    r  =   p_d_tot- D1m(BoxNumberXY)
    D2m(BoxNumberXY)=D1m(BoxNumberXY)+p_clD1D2m
    !      write(LOGUNIT,'(''reBTn, reBTp, reATn, reATp:'',4G15.3)') &
    !            reBTn(BoxNumberXY), reBTp(BoxNumberXY), reATn(BoxNumberXY), reATp(BoxNumberXY)
    else
       call BenDenitriDepthDynamics

       D2m(BoxNumberXY)=D2m(BoxNumberXY) + ShiftD2m(BoxNumberXY)* LocalDelta
       D2m(BoxNumberXY)=min( max( D1m(BoxNumberXY)+p_clD1D2m, D2m(BoxNumberXY) ), &
                     p_d_tot-p_clD1D2m)
    endif

    G2o(BoxNumberXY)= G2o(BoxNumberXY)- (jbotO2o(BoxNumberXY) + &
        rrBTo(BoxNumberXY)+jG2K3o(BoxNumberXY)+jG2K7o(BoxNumberXY)) * LocalDelta

    call BenAmmoniumDynamics   

    call BenNitrateDynamics

    call BenAnoxicDynamics

    call flux(1,iiReset,1,1,0.0D+00)
    enddo


    K4n(BoxNumberXY) =CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, MASS, &
                                            0.0D+00, D1m(BoxNumberXY))
    K14n(BoxNumberXY)=CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, MASS, &
                                   D1m(BoxNumberXY),D2m(BoxNumberXY))
    K24n(BoxNumberXY)=CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, MASS, &
                                   D2m(BoxNumberXY),p_d_tot)

    K3n(BoxNumberXY) =CalculateFromSet( KNO3(BoxNumberXY), INTEGRAL, MASS, &
                                            0.0D+00, D2m(BoxNumberXY))

    K6r(BoxNumberXY)=CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, MASS, &
                                   0.0D+00,p_d_tot)
    !       if ( K14n(BoxNumberXY) > K24n(BoxNumberXY) ) then
    !         write(LOGUNIT,*) 'Unrealistic profile'
    !         write(LOGUNIT,'(''D1m='',F12.3,'' D2m='',F12.3)') D1m(BoxNumberXY), D2m(BoxNumberXY)
    !         write(LOGUNIT,'(''D6m='',F12.3,'' D7m='',F12.3)') D6m(BoxNumberXY), D7m(BoxNumberXY)
    !         write(LOGUNIT,'(''K3n='',F12.3,'' K4n='',F12.3)') K3n(BoxNumberXY), K4n(BoxNumberXY)
    !         write(LOGUNIT,'(''K14n='',F12.3,'' K24n='',F12.3)') K14n(BoxNumberXY), K24n(BoxNumberXY)
    !         write(LOGUNIT,'(''N3n='',F12.3,'' N4n='',F12.3)') N3n_Ben(BoxNumberXY), N4n_Ben(BoxNumberXY)
    !         call set_warning_for_getm
    !       endif

    call BenPhosphateDynamics

    call BenSilicaDynamics


    K1p(BoxNumberXY) =CalculateFromSet( KPO4(BoxNumberXY), INTEGRAL, MASS, &
                                          0.0D+00, D1m(BoxNumberXY))
    K11p(BoxNumberXY)=CalculateFromSet( KPO4(BoxNumberXY), INTEGRAL, MASS, &
                                   D1m(BoxNumberXY),D2m(BoxNumberXY))
    K21p(BoxNumberXY)=CalculateFromSet( KPO4(BoxNumberXY), INTEGRAL, MASS, &
                                   D2m(BoxNumberXY),p_d_tot)

    cD2m  =   min(  max(  D2m(BoxNumberXY),   p_clD2m),  p_chD2m)
    K5s(BoxNumberXY)=CalculateFromSet( KSIO3(BoxNumberXY), INTEGRAL, MASS, &
                                   0.0D+00,cD2m)

    !     call BenQ1TransportDynamics

enddo

end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
