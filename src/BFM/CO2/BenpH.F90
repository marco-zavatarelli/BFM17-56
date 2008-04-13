#include "INCLUDE.h"
#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenpH
!
! DESCRIPTION
!   Computation of pH in sediments according to the carbonate system 
!   equations (see ModuleCO2System.F90)
!
! !INTERFACE
  subroutine BenpHDynamics
!

#ifdef INCLUDE_BENCO2

! !USES:

  ! For the following Benthic-states fluxes are defined: G13c, G3c
  ! The following Benthic-states are used (NOT in fluxes): D1m, Q1c, D6m, D2m
  ! The following global vars are modified: dummy
  ! The following global scalar vars are used: &
  !    NO_BOXES_XY,  &
  !   BoxNumberXY, InitializeModel, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : KCO2, jbotO3c
  ! The following Benthic 1-d global boxvars got a value: DICae, DICan
  ! The following Benthic 1-d global boxvars are used: rrBTo, KQ1, &
  ! irrenh, ETW_Ben, rrATo, O3c_Ben, shiftD1m
  ! The following Benthic 2-d global boxvars  are used: ruHI
  ! The following groupmember vars  are used: iiH1
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global parameters are used: p_d_tot, p_q10diff
  ! The following global constants are used: RLEN
  ! The following constants are used: GET, &
  ! LABDA_1, COEFFICIENT, LAYERS, LAYER1, LAYER2, &
  ! DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DOUBLE_DEFINE, &
  ! ZERO_EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, LINEAR_TERM, CONSTANT_TERM, &
  ! SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, EQUATION, &
  ! INPUT_TERM, PARAMETER, START_ADD_TERM, INPUT_ADD_TERM, &
  ! SET_LAYER_INTEGRAL_UNTIL, LAYER3, SET_LAYER_INTEGRAL, ADD, DERIVATIVE, &
  ! RFLUX, SHIFT, ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
  use mem,  ONLY: G13c, G3c, D1m, Q1c, D6m, D2m, D2STATE
  use mem, ONLY: ppG13c, ppG3c, ppD1m, ppQ1c, ppD6m, ppD2m, &
    dummy,    NO_BOXES_XY,   &
     BoxNumberXY, DICae,  pHAe, pCO2ae, DICan,  pHan, pCO2an,  ETW_Ben, &
    ESW_Ben, ERHO_Ben, M1p, M5s,AcAe, AcAn,M11p,M21p,D1m,D2m
  USE BFM_ERROR_MSG, ONLY: BFM_ERROR
  use CO2System,ONLY: CalcCO2System
  use mem_Param,  ONLY: p_d_tot 
  use mem_CO2, ONLY: DYNAMIC
  IMPLICIT NONE
!  
! !LOCAL VARIABLES
  real(RLEN)  :: r
  real(RLEN)  :: CO2
  real(RLEN)  :: HCO3
  real(RLEN)  :: CO3
  real(RLEN)  :: m1
  integer     :: error
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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

  do BoxNumberXY=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pH value in oxic sediments
      ! Only the iterative solution of the carbonate system can be
      ! used
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       error= CalcCO2System(DYNAMIC,ESW_Ben(BoxNumberXY),& 
                   ETW_Ben(BoxNumberXY),ERHO_Ben(BoxNumberXY),&
                   M1p(BoxNumberXY),M5s(BoxNumberXY),Acae(BoxNumberXY),&
                   CO2,HCO3,CO3,pHae(BoxNumberXY),&
                   DIC_in=DICae(BoxNumberXY),pCO2_out=pCO2ae(BoxNumberXY))
       if ( error > 0 ) then
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ESW_Ben',ESW_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ETW_Ben',ETW_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ERHO_Ben',ERHO_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICae',DICae(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M1p',M1p(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M5s',M5s(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acae',Acae(BoxNumberXY)
            write(LOGUNIT,'('' pHae='',G12.6)') pHae(BoxNumberXY)
            write(LOGUNIT,*) "BenpHDynamics pHae outside range 2-11"
            pHae(BoxNumberXY)=-1
!           call BFM_ERROR("BenpHDynamics","pHae outside range 2-11")
       endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pH value in anoxic sediments
      ! Only the iterative solution of the carbonate system can be
      ! used
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       m1=M11p(BoxNumberXY)*(D2m(BoxNumberXY)-D1m(BoxNumberXY)) + &
       M21p(BoxNumberXY)*(p_d_tot-D2m(BoxNumberXY))/ (p_d_tot-D1m(BoxNumberXY))
       error= CalcCO2System(DYNAMIC,ESW_Ben(BoxNumberXY),& 
                   ETW_Ben(BoxNumberXY),ERHO_Ben(BoxNumberXY),&
                   m1,M5s(BoxNumberXY),Acan(BoxNumberXY),&
                   CO2,HCO3,CO3,pHan(BoxNumberXY),&
                   DIC_in=DICan(BoxNumberXY),pCO2_out=pCO2an(BoxNumberXY))
       if ( error > 0 ) then
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICan',DICan(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acan',Acan(BoxNumberXY)
            write(LOGUNIT,'('' pHan='',G12.6)') pHan(BoxNumberXY)
            write(LOGUNIT,*) "BenpHDynamics:pHan outside range 2-11"
            pHan(BoxNumberXY)=-1
!           call BFM_ERROR("BenpHDynamics","pHan outside range 2-11")
       endif
  end do
#endif

  end subroutine BenpHDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
