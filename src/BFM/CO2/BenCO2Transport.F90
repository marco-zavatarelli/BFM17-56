#include "INCLUDE.h"
#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenCO2Transport
!
! DESCRIPTION
!   Description of the anoxic diagenetic processes in the sediment
!
! !INTERFACE
  subroutine BenCO2TransportDynamics
!

#ifdef INCLUDE_BENCO2
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: G23c,G13c, G3c, D1m, Q1c, D6m, D2m, D2STATE_BEN
  use mem, ONLY: ppG23c, ppG13c, ppG3c, ppD1m, ppQ1c, ppD6m, ppD2m, &
    NO_BOXES_XY,  ERHO_ben, &
    BoxNumberXY_ben, InitializeModel, LocalDelta, KCO2, jbotO3c, DICae, &
    DICan, rrBTo, irrenh, ETW_Ben, rrATo, O3c_Ben, shiftD1m, shiftD2m, ruHI, iiH1, &
    Depth_Ben,iiBen, iiPel, flux
#endif
  use constants, ONLY: MW_C, GET, LABDA_1, COEFFICIENT, &
    LAYERS, LAYER1, LAYER2, DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, DOUBLE_DEFINE, ZERO_EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, &
    LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, &
    EQUATION, INPUT_TERM, PARAMETER, START_ADD_TERM, &
    INPUT_ADD_TERM, SET_LAYER_INTEGRAL_UNTIL, LAYER3, SET_LAYER_INTEGRAL, &
    ADD, DERIVATIVE, RFLUX, SHIFT, ONE_PER_DAY,LAYER4,LAYER5
  use mem_Param,  ONLY: p_poro, p_d_tot,p_d_tot_2, p_q10diff
  use mem_BenCO2Transport
  use mem_BenthicNutrient3, ONLY:p_max_shift_change,p_max_state_change


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:GetInfoFromSet, &
  ! InitializeSet, DefineSet, CompleteSet, CalculateSet, CalculateTau, &
  ! CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: GetInfoFromSet, InitializeSet, DefineSet, &
    CompleteSet, CalculateSet, CalculateTau, CalculateFromSet


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp,insw

!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: r
  real(RLEN)  :: alpha
  real(RLEN)  :: diff
  real(RLEN)  :: Tau
  real(RLEN)  :: cG3c
  real(RLEN)  :: zu
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: jG13G3c
  real(RLEN)  :: jG23G13c
  real(RLEN)  :: jG33G23c
  real(RLEN)  :: jG3O3c
  real(RLEN)  :: Dnew
  real(RLEN)  :: Dx
  real(RLEN)  :: Dy
  real(RLEN)  :: rrQ1c
  real(RLEN)  :: dummy

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations 
      ! from mgC/m2 --> umol/kg
      ! Diagnostic variable used to compute pH
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      DICae(BoxNumberXY_ben) = G3c(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)/( &
        p_p+ ONE)/( D1m(BoxNumberXY_ben))/MW_C/ERHO_Ben(BoxNumberXY_ben)*1000._RLEN
      DICan(BoxNumberXY_ben) = (G13c(BoxNumberXY_ben)+G23c(BoxNumberXY_ben)) &
                                   / p_poro(BoxNumberXY_ben)/(p_p+ ONE)/ &
        (p_d_tot- D1m(BoxNumberXY_ben))/MW_C/ERHO_Ben(BoxNumberXY_ben)*1000._RLEN

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization mmo/m2 --> mgC/m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      zu = ( rrBTo(BoxNumberXY_ben)*MW_C) / D1m(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Temperature Correction (diffusion)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff = p_diff* irrenh(BoxNumberXY_ben)* p_poro(BoxNumberXY_ben)* &
        eTq( ETW_Ben(BoxNumberXY_ben), p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Assume Negative Exponential Distribution of Part.Carb. according D6.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      alpha  =   ONE/ D6m(BoxNumberXY_ben)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! average Anoxic Mineralization in the anaerobic layers  (mgC /m3/d)
      zuD1 = ( rrATo(BoxNumberXY_ben)* MW_C)/ p_poro(BoxNumberXY_ben)/ IntegralExp( &
        -alpha, p_d_tot- D1m(BoxNumberXY_ben))
      ! Anoxic Mineralization at D2.m  (mgC /m3/d)
      zuD2  =   zuD1* exp( - alpha*( D2m(BoxNumberXY_ben)- D1m(BoxNumberXY_ben)))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize and input physical boundaries and forcing:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KCO2(BoxNumberXY_ben)  =   InitializeSet(  KCO2(BoxNumberXY_ben),  5,  14)
      Dx=(D1m(BoxNumberXY_ben) +D2m(BoxNumberXY_ben)) * 0.5_RLEN
      Dy=(D2m(BoxNumberXY_ben) +p_d_tot) * 0.5_RLEN

      call DefineSet( KCO2(BoxNumberXY_ben), LAYERS, LAYER1, LAYER2, D1m(BoxNumberXY_ben),Dx)
      call DefineSet( KCO2(BoxNumberXY_ben), LAYERS, LAYER3, LAYER4, D2m(BoxNumberXY_ben),Dy)

      call DefineSet( KCO2(BoxNumberXY_ben), DIFFUSION, FOR_ALL_LAYERS, 0, diff,dummy)

      call DefineSet( KCO2(BoxNumberXY_ben), POROSITY, FOR_ALL_LAYERS, 0, &
                                                       p_poro(BoxNumberXY_ben), dummy)

      call DefineSet( KCO2(BoxNumberXY_ben), ADSORPTION, FOR_ALL_LAYERS, 0, p_p,dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Give particular solution for all layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 13, QUADRATIC_TERM, dummy, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 14, LINEAR_TERM, dummy, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 21, ZERO_EXPONENTIAL_TERM, &
        -alpha, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 31, ZERO_EXPONENTIAL_TERM, &
        -alpha, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 34, LINEAR_TERM, dummy, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 35, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 41, ZERO_EXPONENTIAL_TERM, &
        -alpha, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 44, LINEAR_TERM, dummy, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 45, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 51, ZERO_EXPONENTIAL_TERM, &
        -alpha, dummy)
      call DefineSet( KCO2(BoxNumberXY_ben), DEFINE, 55, CONSTANT_TERM, dummy, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert boundary conditions:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !1-8 boundary conditions:
      call CompleteSet( KCO2(BoxNumberXY_ben), SET_CONTINUITY, FLAG, MASS, &
        dummy, dummy)

      !9th boundary condition:
      call CompleteSet( KCO2(BoxNumberXY_ben), SET_BOUNDARY, LAYER1, &
        EQUATION, ZERO, O3c_Ben(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! a11 / (labda * labda * diff) = a12 / (labda * labda * diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !10th boundary condition:
      r  =   exp( - alpha*( Dx - D1m(BoxNumberXY_ben)))
      call FixProportionCoeff(KCO2(BoxNumberXY_ben),21,31,ONE,r)

      !11-12 boundary condition:
      select case ( InitializeModel)
        case ( 0 )
          call CompleteSet( KCO2(BoxNumberXY_ben), SET_LAYER_INTEGRAL, &
            LAYER2, LAYER3, dummy, G13c(BoxNumberXY_ben))
          call CompleteSet( KCO2(BoxNumberXY_ben), SET_LAYER_INTEGRAL_UNTIL, &
            LAYER4, LAYER5, p_d_tot_2, G23c(BoxNumberXY_ben))
        case ( 1 )
          call CompleteSet( KCO2(BoxNumberXY_ben), INPUT_TERM, 21, PARAMETER, dummy, zuD1)
          r  =   exp( - alpha*( D2m(BoxNumberXY_ben)-Dx) )
          call FixProportionCoeff(KCO2(BoxNumberXY_ben),31,41,ONE,r)

      end select

      !13th boundary condition:
      r  =   exp( - alpha*( Dy- D2m(BoxNumberXY_ben)))
      call FixProportionCoeff(KCO2(BoxNumberXY_ben),41,51,ONE,r)

      !14th boundary condition:
      call CompleteSet( KCO2(BoxNumberXY_ben), INPUT_TERM, 13, PARAMETER, dummy, zu)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 1. Calculate for above defined set of boundary conditions
      !      the gradient of the nutrient
      ! 2. Replace last condition by an alternative new one
      ! 3. Calculate the value belongin by the alternative condition with
      !  the gradient calculated under 1.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cG3c = CalculateSet( KCO2(BoxNumberXY_ben), SET_LAYER_INTEGRAL, LAYER1, &
        LAYER1, dummy, ZERO)

      if ( InitializeModel== 0) then
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate adaptation time absolute and Delta() relative to Tau:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        Tau  =   CalculateTau(  ZERO,  diff,  p_p,  D1m(BoxNumberXY_ben))

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Estimate the average value of K4n during the next time step:
        ! Hence this value depend as well as on adaptation time,
        ! ''old'' value, and on ''equilibrium value''
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        cG3c = cG3c+( G3c(BoxNumberXY_ben)- cG3c)* IntegralExp( - LocalDelta/ &
          Tau, ONE)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Recalulate gradient, now using cG3c
        ! Complete first the alternative condition by inputting the value
        ! for G3c (cG3c)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        dummy = CalculateSet( KCO2(BoxNumberXY_ben), ADD, 0, 0, dummy, cG3c)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! flux at D1.m
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        jG13G3c = CalculateFromSet( KCO2(BoxNumberXY_ben), DERIVATIVE, &
          RFLUX, D1m(BoxNumberXY_ben), dummy)

        Dnew  =   D1m(BoxNumberXY_ben)+ LocalDelta* shiftD1m(BoxNumberXY_ben)
        jG13G3c = jG13G3c+ CalculateFromSet( KCO2(BoxNumberXY_ben), SHIFT, &
          LAYER1, D1m(BoxNumberXY_ben), Dnew)/ LocalDelta


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Damp for too large fluxes
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        call LimitShift(jG13G3c,G3c(BoxNumberXY_ben),G13c(BoxNumberXY_ben),p_max_shift_change)
        call flux(BoxNumberXY_ben, iiBen, ppG13c, ppG3c,   jG13G3c* insw(  jG13G3c) )
        call flux(BoxNumberXY_ben, iiBen, ppG3c, ppG13c, - jG13G3c* insw( -jG13G3c) )

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! flux at D2m
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        jG23G13c = CalculateFromSet( KCO2(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
          D2m(BoxNumberXY_ben), dummy)

        Dnew  =   D2m(BoxNumberXY_ben)+ LocalDelta* shiftD2m(BoxNumberXY_ben)
        jG23G13c = jG23G13c+ CalculateFromSet( KCO2(BoxNumberXY_ben), SHIFT, &
          LAYER3, D2m(BoxNumberXY_ben), Dnew)/ LocalDelta

        jG23G13c = jG23G13c - zuD2* &
           p_poro(BoxNumberXY_ben)* IntegralExp( -alpha, p_d_tot- D2m(BoxNumberXY_ben));

        call LimitShift(jG23G13c,G13c(BoxNumberXY_ben),G23c(BoxNumberXY_ben),p_max_shift_change)
        call flux(BoxNumberXY_ben, iiBen, ppG23c, ppG13c,   jG23G13c* insw(  jG23G13c) )
        call flux(BoxNumberXY_ben, iiBen, ppG13c, ppG23c, - jG23G13c* insw( -jG23G13c) )

        jG33G23c = CalculateFromSet( KCO2(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
          p_d_tot_2, dummy)
        call LimitChange(1,jG33G23c,G23c(BoxNumberXY_ben),p_max_state_change)
        call flux(BoxNumberXY_ben, iiBen, ppG23c, ppG23c, jG33G23c )


        jG3O3c = CalculateFromSet( KCO2(BoxNumberXY_ben), DERIVATIVE, &
          RFLUX, ZERO, ZERO)

        call LimitShift(jG3O3c,O3c_Ben(BoxNumberXY_ben)*Depth_Ben(BoxNumberXY_ben) ,&
                                                         G3c(BoxNumberXY_ben),p_max_state_change)

        call flux(BoxNumberXY_ben, iiBen, ppG3c, ppG3c, -jG3O3c )
        jbotO3c(BoxNumberXY_ben)=jbotO3c(BoxNumberXY_ben)+jG3O3c

      end if

  end do
#endif

  end subroutine BenCO2TransportDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
