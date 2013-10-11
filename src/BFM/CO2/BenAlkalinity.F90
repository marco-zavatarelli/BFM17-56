#include "INCLUDE.h"
#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAlkalinity
!
! DESCRIPTION
!   Description of the anoxic diagenitic processes in the sediment
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenAlkalinityDynamics
!

#ifdef INCLUDE_BENCO2
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: G23h,G13h, G3h, D1m, Q1c, D6m, D2m, D2STATE_BEN
  use mem, ONLY: ppG23h, ppG13h, ppG3h, ppD1m, ppQ1c, ppD6m, ppD2m,Acae, Acan, &
    NO_BOXES_XY, ERHO_Ben,  &
    BoxNumberXY_ben, InitializeModel, LocalDelta, KALK, jbotO3h, &
    irrenh, ETW_Ben, jK4K3n,jG2K7o,rrATo, O3h_Ben, shiftD1m, shiftD2m, ruHI, iiH1, &
    Depth_Ben,iiBen, iiPel, flux
#endif
  use constants, ONLY: GET, LABDA_1, DIFFUSION, COEFFICIENT, &
    LAYERS, LAYER1, LAYER2, LAYER3,LAYER4,LAYER5, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, DOUBLE_DEFINE, EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, &
    LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, &
    EQUATION, INPUT_TERM, PARAMETER, START_ADD_TERM, ZERO_EXPONENTIAL_TERM, &
    INPUT_ADD_TERM, SET_LAYER_INTEGRAL_UNTIL, SET_LAYER_INTEGRAL, &
    ADD, DERIVATIVE, RFLUX, SHIFT, ONE_PER_DAY,INTEGRAL,STANDARD
  use mem_Param,  ONLY: p_poro, p_d_tot, p_d_tot_2,p_q10diff
  use mem_BenAlkalinity
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
  real(RLEN)  :: s
  real(RLEN)  :: loss
  real(RLEN)  :: sG3
  real(RLEN)  :: lambda
  real(RLEN)  :: gamma
  real(RLEN)  :: alpha
  real(RLEN)  :: diff
  real(RLEN)  :: Tau
  real(RLEN)  :: cG3h
  real(RLEN)  :: zu
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: jG3O3h
  real(RLEN)  :: jG13G3h
  real(RLEN)  :: jG23G13h
  real(RLEN)  :: jG33G23h
  real(RLEN)  :: Dnew
  real(RLEN)  :: Dx,Dy
  real(RLEN)  :: dummy

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore water Alkalinity
      ! from mmol eq/m2 --> umol eq/kg
      ! Diagnostic variable used to compute pH
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      Acae(BoxNumberXY_ben) = &
           G3h(BoxNumberXY_ben)/ &
           p_poro(BoxNumberXY_ben)/( &
           p_p+ ONE)/( D1m(BoxNumberXY_ben))/ &
           ERHO_Ben(BoxNumberXY_ben)*1000._RLEN
      Acan(BoxNumberXY_ben) = &
           (G23h(BoxNumberXY_ben)+&
           G13h(BoxNumberXY_ben))/ &
           p_poro(BoxNumberXY_ben)/( &
           p_p+ ONE)/( p_d_tot- D1m(BoxNumberXY_ben))/ &
           ERHO_Ben(BoxNumberXY_ben)*1000._RLEN

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! H+-loss  due to nitrification and deoxidation of OH-
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      loss=p_qnh * jK4K3n(BoxNumberXY_ben) + p_qoh * jG2K7o(BoxNumberXY_ben) 
      sG3 = max(0.01_RLEN,loss/G3h(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients describing ammonium in the oxic layer :
      ! 1. labda of the exponential curve
      ! 2. parameter of the nitrification term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Temperature Correction (diffusion)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY_ben)* p_poro(BoxNumberXY_ben)* &
        eTq( ETW_Ben(BoxNumberXY_ben), p_q10diff)

      lambda= sqrt(sG3/diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Assume Negative Exponential Distribution of Part.Carb. according D6.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   ONE/ D6m(BoxNumberXY_ben)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! anoxidation rate at interface D1

      zuD1 = p_qoh * rrATo(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)/ IntegralExp( & 
                                          -alpha, p_d_tot- D1m(BoxNumberXY_ben))
      zuD2  =   zuD1* exp( - alpha*( D2m(BoxNumberXY_ben)- D1m(BoxNumberXY_ben)))


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize and input physical boundaries and forcing:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KALK(BoxNumberXY_ben)  =   InitializeSet(  KALK(BoxNumberXY_ben),  5,  14)

      Dx=(D1m(BoxNumberXY_ben) +D2m(BoxNumberXY_ben)) * 0.5_RLEN
      Dy=(D2m(BoxNumberXY_ben) +p_d_tot) * 0.5_RLEN
      call DefineSet( KALK(BoxNumberXY_ben), LAYERS, LAYER1, &
                                         LAYER2, D1m(BoxNumberXY_ben),Dx) 
      call DefineSet( KALK(BoxNumberXY_ben), LAYERS, LAYER3, LAYER4, D2m(BoxNumberXY_ben), Dy)
      call DefineSet( KALK(BoxNumberXY_ben), DIFFUSION, FOR_ALL_LAYERS, 0, diff, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), POROSITY, FOR_ALL_LAYERS, 0, &
                                                          p_poro(BoxNumberXY_ben), dummy)
      call DefineSet( KALK(BoxNumberXY_ben), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Give particular solution for all layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet( KALK(BoxNumberXY_ben), DOUBLE_DEFINE, 11, &
                                              EXPONENTIAL_TERM, lambda, sG3)
      call DefineSet( KALK(BoxNumberXY_ben), DOUBLE_DEFINE, 12, &
                                              EXPONENTIAL_TERM, - lambda, sG3)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 21, ZERO_EXPONENTIAL_TERM, &
                                                                  -alpha, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 31, ZERO_EXPONENTIAL_TERM, &
                                                                  -alpha, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 34, LINEAR_TERM, dummy, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 35, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 41, ZERO_EXPONENTIAL_TERM, &
                                                                  -alpha, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 44, LINEAR_TERM, dummy, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 45, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 51, ZERO_EXPONENTIAL_TERM, &
                                                                  -alpha, dummy)
      call DefineSet( KALK(BoxNumberXY_ben), DEFINE, 55, CONSTANT_TERM, dummy, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert boundary conditions:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !1-8 boundary conditions:
      call CompleteSet( KALK(BoxNumberXY_ben), SET_CONTINUITY, FLAG, MASS, dummy, dummy)

      !9th boundary condition:
      call CompleteSet( KALK(BoxNumberXY_ben), SET_BOUNDARY, LAYER1, &
                                               EQUATION, ZERO, O3h_Ben(BoxNumberXY_ben))

      !10-11th boundary condition:
      select case  ( InitializeModel ) 
        case (0)
          call CompleteSet( KALK(BoxNumberXY_ben), SET_LAYER_INTEGRAL, &
                                            LAYER2, LAYER3, dummy, G13h(BoxNumberXY_ben))
          call CompleteSet( KALK(BoxNumberXY_ben), SET_LAYER_INTEGRAL_UNTIL, &
                                         LAYER4, LAYER5, p_d_tot_2, G23h(BoxNumberXY_ben))
        case(1)
          r  =   exp( - alpha*( Dx - D1m(BoxNumberXY_ben)))
          call FixProportionCoeff(KALK(BoxNumberXY_ben),21,31,ONE,r)
          r  =   exp( - alpha*( D2m(BoxNumberXY_ben)-Dx) )
          call FixProportionCoeff(KALK(BoxNumberXY_ben),31,41,ONE,r)
      end select

      !12th boundary condition:
      r  =   exp( - alpha*( Dy- D2m(BoxNumberXY_ben)))
      call FixProportionCoeff(KALK(BoxNumberXY_ben),41,51,ONE,r)

      !13th boundary condition:
      call CompleteSet( KALK(BoxNumberXY_ben), INPUT_TERM, 21, PARAMETER, dummy, zuD1)

      r=-loss /D1m(BoxNumberXY_ben)/p_poro(BoxNumberXY_ben)
      call CompleteSet( KALK(BoxNumberXY_ben), INPUT_TERM, 15, STANDARD, &
             dummy, value=r/sG3)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! 1. Calculate for above defined set of boundary conditions
      !      the gradient of the nutrient
      ! 2. Replace last condition by an alternative new one
      ! 3. Calculate the value belongin by the alternative condition with
      !  the gradient calculated under 1.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cG3h = CalculateSet( KALK(BoxNumberXY_ben), SET_LAYER_INTEGRAL, LAYER1, &
        LAYER1, dummy, ZERO)

      if ( InitializeModel== 0) then
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate adaptation time absolute and Delta() relative to Tau:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        Tau  =   CalculateTau(  sG3,  diff,  p_p,  D1m(BoxNumberXY_ben))

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Estimate the average value of K4n during the next time step:
        ! Hence this value depend as well as on adaptation time,
        ! ''old'' value, and on ''equilibrium value''
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        cG3h = cG3h+( G3h(BoxNumberXY_ben)- cG3h)* IntegralExp( - LocalDelta/ &
          Tau, ONE)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Recalulate gradient, now using cG3h
        ! Complete first the alternative condition by inputting the value
        ! for G3h (cG3h)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        dummy = CalculateSet( KALK(BoxNumberXY_ben), ADD, 0, 0, dummy, cG3h)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! loss and gain terms
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jG3O3h = CalculateFromSet(KALK(BoxNumberXY_ben),DERIVATIVE,RFLUX, ZERO, ZERO) 
        s= sG3* CalculateFromSet( KALK(BoxNumberXY_ben), INTEGRAL, &
             RFLUX, ZERO, D1m(BoxNumberXY_ben))

        r=-s-jG3O3h

        call LimitChange(2,r,G3h(BoxNumberXY_ben),p_max_state_change)
        r=-r/(s+jG3O3h)
!       if ( loss+jG3o3h .gt.G3h(BoxNumberXY_ben)) then
!         write(LOGUNIT,*), '>>>',loss+jG3O3h,G3h(BoxNumberXY_ben),r
!       endif
        s=r*s;jG3O3h=jG3O3h*r;
        call flux(BoxNumberXY_ben, iiBen, ppG3h, ppG3h,  -s  )
        call flux(BoxNumberXY_ben, iiBen, ppG3h, ppG3h, -jG3O3h)
        ! to keep a closed budget : all H+ wich cannot removed in the benthic will be removed to pelagic
        jbotO3h(BoxNumberXY_ben)=jbotO3h(BoxNumberXY_ben)+jG3O3h+(loss-s) 

        r= zuD2 * IntegralExp( -alpha, p_d_tot- D2m(BoxNumberXY_ben))
        call flux(BoxNumberXY_ben, iiBen, ppG13h, ppG13h,  p_qoh * rrATo(BoxNumberXY_ben)-r)
        call flux(BoxNumberXY_ben, iiBen, ppG23h, ppG23h,  r)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! flux at D1m
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        jG13G3h = CalculateFromSet( KALK(BoxNumberXY_ben), DERIVATIVE, &
          RFLUX, D1m(BoxNumberXY_ben), dummy)

        Dnew  =   D1m(BoxNumberXY_ben)+ LocalDelta* shiftD1m(BoxNumberXY_ben)
        jG13G3h = jG13G3h+ CalculateFromSet( KALK(BoxNumberXY_ben), SHIFT, &
          LAYER1, D1m(BoxNumberXY_ben), Dnew)/ LocalDelta

        ! Damp for too large fluxes

        call LimitShift(jG13G3h, &
             G3h(BoxNumberXY_ben)- &
             loss-jG3O3h, &
             G13h(BoxNumberXY_ben), &
             p_max_shift_change)
        call flux(BoxNumberXY_ben, iiBen, ppG13h, ppG3h,   jG13G3h* insw(  jG13G3h) )
        call flux(BoxNumberXY_ben, iiBen, ppG3h, ppG13h, - jG13G3h* insw( -jG13G3h) )


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! flux at D2m
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        jG23G13h =  CalculateFromSet( KALK(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
          D2m(BoxNumberXY_ben), dummy)

        Dnew  =   D2m(BoxNumberXY_ben)+ LocalDelta* shiftD2m(BoxNumberXY_ben)
        jG23G13h = jG23G13h+ CalculateFromSet( KALK(BoxNumberXY_ben), SHIFT, &
          LAYER3, D2m(BoxNumberXY_ben), Dnew)/ LocalDelta

        ! Damp for too large fluxes

        call LimitShift(jG23G13h, &
             G13h(BoxNumberXY_ben), &
             G23h(BoxNumberXY_ben), &
             p_max_shift_change)
        call flux(BoxNumberXY_ben, iiBen, ppG23h, ppG13h,   jG23G13h* insw(  jG23G13h) )
        call flux(BoxNumberXY_ben, iiBen, ppG13h, ppG23h, - jG23G13h* insw( -jG23G13h) )

        jG33G23h = CalculateFromSet( KALK(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
          p_d_tot_2, dummy)

        call LimitChange(1,jG33G23h,G23h(BoxNumberXY_ben),p_max_state_change)
        call flux(BoxNumberXY_ben, iiBen, ppG23h, ppG23h, jG33G23h )

      end if

  end do
#endif

  end subroutine BenAlkalinityDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
