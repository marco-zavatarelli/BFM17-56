#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenQ1Transport
!
! DESCRIPTION
!   Description of the DOM dynamics in the sediments
!
! !INTERFACE
  subroutine BenQ1TransportDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ONE,ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: Q11c, Q1c, Q11n, Q1n, Q11p, Q1p, D1m, D6m, D2m, D2STATE_BEN
  use mem, ONLY: ppQ11c, ppQ1c, ppQ11n, ppQ1n, ppQ11p, ppQ1p, &
    ppD1m, ppD6m, ppD2m, NO_BOXES_XY, &
    BoxNumberXY_ben, LocalDelta, KQ1, irrenh, &
    ETW_Ben, shiftD1m, ruHI, reHI, iiH1, iiH2, iiBen, iiPel, flux
#endif
  use constants, ONLY: LAYERS, LAYER1, LAYER2, &
    DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DEFINE, EXPONENTIAL_TERM, &
    CONSTANT_TERM, ZERO_EXPONENTIAL_TERM, LINEAR_TERM, SET_CONTINUITY, FLAG, &
    MASS, SET_BOUNDARY, DERIVATIVE, SET_LAYER_INTEGRAL, LAYER3, INPUT_TERM, &
    STANDARD, START_ADD_TERM, SET_LAYER_INTEGRAL_UNTIL, ADD, INPUT_ADD_TERM,&
    SHIFT, RFLUX, ONE_PER_DAY
  use mem_Param,  ONLY: p_poro, p_q10diff, p_clDxm, p_d_tot,p_small
  use mem_BenQ1Transport
  use mem_BenthicNutrient3, ONLY:p_max_shift_change
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet, &
    CalculateSet, CalculateTau, CalculateFromSet
  use global_interface,   ONLY: eTq
  use mem_globalfun,   ONLY: IntegralExp, insw
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
!
! !REVISION_HISTORY
!   by Piet Ruardij  *:0 at Sun Dec 04 23:09:55 CET 2005
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
  real(RLEN)  :: M
  real(RLEN)  :: a15
  real(RLEN)  :: r
  real(RLEN)  :: alpha
  real(RLEN)  :: gamma
  real(RLEN)  :: diff
  real(RLEN)  :: Tau
  real(RLEN)  :: cQ1c
  real(RLEN)  :: zu
  real(RLEN)  :: sQ1
  real(RLEN)  :: zuD1
  real(RLEN)  :: rQ11
  real(RLEN)  :: jQ11Q1c
  real(RLEN)  :: jQ1Q11c
  real(RLEN)  :: flow
  real(RLEN)  :: Dnew
  real(RLEN)  :: dummy

  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get total Net Benthic DOC (Q1.c)
      ! production/consumption in the oxic layer (m2 --> m3 porewater)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      sQ1 = max( 0.001_RLEN, ruHI(iiH1, BoxNumberXY_ben)/( p_small+ Q1c(BoxNumberXY_ben)))
      M  =   reHI(iiH1,BoxNumberXY_ben)/D1m(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)
!     sQ1 = max( 0.001D+00, ruHI(iiH1, BoxNumberXY_ben)/ &
!       D1m(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben))/( 1.0D-80+ Q1c(BoxNumberXY_ben))
!     M  =   reHI(iiH1,BoxNumberXY_ben)

      a15  =   M/ sQ1

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff = p_diff* irrenh(BoxNumberXY_ben)* p_poro(BoxNumberXY_ben)* &
        eTq( ETW_Ben(BoxNumberXY_ben), p_q10diff)
      gamma  =   sqrt(  sQ1/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate specific bacterial consumption rate in anoxic layers (limited)
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      alpha  =   ONE/ max(  p_clDxm,  D6m(BoxNumberXY_ben))
      rQ11 = ( reHI(iiH2,BoxNumberXY_ben)- ruHI(iiH2,BoxNumberXY_ben))/ &
                  p_poro(BoxNumberXY_ben)/( p_d_tot- D1m(BoxNumberXY_ben))
      zuD1 = max( 1.D-20, rQ11)/ IntegralExp( - alpha, p_d_tot- &
        D1m(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      KQ1(BoxNumberXY_ben) = InitializeSet( KQ1(BoxNumberXY_ben), N_layers, &
        N_coeff)

      call DefineSet( KQ1(BoxNumberXY_ben), LAYERS, LAYER1, &
        LAYER2, D1m(BoxNumberXY_ben), D2m(BoxNumberXY_ben))

      call DefineSet( KQ1(BoxNumberXY_ben), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet( KQ1(BoxNumberXY_ben), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY_ben), dummy)

      call DefineSet( KQ1(BoxNumberXY_ben), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! Q(z) = q13*z^2 + q14*z + q15
      ! 2nd layer:
      ! Q(z) = q21*exp(gamma*z) + q22*exp(-gamma*z)
      ! 3rd layer:
      ! Q(z) = q31*exp(gamma*z) + q32*exp(-gamma*z) `
      !    q32 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 11, EXPONENTIAL_TERM, gamma, &
        dummy)

      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 12, EXPONENTIAL_TERM,  &
        -gamma, dummy)
      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 21, ZERO_EXPONENTIAL_TERM, &
        -alpha, dummy)
      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 31, ZERO_EXPONENTIAL_TERM, &
         -alpha, dummy)
      call DefineSet( KQ1(BoxNumberXY_ben), DEFINE, 35, CONSTANT_TERM, dummy, dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !4
      call CompleteSet( KQ1(BoxNumberXY_ben), SET_CONTINUITY, FLAG, MASS, dummy)
      !5
      call CompleteSet( KQ1(BoxNumberXY_ben), SET_BOUNDARY, LAYER1, DERIVATIVE, &
        ZERO, value=ZERO)
      !6:
      call CompleteSet( KQ1(BoxNumberXY_ben), SET_LAYER_INTEGRAL, LAYER2, &
        LAYER3, dummy, value=Q11c(BoxNumberXY_ben))
      !7:
      r  =   exp( - alpha*( D1m(BoxNumberXY_ben)- D2m(BoxNumberXY_ben)))

      select case ( r> p_small)

        case( .FALSE. )
          call CompleteSet( KQ1(BoxNumberXY_ben), INPUT_TERM, 31, STANDARD, &
            dummy, value=ZERO)

        case( .TRUE. )
          call CompleteSet( KQ1(BoxNumberXY_ben), START_ADD_TERM, 31, STANDARD, &
            dummy, mfac=ONE/r)

          call CompleteSet( KQ1(BoxNumberXY_ben), INPUT_ADD_TERM, 21, STANDARD, &
            dummy, mfac=-ONE)

      end select

      !8:
      call CompleteSet( KQ1(BoxNumberXY_ben), INPUT_TERM, 15, STANDARD, dummy, &
        value=a15)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      cQ1c = CalculateSet( KQ1(BoxNumberXY_ben), SET_LAYER_INTEGRAL_UNTIL, &
        LAYER1, LAYER1, D1m(BoxNumberXY_ben), ZERO)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      Tau  =   CalculateTau(  sQ1,  diff,  p_p,  D1m(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of Q11 over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      cQ1c = cQ1c+( Q1c(BoxNumberXY_ben)- cQ1c)* IntegralExp( - LocalDelta/ &
        Tau, ONE)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cQ11c as new constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      dummy  =   CalculateSet(  KQ1(BoxNumberXY_ben),  ADD,  0,  0,  dummy,  cQ1c)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Start calculation of fluxes:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate new depth of the oxygen horizon and the flux related
      ! to the shifting
      ! Add the flux at D1.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      Dnew  =   D1m(BoxNumberXY_ben)+ LocalDelta* shiftD1m(BoxNumberXY_ben)
      flow = CalculateFromSet( KQ1(BoxNumberXY_ben), SHIFT, LAYER1, &
        D1m(BoxNumberXY_ben), Dnew)/ LocalDelta

      flow = flow+ CalculateFromSet( KQ1(BoxNumberXY_ben), DERIVATIVE, &
        RFLUX, D1m(BoxNumberXY_ben), dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Limit for too large fluxes
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      flow=flow*p_max_shift_change/(abs(flow/r)+p_max_shift_change);
      call LimitShift(flow,Q1c(BoxNumberXY_ben),Q11c(BoxNumberXY_ben),p_max_shift_change)

      jQ1Q11c=-flow*insw(-flow)
      jQ11Q1c= flow*insw( flow)
      call flux(BoxNumberXY_ben, iiBen, ppQ11c, ppQ1c, jQ11Q1c) 
      call flux(BoxNumberXY_ben, iiBen, ppQ1c, ppQ11c, jQ1Q11c )


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! One of the 2 fluxes between the some constituents is 0!
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux(BoxNumberXY_ben, iiBen, ppQ1n, ppQ11n, jQ1Q11c/ Q1c(BoxNumberXY_ben)* &
        Q1n(BoxNumberXY_ben) )
      call flux(BoxNumberXY_ben, iiBen, ppQ11n, ppQ1n, jQ11Q1c/ Q11c(BoxNumberXY_ben)* &
        Q11n(BoxNumberXY_ben) )

      call flux(BoxNumberXY_ben, iiBen, ppQ1p, ppQ11p, jQ1Q11c/ Q1c(BoxNumberXY_ben)* &
        Q1p(BoxNumberXY_ben) )
      call flux(BoxNumberXY_ben, iiBen, ppQ11p, ppQ1p, jQ11Q1c/ Q11c(BoxNumberXY_ben)* &
        Q11p(BoxNumberXY_ben) )


  end do

  end subroutine BenQ1TransportDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
