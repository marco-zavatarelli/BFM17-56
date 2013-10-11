#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenNitrate
!
! DESCRIPTION
!   Description of the diagenetic nitrate processes in the sediment
!       Details on the equations and the method used to calculate
!       the equilibrium and transient profiles can be found in
!       Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483   
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenNitrateDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: K3n, G4n, D2m, D6m, D1m, K16r,K26r, D2STATE_BEN
  use mem, ONLY: ppK3n, ppG4n, ppD2m, ppD6m, ppD1m, ppK16r, &
    NO_BOXES_XY,   &
    BoxNumberXY_ben, InitializeModel, LocalDelta, KNO3, jbotN3n, M3n, jK3G4n,&
    jK4K3n, rrATo, irrenh, ETW_Ben, KNH4, N3n_Ben, iiBen, iiPel, flux,Depth_Ben
#endif
  use constants, ONLY: GET, LABDA_1, LABDA_2, &
    COEFFICIENT, LAYERS, LAYER1, DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, DOUBLE_DEFINE, ZERO_EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, &
    LINEAR_TERM, CONSTANT_TERM, EXPONENTIAL_TERM, SET_CONTINUITY, FLAG, MASS, &
    SET_BOUNDARY, EQUATION, INPUT_TERM, PARAMETER, START_ADD_TERM, STANDARD, &
    SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, DERIVATIVE, RFLUX, INPUT_ADD_TERM,&
    INTEGRAL, MIN_VAL_EXPFUN,LAYER3
  use mem_Param,  ONLY: p_poro, p_d_tot, p_q10diff, p_qro, p_qon_dentri,p_small
  use mem_BenNitrate
  use mem_BenthicNutrient3, ONLY:p_max_state_change


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
! !REVISION_HISTORY
!   September 1999 by M. Vichi  Commented version 
!
!
! COPYING
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
  real(RLEN)  :: Dxm
  real(RLEN)  :: sK4K3
  real(RLEN)  :: sK3G4
  real(RLEN)  :: diff
  real(RLEN)  :: gamma
  real(RLEN)  :: labda
  real(RLEN)  :: a11
  real(RLEN)  :: a12
  real(RLEN)  :: a15
  real(RLEN)  :: n12
  real(RLEN)  :: cK3n
  real(RLEN)  :: Tau
  real(RLEN)  :: zATo
  real(RLEN)  :: alpha
  real(RLEN)  :: s
  real(RLEN)  :: r
  real(RLEN)  :: dummy

  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M3n(BoxNumberXY_ben) = K3n(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)/( p_p+ &
        ONE)/( D2m(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D6.m is the average penetration depth for C-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      alpha  =   ONE/ D6m(BoxNumberXY_ben)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      ! Calculate the total anoxic mineralization in mmol O/m3/d
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      zATo = rrATo(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)/ IntegralExp( &
        -alpha, p_d_tot- D1m(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! denitrification: temperature and coupling with anoxic mineralization
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff = p_diff* irrenh(BoxNumberXY_ben)* p_poro(BoxNumberXY_ben)* &
        eTq( ETW_Ben(BoxNumberXY_ben), p_q10diff)
      sK3G4  =   p_sK3G4* zATo/ p_zATo

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Technical correction when K6r is nearly zero. This denitrification
      ! is limited because not enough red. equiv. are present to oxidize &
      ! material
      ! Calculate net consumption of reduction equivalents and limit &
      ! denitrifaction rate
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if ( InitializeModel== 0) then
        s=K16r(BoxNumberXY_ben)+K26r(BoxNumberXY_ben)
        r = max( p_small, -rrATo(BoxNumberXY_ben)* p_qro &
             +sK3G4* p_qro* p_qon_dentri* 0.5_RLEN * K3n(BoxNumberXY_ben))
        sK3G4  = max( 0.001_RLEN, sK3G4* s/( r+s))
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      gamma  =   sqrt(  sK3G4/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients describing ammonium in the oxic layer :
      ! 1. labda of the exponential curve
      ! 2. parameter of the nitrification term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      labda = GetInfoFromSet( KNH4(BoxNumberXY_ben), GET, LABDA_1, 11)
      sK4K3 = GetInfoFromSet( KNH4(BoxNumberXY_ben), GET, LABDA_2, 12)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients of all terms of equation valid for the
      ! first layer of ammonium (integration constants)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      a11 = GetInfoFromSet( KNH4(BoxNumberXY_ben), GET, COEFFICIENT, 11)
      a12 = GetInfoFromSet( KNH4(BoxNumberXY_ben), GET, COEFFICIENT, 12)
      a15 = GetInfoFromSet( KNH4(BoxNumberXY_ben), GET, COEFFICIENT, 15)
      n12  =   sK4K3* a12

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      Dxm=D1m(BoxNumberXY_ben)+D2m(BoxNumberXY_ben)
      KNO3(BoxNumberXY_ben) = InitializeSet( KNO3(BoxNumberXY_ben), 3, 8)
      call DefineSet(KNO3(BoxNumberXY_ben), LAYERS, LAYER1, LAYER2, &
        D1m(BoxNumberXY_ben), Dxm)

      call DefineSet(KNO3(BoxNumberXY_ben), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet(KNO3(BoxNumberXY_ben), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY_ben), dummy)

      call DefineSet( KNO3(BoxNumberXY_ben), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, &
        dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! N(z) = n11*exp(labda*z) + n12*exp(-labda*z) + n13*z^2 + n14*z + n15
      ! 2nd layer:
      ! N(z) = n21*exp(gamma*z) + n22*exp(-gamma*z)
      !    n22 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call DefineSet( KNO3(BoxNumberXY_ben), DOUBLE_DEFINE, 11, &
        ZERO_EXPONENTIAL_TERM, labda, sK4K3)

      call DefineSet( KNO3(BoxNumberXY_ben), DOUBLE_DEFINE, 12, &
        ZERO_EXPONENTIAL_TERM, - labda, sK4K3)

      call DefineSet( KNO3(BoxNumberXY_ben), DEFINE, 13, QUADRATIC_TERM, dummy, &
        dummy)

      call DefineSet(KNO3(BoxNumberXY_ben), DEFINE, 14, LINEAR_TERM, dummy,dummy)
      call DefineSet(KNO3(BoxNumberXY_ben), DEFINE, 15, CONSTANT_TERM, dummy,dummy)

      call DefineSet( KNO3(BoxNumberXY_ben), DOUBLE_DEFINE, 21, EXPONENTIAL_TERM, - &
        gamma, sK3G4)

      call DefineSet( KNO3(BoxNumberXY_ben), DOUBLE_DEFINE, 22, EXPONENTIAL_TERM, &
        gamma, sK3G4)

      call DefineSet( KNO3(BoxNumberXY_ben), DOUBLE_DEFINE, 31, EXPONENTIAL_TERM, - &
        gamma, sK3G4)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call CompleteSet( KNO3(BoxNumberXY_ben), SET_CONTINUITY, FLAG, MASS, dummy)

      call CompleteSet( KNO3(BoxNumberXY_ben), SET_BOUNDARY, LAYER1, &
        EQUATION, ZERO, value=N3n_Ben(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! a11 / (labda * labda * diff) = a12 / (labda * labda * diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call FixProportionCoeff(KNO3(BoxNumberXY_ben),12,11,a12,a11)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! a11 / (labda * labda * diff) = a15 / (2 * diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call FixProportionCoeff(KNO3(BoxNumberXY_ben),11,13,a11,a15)

      if ( InitializeModel== 0) then
        call CompleteSet( KNO3(BoxNumberXY_ben), INPUT_TERM, 12, PARAMETER, dummy, &
           value=n12)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        cK3n = CalculateSet( KNO3(BoxNumberXY_ben), SET_LAYER_INTEGRAL_UNTIL, &
          LAYER1, LAYER3, D2m(BoxNumberXY_ben), ZERO)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate the adaptation time to the steady-state profile
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        Tau  =   CalculateTau(  ZERO,  diff,  p_p,  D2m(BoxNumberXY_ben))

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Estimate the average value of K3n over the actual time step
        ! (transient value).
        ! This value depends on the adaptation time, the actual time step,
        ! the ''old'' value and the ''equilibrium value''
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        cK3n = cK3n+( K3n(BoxNumberXY_ben)- cK3n)* IntegralExp( - LocalDelta/ &
          Tau, ONE)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Derive the equations for the transient profiles, assuming the same
        ! solution as for the steady-state case and using cK3n as new &
        ! constraint.
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        dummy = CalculateSet( KNO3(BoxNumberXY_ben), ADD, 0, 0, dummy, &
          cK3n)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Start calculation of fluxes:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Flux at the lower boundary
        ! All transport between layers are done in BenNitrogenShifting:
        !
        ! jK13K3n = CalculateFromSet(KNO3, DERIVATIVE, RFLUX, D2.m, dummy);
        !
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Nitrification is already calculated by BenAmmonium:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        !jK4K3n = sK4K3* CalculateFromSet( KNH4(BoxNumberXY_ben), INTEGRAL, &
        !    RFLUX, ZERO, D1m(BoxNumberXY_ben))

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Denitrification:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jK3G4n(BoxNumberXY_ben) = CalculateFromSet( KNO3(BoxNumberXY_ben), INTEGRAL, &
          RFLUX, D1m(BoxNumberXY_ben), D2m(BoxNumberXY_ben))* sK3G4
        if (  jK3G4n(BoxNumberXY_ben) > ZERO ) then
           call LimitChange(1,jK3G4n(BoxNumberXY_ben),K3n(BoxNumberXY_ben),p_max_state_change)
           call flux(BoxNumberXY_ben, iiBen, ppK3n, ppG4n, jK3G4n(BoxNumberXY_ben) )
        else
           call PrintSet(KNO3(BoxNumberXY_ben),"Negative (or Nan) flux for benthic nitrate")
           write(LOGUNIT,'(''D1m='',F10.3)') D1m(BoxNumberXY_ben)
           write(LOGUNIT,'(''D2m='',F10.3)') D2m(BoxNumberXY_ben)
           write(LOGUNIT,'(''nitrate/m2 (K3n)='',F10.3)') K3n(BoxNumberXY_ben)
        endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Flux at the water/sediment interface
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jbotN3n(BoxNumberXY_ben) = CalculateFromSet( KNO3(BoxNumberXY_ben), DERIVATIVE, &
          RFLUX, ZERO, dummy)

        call LimitShift(jbotN3n(BoxNumberXY_ben), &
             N3n_Ben(BoxNumberXY_ben)*Depth_Ben(BoxNumberXY_ben) ,&
             K3n(BoxNumberXY_ben),p_max_state_change)
        jbotN3n(BoxNumberXY_ben)=min(jbotN3n(BoxNumberXY_ben), K3n(BoxNumberXY_ben)  &
                  * p_max_state_change+ jK4K3n(BoxNumberXY_ben)-jK3G4n(BoxNumberXY_ben))
        call flux(BoxNumberXY_ben, iiBen, ppK3n, ppK3n, -( jbotN3n(BoxNumberXY_ben)) )

      else
        call CompleteSet( KNO3(BoxNumberXY_ben), INPUT_TERM, 12, PARAMETER, dummy, &
           value=n12)
        dummy = CalculateSet( KNO3(BoxNumberXY_ben), 0, 0,  0, dummy, dummy)
      end if

  end do

  end subroutine BenNitrateDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
