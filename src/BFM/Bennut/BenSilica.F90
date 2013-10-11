#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenSilica
!
! DESCRIPTION
!   Description of the diagenitic processes in the sediment
!       Details on the equations and the method used to calculate
!       the equilibrium and transient profiles can be found in
!       Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
!
! !INTERFACE
  subroutine BenSilicaDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: K5s, Q6s, D9m, D1m, D2m, D2STATE_BEN
  use mem, ONLY: ppK5s, ppQ6s, ppD9m, ppD1m, ppD2m, &
    NO_BOXES_XY,   &
    BoxNumberXY_ben, InitializeModel, LocalDelta, M5s, KSIO3, Depth_Ben, &
    KSIO3E, jbotN5s, jK15K5s, irrenh, ETW_Ben, N5s_Ben, shiftD2m, iiBen, iiPel, flux
#endif
  use constants, ONLY: LAYERS, LAYER1, DIFFUSION, &
    FOR_ALL_LAYERS, POROSITY, ADSORPTION, DEFINE, QUADRATIC_TERM, LINEAR_TERM, &
    CONSTANT_TERM, PARAMETER_DEFINE, BESSELI_EXP_TERM, SET_CONTINUITY, STANDARD, &
    SET_BOUNDARY, EQUATION, INPUT_TERM, PARAMETER, SET_LAYER_INTEGRAL_UNTIL, &
    LAYER2, ADD, INTEGRAL, DERIVATIVE, RFLUX, MASS, EXPONENTIAL_INTEGRAL
  use mem_Param,  ONLY: p_poro, p_clD1D2m, p_q10diff, p_clDxm, p_d_tot,p_small
  use mem_BenSilica
  use mem_BenthicNutrient3, ONLY:p_max_state_change,p_max_shift_change

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CopySet, &
  ! CalculateFromSet, GetInfoFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet, &
    CalculateSet, CalculateTau, CopySet, CalculateFromSet, GetInfoFromSet

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp, insw
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! !REVISION_HISTORY
!       September 1999 by M. Vichi     Commented version
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij & M. Vichi
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
  real(RLEN)  :: cD1m
  real(RLEN)  :: cD2m
  real(RLEN)  :: cD2mNew
  real(RLEN)  :: cShiftD2m
  real(RLEN)  :: chM5s
  real(RLEN)  :: cM5s
  real(RLEN)  :: cmm
  real(RLEN)  :: Tau
  real(RLEN)  :: alpha
  real(RLEN)  :: diff
  real(RLEN)  :: M5b0
  real(RLEN)  :: M5b_0_d1
  real(RLEN)  :: M5bD1
  real(RLEN)  :: zuBT
  real(RLEN)  :: suD1
  real(RLEN)  :: rmQ6s
  real(RLEN)  :: shiftmass
  real(RLEN)  :: jQ6K5s
  real(RLEN)  :: jQ6K15s
  real(RLEN)  :: smQ6
  real(RLEN)  :: dummy
  integer     :: idummy

  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! saturation value: temperature
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      chM5s = p_chM5s+ p_cvM5s*( eTq( ETW_Ben(BoxNumberXY_ben), p_q10)- &
        ONE)

      cD1m  =   min(  max(  D1m(BoxNumberXY_ben),   p_clD1m),  p_chD2m- p_clD1D2m)
      cD2m  =   min(  max(  D2m(BoxNumberXY_ben),   p_clD2m),  p_chD2m)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! Here M5s is used in the calculation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M5s(BoxNumberXY_ben)= chM5s- K5s(BoxNumberXY_ben)/cD2m &
                                 / ( 1.0+p_p) / p_poro(BoxNumberXY_ben)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! dissolution rate: temperature
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff = p_diff* irrenh(BoxNumberXY_ben)* p_poro(BoxNumberXY_ben)* &
        eTq( ETW_Ben(BoxNumberXY_ben), p_q10diff)
      smQ6  =   p_smQ6* eTq(  ETW_Ben(BoxNumberXY_ben),  p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D9.m is the average penetration depth for biogenic Si
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      alpha  =   ONE/ max(  p_clDxm,  D9m(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate total biogenic silica from m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M5b0 = Q6s(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)/ IntegralExp( - alpha, &
        p_d_tot)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Average content of Biogenic silica in the oxic layer
      ! and calculation of the zero-order dissolution term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M5b_0_d1  =   M5b0* IntegralExp( - alpha,  cD1m)/ cD1m
      zuBT  =   smQ6* M5b_0_d1

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Biogenic silica at cD1m and calculation of the dissolution rate
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M5bD1  =   M5b0* exp( - alpha* cD1m)
      suD1  =   smQ6* M5bD1/ chM5s

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      KSIO3(BoxNumberXY_ben) = InitializeSet( KSIO3(BoxNumberXY_ben), N_layers, N_coeff)

      call  DefineSet(  KSIO3(BoxNumberXY_ben),  LAYERS,  LAYER1,  0,  cD1m,  dummy)

      call DefineSet( KSIO3(BoxNumberXY_ben), DIFFUSION, FOR_ALL_LAYERS, idummy, &
        diff, dummy)

      call DefineSet( KSIO3(BoxNumberXY_ben), POROSITY, FOR_ALL_LAYERS, &
        idummy, p_poro(BoxNumberXY_ben), dummy)

      call DefineSet( KSIO3(BoxNumberXY_ben), ADSORPTION, FOR_ALL_LAYERS, idummy, &
        p_p, dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! C = Ssat - S
      !
      ! 1st layer:
      ! C(z) = c13*z^2 + c14*z + c15
      ! 2nd layer:
      ! C(z) = c21*I0*exp[-alpha*(z-cD1m
      !    I0 = modified Bessel function of 0-order
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call DefineSet( KSIO3(BoxNumberXY_ben), DEFINE, 13, QUADRATIC_TERM, dummy, dummy)
      call DefineSet( KSIO3(BoxNumberXY_ben), DEFINE, 14, LINEAR_TERM, dummy, dummy)
      call DefineSet( KSIO3(BoxNumberXY_ben), DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KSIO3(BoxNumberXY_ben), PARAMETER_DEFINE, 21, &
        BESSELI_EXP_TERM, - alpha, suD1)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call CompleteSet( KSIO3(BoxNumberXY_ben), SET_CONTINUITY, STANDARD, idummy, dummy)

      call CompleteSet( KSIO3(BoxNumberXY_ben), SET_BOUNDARY, LAYER1, &
        EQUATION, ZERO, value=chM5s- N5s_Ben(BoxNumberXY_ben))


      if ( InitializeModel== 0) then
         call CompleteSet( KSIO3(BoxNumberXY_ben), INPUT_TERM, 13, PARAMETER, dummy, &
                value=- zuBT)
      else
         call CompleteSet( KSIO3(BoxNumberXY_ben), INPUT_TERM, 13, PARAMETER, dummy, &
                value=- zuBT)
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the average concentration
      ! in the oxic and denitrification layers
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      cM5s = CalculateSet( KSIO3(BoxNumberXY_ben), SET_LAYER_INTEGRAL_UNTIL, LAYER1, &
        LAYER2, cD2m, ZERO)/ cD2m

      if ( InitializeModel== 0) then
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate the adaptation time to the steady-state profile
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        Tau  =   CalculateTau(  ZERO,  diff,  p_p,  cD2m)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Estimate the average value of M5s over the actual time step
        ! (transient value).
        ! This value depends on the adaptation time, the actual time step,
        ! the ''old'' value and the ''equilibrium value''
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        cM5s = cM5s+( M5s(BoxNumberXY_ben)- cM5s)* IntegralExp( - &
          LocalDelta/ Tau, ONE)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! 1.Store equilibrium profile
        ! 2.Derive the equations for the transient profiles, assuming the same
        ! solution as for the steady-state case and using cM5s as new &
        ! constraint.
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        KSIO3E(BoxNumberXY_ben) = CopySet( KSIO3(BoxNumberXY_ben), &
          KSIO3E(BoxNumberXY_ben))
        dummy = CalculateSet( KSIO3(BoxNumberXY_ben), ADD, 0, 0, dummy, cD2m* cM5s)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Recalculate the pore-water average concentrations for the standard &
        ! ''D2.n''
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Start calculation of fluxes:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate flux at the sediment/water interface:
        ! Flux limitation at very low values of N5s
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jbotN5s(BoxNumberXY_ben) = - CalculateFromSet( KSIO3(BoxNumberXY_ben), &
          DERIVATIVE, RFLUX, ZERO, dummy)

        call LimitShift(jbotN5s(BoxNumberXY_ben), &
             N5s_Ben(BoxNumberXY_ben)*Depth_Ben(BoxNumberXY_ben) ,&
             K5s(BoxNumberXY_ben),p_max_state_change)
        call flux(BoxNumberXY_ben, iiBen, ppK5s, ppK5s, -jbotN5s(BoxNumberXY_ben) )

        jK15K5s(BoxNumberXY_ben) = - CalculateFromSet( KSIO3(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
          cD2m, dummy)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate new depth of the sulphide horizon
        ! and the flux of silicate related to this shifting
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        shiftmass=ZERO
        if ( abs(shiftD2m(BoxNumberXY_ben))> ZERO) then

          cD2mNew = min( p_chD2m, max( p_clD2m, &
            D2m(BoxNumberXY_ben)+ shiftD2m(BoxNumberXY_ben)* LocalDelta))

          cShiftD2m  =   cD2mNew- cD2m

          if ( abs(cShiftD2m)> ZERO) then
            shiftmass = (chM5s* cShiftD2m* p_poro(BoxNumberXY_ben)*( ONE+ p_p) &
             - CalculateFromSet( KSIO3(BoxNumberXY_ben), INTEGRAL, MASS,cD2m, cD2mNew)) &
             /LocalDelta
          end if
        endif

        jK15K5s(BoxNumberXY_ben)  =   jK15K5s(BoxNumberXY_ben)+ shiftmass

        call LimitChange(2, &
             jK15K5s(BoxNumberXY_ben), &
             K5s(BoxNumberXY_ben), &
             p_max_shift_change)


        call flux(BoxNumberXY_ben, iiBen, ppK5s, ppK5s, jK15K5s(BoxNumberXY_ben) )

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! the dissolution fluxes:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jQ6K15s = max( ZERO, suD1* CalculateFromSet( &
          KSIO3E(BoxNumberXY_ben), EXPONENTIAL_INTEGRAL, RFLUX, cD2m, p_d_tot))
        call flux(BoxNumberXY_ben, iiBen, ppQ6s, ppQ6s, -( jQ6K15s) )

        jQ6K5s = max( ZERO, suD1* CalculateFromSet( &
          KSIO3E(BoxNumberXY_ben), EXPONENTIAL_INTEGRAL, RFLUX, cD1m, cD2m))

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Determine dissolution rate in oxidized layer (mMol/m3).
        ! Zero order process:
        ! Maximalization: this important source can not cause higher values
        ! than equilibrium flux of jQ6M5s is limited in such a way that M5s can
        ! never reach a value higher than the chM5s:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        rmQ6s = - GetInfoFromSet( KSIO3E(BoxNumberXY_ben), INTEGRAL, PARAMETER, &
          13, at_x=ZERO, to_x=cD1m)

        rmQ6s = max( ZERO, min( &
          M5s(BoxNumberXY_ben)* cD2m* p_poro(BoxNumberXY_ben)*( ONE+ p_p) &
          + jbotN5s(BoxNumberXY_ben)- jK15K5s(BoxNumberXY_ben)- jQ6K5s, rmQ6s))

        call flux(BoxNumberXY_ben, iiBen, ppQ6s, ppK5s, jQ6K5s+ rmQ6s )

        ! Determine where the median is of Q6 in the range from clm to D1m
         cmm  =  -log(0.5_RLEN*(ONE+exp(- cD2m/D9m(BoxNumberXY_ben))))*D9m(BoxNumberXY_ben)*0.5_RLEN

         call flux(BoxNumberXY_ben, iiBen, ppD9m, ppD9m, ( cmm &
            - D9m(BoxNumberXY_ben))*( jQ6K5s+rmQ6s)/( p_small+ Q6s(BoxNumberXY_ben)) )

        ! Determine where the median is of Q6 in the range from clm to D1m
         cmm  =   cD2m-log(0.5_RLEN*(ONE+exp(- (p_d_tot-cD2m)/D9m(BoxNumberXY_ben)))) &
                                                       *D9m(BoxNumberXY_ben)*0.5_RLEN
        call flux(BoxNumberXY_ben, iiBen, ppD9m, ppD9m, ( cmm  &
            - D9m(BoxNumberXY_ben))*( jQ6K15s)/( p_small+Q6s(BoxNumberXY_ben)) )

        M5s(BoxNumberXY_ben) = max(ZERO,chM5s- CalculateFromSet( KSIO3E(BoxNumberXY_ben), &
          INTEGRAL, STANDARD, ZERO, D2m(BoxNumberXY_ben))/ D2m(BoxNumberXY_ben))
      end if

  end do

  end subroutine BenSilicaDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
