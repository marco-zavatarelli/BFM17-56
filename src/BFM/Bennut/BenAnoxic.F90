#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAnoxic
!
! DESCRIPTION
!   Description of the anoxic diagenetic processes in the sediment
!       Details on the equations and the method used to calculate
!       the equilibrium and transient profiles can be found in
!       Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
! !INTERFACE
  subroutine BenAnoxicDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: K26r, K16r, K6r, G2o, D6m, D1m, D2m, D2STATE
  use mem, ONLY: ppK26r, ppK16r, ppK6r, ppG2o, ppD6m, ppD1m, dummy,  &
      NO_BOXES_XY,    &
    BoxNumberXY, LocalDelta, InitializeModel, M6r, KRED, jbotN6r, jG2K7o, rrATo, &
    rrBTo, jK36K26r, irrenh, ETW_Ben, KNO3, N6r_Ben,Depth_Ben, &
    shiftD1m, shiftD2m,iiBen, iiPel, flux,jK3G4n
#endif
  use constants, ONLY: GET, LABDA_1, LABDA_2, COEFFICIENT, SET_LAYER_INTEGRAL, &
    LAYERS, LAYER1, DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, LAYER5, STANDARD, &
    DEFINE, EXPONENTIAL_TERM, ZERO_EXPONENTIAL_TERM, DOUBLE_DEFINE,LINEAR_TERM, &
    CONSTANT_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, EQUATION, LAYER4, &
    INPUT_TERM, PARAMETER, SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, DERIVATIVE, &
    RFLUX, INTEGRAL, ONE_PER_DAY, START_ADD_TERM, INPUT_ADD_TERM, LAYER3, SHIFT
  use mem_Param,  ONLY: p_poro, p_d_tot,p_d_tot_2, p_clDxm, p_qro, p_q10diff, p_qon_dentri
  use mem_BenthicNutrient3, ONLY:p_max_state_change, p_max_shift_change
  use mem_BenAnoxic


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
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp, insw
!
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! !REVISION_HISTORY
!   September 1999 by M. Vichi Commented version
!
!
! COPYING
!
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij & M.VIchi
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
  real(RLEN)  :: gamma
  real(RLEN)  :: alpha
  real(RLEN)  :: zuD1
  real(RLEN)  :: diff
  real(RLEN)  :: labda
  real(RLEN)  :: sK3G4
  real(RLEN)  :: n21
  real(RLEN)  :: n22
  real(RLEN)  :: n31
  real(RLEN)  :: Tau
  real(RLEN)  :: cK6r
  real(RLEN)  :: jBTK6r
  real(RLEN)  :: jATK6r
  real(RLEN)  :: jK26K16r
  real(RLEN)  :: jK6BTr
  real(RLEN)  :: jK6G4r
  real(RLEN)  :: jK16K6r
  real(RLEN)  :: shiftmass
  real(RLEN)  :: Dnew
  real(RLEN)  :: r
  real(RLEN)  :: r2
  real(RLEN)  :: Dxm,Dym

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M6r(BoxNumberXY) = K6r(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ 1.0D+00)/( &
        D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D6.m is the average penetration depth for C-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   1.0D+00/ max(  p_clDxm,  D6m(BoxNumberXY))


      if ( InitializeModel == 0 ) then
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Convert anoxic mineralization (mmol S/m2/d)
        ! This rate is already assigned to the dynamical equation for K6.r
        ! in BenBacDyanmics for H2:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        jATK6r  =   p_qro* rrATo(BoxNumberXY)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Recalculate Mineralization m2 --> m3 porewater
        ! Anoxic mineralization at D1.m
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        zuD1 = jATK6r/ p_poro(BoxNumberXY)/ IntegralExp( - alpha, &
             p_d_tot- D1m(BoxNumberXY))
      else
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! In case of no info about anoxi mineralization at the start
        ! Reconstruct using detritus distribution alpha and oxic mineralization
        ! an anoxic mineralization at the upperside of the denitrification layer
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jBTK6r= p_qro * rrBTo(BoxNumberXY)
        zuD1 = jBTK6r / p_poro(BoxNumberXY)/   &
            IntegralExp( - alpha, D1m(BoxNumberXY)) *exp(-alpha * D1m(BoxNUmberXY))
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)* &
        eTq( ETW_Ben(BoxNumberXY), p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      gamma  =   sqrt(p_sOS/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients describing Nitrate in anoxic layer :
      ! 1. labda of the exponential curve, and the denitrification rate
      ! 2. parameter of the denitrification term (integration constant)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      labda = GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_1, 21)
      sK3G4 = GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_2, 21) &
                                               * p_qro* p_qon_dentri
      n21 = - GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 21)* sK3G4
      n22 = - GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 22)* sK3G4
      n31 = - GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 31)* sK3G4

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dxm=(D1m(BoxNumberXY)+D2m(BoxNumberXY)) *0.5
      Dym=(D2m(BoxNumberXY)+p_d_tot_2) *0.5

      KRED(BoxNumberXY) = InitializeSet( KRED(BoxNumberXY), 5, 16)
      call DefineSet(KRED(BoxNumberXY), LAYERS,LAYER1,LAYER2, D1m(BoxNumberXY),Dxm)
      call DefineSet(KRED(BoxNumberXY), LAYERS,LAYER3,LAYER4, D2m(BoxNumberXY),Dym)

      call DefineSet(KRED(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, dummy)
      call DefineSet(KRED(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet(KRED(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! R(z) = r11*exp(gamma*z) + r12*exp(-gamma*z)
      ! 2nd layer:
      ! R(z) = r21*exp[-alpha*(z-D1.m)] + r22*exp(-labda*z) + r23*z^2 + r24*z + &
      ! r25
      !    r23, r24 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet(KRED(BoxNumberXY), DEFINE, 11, EXPONENTIAL_TERM, -gamma, dummy)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 12, EXPONENTIAL_TERM, gamma, dummy)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet(KRED(BoxNumberXY), DEFINE, 21, ZERO_EXPONENTIAL_TERM,-alpha, dummy)
      call DefineSet(KRED(BoxNumberXY), DOUBLE_DEFINE, 22, &
                                           ZERO_EXPONENTIAL_TERM, labda, sK3G4)
      call DefineSet(KRED(BoxNumberXY), DOUBLE_DEFINE, 23, &
                                           ZERO_EXPONENTIAL_TERM, -labda, sK3G4)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      call DefineSet(KRED(BoxNumberXY), DEFINE, 31, ZERO_EXPONENTIAL_TERM, - &
        alpha, dummy)
      call DefineSet(KRED(BoxNumberXY), DOUBLE_DEFINE, 32, &
                                           ZERO_EXPONENTIAL_TERM, labda, sK3G4)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 35, CONSTANT_TERM, dummy, dummy)

      call DefineSet(KRED(BoxNumberXY), DEFINE, 41, ZERO_EXPONENTIAL_TERM,-alpha, dummy)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 44, LINEAR_TERM, dummy, dummy)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 45, CONSTANT_TERM, dummy, dummy)

      call DefineSet(KRED(BoxNumberXY), DEFINE, 51, ZERO_EXPONENTIAL_TERM,-alpha, dummy)
      call DefineSet(KRED(BoxNumberXY), DEFINE, 55, CONSTANT_TERM, dummy, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !1-8:
      call CompleteSet( KRED(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, dummy)

      !9:
      call CompleteSet( KRED(BoxNumberXY), SET_BOUNDARY, LAYER1, &
        EQUATION, 0.0D+00, value=N6r_Ben(BoxNumberXY))

      !10-11:
      call FixProportionCoeff(KRED(BoxNumberXY),22,32,n21,n31)
      call FixProportionCoeff(KRED(BoxNumberXY),22,23,n21,n22)

      select case (InitializeModel)

        !12 -13
        case(0)
          call CompleteSet(KRED(BoxNumberXY), SET_LAYER_INTEGRAL, &
                                   LAYER2, LAYER3,dummy ,     K16r(BoxNumberXY))
          call CompleteSet(KRED(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
                                   LAYER4, LAYER5, p_d_tot_2, K26r(BoxNumberXY))
        case(1)
           call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 31, PARAMETER, &
                dummy, value=zuD1*exp(-alpha *(Dxm-D1m(BoxNumberXY))))
           call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 41, PARAMETER, &
                dummy, value=zuD1*exp(-alpha * (D2m(BoxNumberXY)-D1m(BoxNumberXY))))
      end select

      !14
      r= exp(-alpha*(Dxm- D1m(BoxNUmberXY)))
      call FixProportionCoeff(KRED(BoxNumberXY),21,31,1.0D+00,r)
      !15
      r= exp(-alpha*(Dym- D2m(BoxNUmberXY)))
      call FixProportionCoeff(KRED(BoxNumberXY),41,51,1.0D+00,r)

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Calculate for the above defined set of boundary conditions
     ! the steady-state profiles and return the vertically integrated
     ! concentration.
     !
     ! Technical improvements: in case of utlimate low mineralization rates and
     ! a nearly empty reduction equivalent pool. There is a chance the
     ! estimated equilibrium value is negative. Therfore cK6r is limited to &
     ! values >=0
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     if ( InitializeModel== 0) then
         !16
         ! max. 0 order input available
         r=-p_sOS* p_qro*G2o(BoxNumberXY) &
                                       /D1m(BoxNumberXY)/p_poro(BoxNumberXY)
         call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 15, STANDARD, &
             dummy, value=r/p_sOS)

         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! Calculate the adaptation time to the steady-state profile
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          cK6r = max( 0.0D+00, CalculateSet( KRED(BoxNumberXY), &
                 SET_LAYER_INTEGRAL, LAYER1, LAYER1, 0.0D+00 , 0.0D+00))

          Tau  =   CalculateTau(  p_sOS,  diff,  p_p,  D1m(BoxNumberXY))

         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! Estimate the average value of K6r over the actual time step
         ! (transient value).
         ! This value depends on the adaptation time, the actual time step,
         ! the ''old'' value and the ''equilibrium value''
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

         cK6r = cK6r+( K6r(BoxNumberXY)- cK6r)* IntegralExp( -LocalDelta/Tau, 1.0D+00)

         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! Derive the equations for the transient profiles, assuming the same
         ! solution as for the steady-state case and using cK6r as new &
         ! constraint.
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          dummy = CalculateSet( KRED(BoxNumberXY), ADD, 0, 0, dummy, cK6r)

         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! Estimation of the Vertical flux at surface
         ! from the set of transient solutions:
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! loss flux due to denitrification
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           jK6G4r= p_qro* p_qon_dentri* jK3G4n(BoxNumberXY)
           call LimitChange(1,jK6G4r,K16r(BoxNumberXY),1.0D+00)
           call flux(BoxNumberXY, iiBen, ppK16r, ppK16r, -jK6G4r )
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! flux at the D1m
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           jK16K6r= CalculateFromSet( KRED(BoxNumberXY), DERIVATIVE, RFLUX, D1m(BoxNumberXY), dummy)

           ! Calculate mass shifted in upwards direction:
           Dnew  =   D1m(BoxNumberXY)+ LocalDelta* shiftD1m(BoxNumberXY)
           jK16K6r = jK16K6r+ CalculateFromSet( KRED(BoxNumberXY), SHIFT, LAYER1, &
                                   D1m(BoxNumberXY), Dnew)/ LocalDelta

           call LimitShift(jK16K6r,K6r(BoxNumberXY),K16r(boxNumberXY)-jK6G4r ,p_max_shift_change)

           call flux(BoxNumberXY, iiBen, ppK16r, ppK6r,   jK16K6r* insw(  jK16K6r) )
           call flux(BoxNumberXY, iiBen, ppK6r, ppK16r, - jK16K6r* insw( -jK16K6r) )


           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! flux at the 0.0
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           jbotN6r(BoxNumberXY) = CalculateFromSet( KRED(BoxNumberXY), DERIVATIVE, &
               RFLUX, 0.0D+00, dummy)

           call LimitChange(jbotN6r(BoxNumberXY),N6r_Ben(BoxNumberXY)*Depth_Ben(BoxNumberXY) ,&
                                                         p_max_state_change)

           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! Reoxidation Flux from (S2-, Fe3+, Mg3+ to SO4-, FE2+, Mg2+
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           jK6BTr = p_sOS* CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, &
             RFLUX, 0.0D+00, D1m(BoxNumberXY))


           r=jbotN6r(BoxNumberXY)+jK6BTr ; r2=r
           call LimitChange(1,r,K6r(BoxNumberXY)+jK16K6r,p_max_state_change)
           jbotN6r(BoxNumberXY) =jbotN6r(BoxNumberXY) *r/r2
           jK6BTr=jK6Btr * r/r2;
           call flux(BoxNumberXY, iiBen, ppK6r, ppK6r, -jK6BTr )

           call flux(BoxNumberXY, iiBen, ppK6r, ppK6r, -jbotN6r(BoxNumberXY) )

           jG2K7o(BoxNumberXY)  =   jK6BTr/ p_qro
           call flux(BoxNumberXY, iiBen, ppG2o, ppG2o, -jG2K7o(BoxNumberXY) )


           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! flux at the D2m
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           jK26K16r = -zuD1*exp(-alpha * ( D2m(BoxNumberXY)-D1m(BoxNumberXY))) &
                         * IntegralExp( -alpha, p_d_tot- D2m(BoxNumberXY))

           ! Calculate mass shifted in upwards direction:
           Dnew  =   D2m(BoxNumberXY)+ LocalDelta* shiftD2m(BoxNumberXY)
           jK26K16r = jK26K16r+ CalculateFromSet( KRED(BoxNumberXY), SHIFT, LAYER3, &
                                   D2m(BoxNumberXY), Dnew)/ LocalDelta

           jK26K16r = jK26K16r+ &
                   CalculateFromSet( KRED(BoxNumberXY), DERIVATIVE, RFLUX, D2m(BoxNUmberXY), dummy)
           call LimitShift(jK26K16r,K16r(BoxNumberXY)-jK6G4r,K26r(boxNumberXY),p_max_shift_change)

           call flux(BoxNumberXY, iiBen, ppK26r, ppK16r,   jK26K16r* insw(  jK26K16r) )
           call flux(BoxNumberXY, iiBen, ppK16r, ppK26r, - jK26K16r* insw( -jK26K16r) )

           jK36K26r(BoxNumberXY) = &
                   CalculateFromSet( KRED(BoxNumberXY), DERIVATIVE, RFLUX, p_d_tot_2, dummy)
           call flux(BoxNumberXY, iiBen, ppK26r, ppK26r,  jK36K26r(BoxNumberXY))
      else
          call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 21, PARAMETER, &
                dummy, value=zuD1)

          dummy = CalculateSet( KRED(BoxNumberXY), 0, 0,  0, dummy, dummy)

          jK6BTr = p_sOS* CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, &
                                         RFLUX, 0.0D+00, D1m(BoxNumberXY))
          jG2K7o(BoxNumberXY)  =   jK6BTr/ p_qro
      endif

  end do

  end subroutine BenAnoxicDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
