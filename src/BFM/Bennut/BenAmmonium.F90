#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAmmonium
!
! DESCRIPTION
!   Description of the diagenetic ammonium processes in the sediment
!   Details on the equations and the method used to calculate
!   the equilibrium and transient profiles can be found in
!   Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenAmmoniumDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: K4n, K3n, G2o, D1m, K14n, D2m, K24n, D7m, D2STATE
  use mem, ONLY: ppK4n, ppK3n, ppG2o, ppD1m, ppK14n, ppD2m, &
    ppK24n, ppD7m, NO_BOXES_XY, &
       BoxNumberXY, InitializeModel, LocalDelta, &
    M4n, KNH4, jG2K3o, jbotN4n, M14n, M24n, reBTn, reATn, irrenh, ETW_Ben, &
    jK4K3n, N4n_Ben, Depth_Ben, iiBen, iiPel, flux
#endif
  use constants, ONLY: LAYERS, LAYER1, LAYER2, &
    DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DOUBLE_DEFINE, &
    EXPONENTIAL_TERM, DEFINE, CONSTANT_TERM, ZERO_EXPONENTIAL_TERM, LINEAR_TERM, &
    SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, &
    SET_LAYER_INTEGRAL_UNTIL, LAYER3, INPUT_TERM, PARAMETER, STANDARD, ADD, &
    INTEGRAL, RFLUX, DERIVATIVE,INPUT_ADD_TERM,START_ADD_TERM
  use mem_Param,  ONLY: p_poro, p_d_tot, p_d_tot_2,p_clDxm, p_q10diff, p_qon_nitri
  use mem_BenthicNutrient3, ONLY:p_max_state_change, p_InitCondition

  use mem_BenAmmonium


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet, &
    CalculateSet, CalculateTau, CalculateFromSet


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: insw,IntegralExp

!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! !REVISION_HISTORY
!   September 1999 by M. Vichi !               Commented version
!
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, the BFM team 
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
  real(RLEN)  :: zuBT
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: alpha
  real(RLEN)  :: diff
  real(RLEN)  :: Tau
  real(RLEN)  :: labda
  real(RLEN)  :: a15
  real(RLEN)  :: cK4n
  real(RLEN)  :: jK4N4n
  real(RLEN)  :: sK4K3
  real(RLEN)  :: cO2
  real(RLEN)  :: r
  real(RLEN)  :: dummy

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      ! Calculate pore-water oxygen concentration in the oxic layer
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M4n(BoxNumberXY) = K4n(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ &
        ONE)/( D1m(BoxNumberXY))
      M14n(BoxNumberXY) = K14n(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ &
        ONE)/( D2m(BoxNumberXY)- D1m(BoxNumberXY))
      M24n(BoxNumberXY) = K24n(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p+ &
        ONE)/( p_d_tot_2- D2m(BoxNumberXY))

      cO2  =   max(1.D-20,G2o(BoxNumberXY)/ p_poro(BoxNumberXY)/ D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D7.m is the average penetration depth for N-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      alpha  =   ONE/ max(  p_clDxm,  D7m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Average in the oxic layer:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        zuBT  =   reBTn(BoxNumberXY)/ p_poro(BoxNumberXY)/ D1m(BoxNumberXY)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Anoxic Mineralization at D1.m, using the exponential distribution
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        zuD1 = max( 1.D-20, reATn(BoxNumberXY))/ p_poro(BoxNumberXY)/ IntegralExp( &
          - alpha, p_d_tot_2- D1m(BoxNumberXY))


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! nitrification rate: temperature and oxygen
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff = p_diff* p_poro(BoxNumberXY)* irrenh(BoxNumberXY)* &
      eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      sK4K3 = p_sK4K3* eTq( ETW_Ben(BoxNumberXY), p_q10) 

      if (InitializeModel == 0 ) &
        sK4K3=sK4K3* cO2/( cO2+ p_clO2)* M4n(BoxNumberXY)/( M4n(BoxNumberXY)+ p_clM4)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! if availability of carbon for degradation is low , the nitrification &
      ! will be hampered
      ! by the lack of carbon for nitrification bacteria. As proxy for &
      ! the degrdability
      ! of the carbon the respiration per mass of bacteria is caclulated.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      labda  =   sqrt(  sK4K3/ diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient of the zero order term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      a15  =   zuBT/ sK4K3

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, p_porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      KNH4(BoxNumberXY) = InitializeSet( KNH4(BoxNumberXY), 3, 8)

      call DefineSet( KNH4(BoxNumberXY), LAYERS, LAYER1, &
        LAYER2, D1m(BoxNumberXY), D2m(BoxNumberXY))

      call DefineSet( KNH4(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
                                                                        dummy)

      call DefineSet( KNH4(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet( KNH4(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, &
        dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! A(z) = a11*exp(labda*z) + a12*exp(-labda*z) + a13*z^2 + a14*z + a15
      ! 2nd layer:
      ! A(z) = a21*exp[-alpha*(z-D1.m)] + a24*z + a25
      ! 3rd layer:
      ! A(z) = a31*exp[-alpha*(z-D2.m)] + a34*z + a35
      !    a34 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call DefineSet( KNH4(BoxNumberXY), DOUBLE_DEFINE, 11, EXPONENTIAL_TERM, &
        labda, sK4K3)
      call DefineSet( KNH4(BoxNumberXY), DOUBLE_DEFINE, 12, EXPONENTIAL_TERM, - &
        labda, sK4K3)
      call DefineSet( KNH4(BoxNumberXY), DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 21, ZERO_EXPONENTIAL_TERM, -alpha, dummy)
      call DefineSet( KNH4(BoxNumberXY), DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet( KNH4(BoxNumberXY), DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KNH4(BoxNumberXY), DEFINE, 31, ZERO_EXPONENTIAL_TERM, -alpha, dummy)
      call DefineSet( KNH4(BoxNumberXY), DEFINE, 35, CONSTANT_TERM, dummy, dummy)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !4
      call CompleteSet( KNH4(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, dummy)
      !5
      call CompleteSet( KNH4(BoxNumberXY), SET_BOUNDARY, LAYER1, &
        EQUATION, ZERO, value=N4n_Ben(BoxNumberXY))

      select case ( InitializeModel)
        case ( 0 )
          call CompleteSet( KNH4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER2, &
            LAYER2, dummy, value=K14n(BoxNumberXY))

          call CompleteSet( KNH4(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
            LAYER3, LAYER3, p_d_tot_2, value=K24n(BoxNumberXY))
 
          call CompleteSet( KNH4(BoxNumberXY), INPUT_TERM, 15, STANDARD, dummy, &
            value=a15)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Calculate for the above defined set of boundary conditions
          ! the steady-state profiles and return the vertically integrated
          ! concentration in the first layer.
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          cK4n = CalculateSet( KNH4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
            LAYER1, dummy, ZERO)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Calculate the adaptation time to the steady-state profile
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          Tau  =   CalculateTau(  sK4K3,  diff,  p_p,  D1m(BoxNumberXY))

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Estimate the average value of K4n over the actual time step
          ! (transient value).
          ! This value depends on the adaptation time, the actual time step,
          ! the ''old'' value and the ''equilibrium value''
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          cK4n = cK4n+( K4n(BoxNumberXY)- cK4n)* IntegralExp( - LocalDelta/ &
            Tau, ONE)
  
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Derive the equations for the transient profiles, assuming the same
          ! solution as for the steady-state case and using cK4n as new &
          ! constraint.
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          dummy = CalculateSet( KNH4(BoxNumberXY), ADD, 0, 0, dummy, cK4n)
  
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Start calculation of fluxes:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Next calculation is done in BenNitrogenShiftin which
          ! include all fluxes between layers!
          ! All the nutrient mineralization source term in the anoxic layer
          ! has been added to K14.n in BenBacDynamics
          ! However in the model this layer is subdivided and hence a partition
          ! flux is here calculated according to the exponential distribution.
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          !   Calculate Anoxic Mineralization at D2.m
          !
          !  zuD2 = zuD1 * exp( -alpha * (D2.m - D1.m));
          !
          ! K14.n -> K24.n = zuD2 * p_poro * IntegralExp(-alpha, p_d_tot - D2.m);
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Calculate Nitrification flux in the first layer and the related
          ! oxygen consumption flux:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          jK4K3n(BoxNumberXY)= sK4K3* CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, &
            RFLUX, ZERO, D1m(BoxNumberXY))

          call flux(BoxNumberXY, iiBen, ppK4n, ppK3n, jK4K3n(BoxNumberXY) )

          jG2K3o(BoxNumberXY)  =   jK4K3n(BoxNumberXY)* p_qon_nitri
          call flux(BoxNumberXY, iiBen, ppG2o, ppG2o, -( jG2K3o(BoxNumberXY)) )

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Estimation of the Vertical fluxes from the set of transient solutions:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          jK4N4n = CalculateFromSet( KNH4(BoxNumberXY), &
                                         DERIVATIVE, RFLUX, ZERO, ZERO)
          !avoid too large flxues leading to negative concentrations
          call LimitShift(jK4N4n,N4n_Ben(BoxNumberXY)*Depth_Ben(BoxNumberXY) ,&
                          K4n(boxNumberXY),p_max_state_change)

          call flux(BoxNumberXY, iiBen, ppK4n, ppK4n, -( jK4N4n) )

          jbotN4n(BoxNumberXY)=jbotN4n(BoxNumberXY)+jK4N4n

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          !  All transport between layers are done in BenNitrogenShifting:
          !
          !  jK14K4n = CalculateFromSet(KNH4, DERIVATIVE, RFLUX, D1.m, 0.0);
          !  jK24K14n = CalculateFromSet(KNH4, DERIVATIVE, RFLUX, D2.m, 0.0);
          !  jK34K24n = CalculateFromSet(KNH4, DERIVATIVE, RFLUX, p_d_tot, 0.0);
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        case ( 1 )    ! Determination of Initial condition:
          ! The mineralization at D1m equal to the oxic mineralization at D1m under
          ! assumption that the mineralization distribution in oxic layer is distributed
          ! according detritus distribution alpha 
          zuD1= reBTn(BoxNumberXY) / p_poro(BoxNumberXY)/   &
            IntegralExp( - alpha, D1m(BoxNumberXY)) *exp(-alpha * D1m(BoxNumberXY))

          !6 
          r  =   exp( - alpha*( D2m(BoxNUmberXY)- D1m(BoxNumberXY)))
          call FixProportionCoeff(KNH4(BoxNumberXY),21,31,ONE,r)

          !7
          call CompleteSet( KNH4(BoxNumberXY), INPUT_TERM, 15, STANDARD, dummy, &
            value=a15)

          select case ( p_InitCondition )
          case(2)
            !8
            call CompleteSet( KNH4(BoxNumberXY), INPUT_TERM, 21, PARAMETER, &
            dummy, value=zuD1)
            !9
!           call CompleteSet( KNH4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
!           LAYER1, dummy, value=16.0 * N4n_Ben(BoxNumberXY)* p_poro(boxNumberXY) * 0.02 )
          end select

          dummy = CalculateSet( KNH4(BoxNumberXY), 0, 0,  0, dummy, dummy)

          jK4K3n(BoxNumberXY) = sK4K3* CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, &
            RFLUX, ZERO, D1m(BoxNumberXY))
          jG2K3o(BoxNumberXY)  =   jK4K3n(BoxNumberXY)* p_qon_nitri

          if ( jK4K3n(BoxNumberXY) < ZERO) then
              write(LOGUNIT,'(''Negative value calculated for K4n'')') 
              write(LOGUNIT,'(''Minrealisation (reBTT) ='',F10.3)') reBTn(BoxNumberXY)
              write(LOGUNIT,'(''Temperature (ETW) ='',F10.3)') ETW_Ben(BoxNumberXY)
              write(LOGUNIT,'(''Amommonium at s/w/interface (N4n_Ben) ='',F10.3)') N4N_Ben(BoxNumberXY)
              call PrintSet(KNH4(boxNumberXY),"Error in Initialization of ammonium gradient")
          endif

      end select

  end do

  end subroutine BenAmmoniumDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
