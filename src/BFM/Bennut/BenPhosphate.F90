#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenPhosphate
!
! DESCRIPTION
!   Description of the phosphate diagenitic processes in the sediment
!       Details on the equations and the method used to calculate
!       the equilibrium and transient profiles can be found in
!       Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
!
! !INTERFACE
  subroutine BenPhosphateDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: K1p, K11p, K21p, D1m, D2m, D8m, D2STATE_BEN
  use mem, ONLY: ppK1p, ppK11p, ppK21p, ppD1m, ppD2m, ppD8m, &
    NO_BOXES_XY,   &
    BoxNumberXY_ben, InitializeModel, LocalDelta, M1p, M11p, M21p, KPO4,KPO4_2, &
    jbotN1p, reBTp, reATp, irrenh, ETW_Ben, N1p_Ben, Depth_Ben, shiftD1m, shiftD2m, &
    jK31K21p, iiBen, iiPel, flux
#endif
  use constants, ONLY: QUADRATIC_TERM, ZERO_EXPONENTIAL_TERM, LAYERS, &
    LAYER1, LAYER2, LAYER3 ,DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, LAYER4, DEFINE, LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, &
    FLAG, MASS, SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, &
    SET_LAYER_INTEGRAL_UNTIL, INPUT_TERM,INPUT_ADD_TERM, PARAMETER, START_ADD_TERM, &
    STANDARD, ADD, DERIVATIVE, RFLUX, SHIFT, ONE_PER_DAY,INTEGRAL
  use mem_Param,  ONLY: p_poro, p_p_ae, p_d_tot,p_d_tot_2, p_clDxm, p_q10diff,p_small
  use mem_BenPhosphate
  use mem_BenthicNutrient3, ONLY:p_max_shift_change,p_InitCondition,p_max_state_change

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
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp, insw
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! !REVISION_HISTORY
!       September 1999 by M. Vichi !               Commented version
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij & M.Vichi
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
  integer  :: KSHIFT=0
  real(RLEN)  :: diff
  real(RLEN)  :: jK1N1p
  real(RLEN)  :: jK11K1p
  real(RLEN)  :: jK21K11p
  real(RLEN)  :: zuBT
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: alpha
  real(RLEN)  :: cK1p
  real(RLEN)  :: Tau
  real(RLEN)  :: clM1p
  real(RLEN)  :: s
  real(RLEN)  :: r
  real(RLEN)  :: dn
  real(RLEN)  :: dummy
  real(RLEN)  :: Dnew
  real(RLEN)  :: K11p_C
  real(RLEN)  :: K21p_C
  integer  :: i,j,mode
  real(RLEN)  :: ds1,ds2,layer_shift,p_shift,m_shift

  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D8.m is the average penetration depth for P-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      alpha  =   ONE/ max(  p_clDxm,  D8m(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Average in the oxic layer:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      zuBT = max( 1.079E-6_RLEN, reBTp(BoxNumberXY_ben))/ &
            p_poro(BoxNumberXY_ben)/ D1m(BoxNumberXY_ben)


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Anoxic Mineralization at D1.m, using the exponential distribution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if ( InitializeModel== 0) then
        zuD1 = max( p_small, reATp(BoxNumberXY_ben))/ p_poro(BoxNumberXY_ben)/ IntegralExp( &
          - alpha, p_d_tot- D1m(BoxNumberXY_ben))
      else
        zuD1 = max( p_small, reBTp(BoxNumberXY_ben))/ p_poro(BoxNumberXY_ben)/ &
            IntegralExp( - alpha, D1m(BoxNumberXY_ben))*exp(-alpha * D1m(BoxNumberXY_ben))
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Anoxic Mineralization at D2.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      dn=    0.75_RLEN* D2m(BoxNumberXY_ben) +0.25_RLEN*p_d_tot

      zuD2  =   zuD1* exp( - alpha*( D2m(BoxNumberXY_ben)- D1m(BoxNumberXY_ben)))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff = p_diff* irrenh(BoxNumberXY_ben)* p_poro(BoxNumberXY_ben)* &
        eTq( ETW_Ben(BoxNumberXY_ben), p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      mode=0;
      call  BenPhosphateEquation(KPO4(BoxNumberXY_ben),-InitializeModel, &
         N1p_Ben(BoxNumberXY_ben), K1p(BoxNumberXY_ben),K11p(BoxNumberXY_ben), &
         K21p(BoxNumberXY_ben),D1m(BoxNumberXY_ben),D2m(BoxNumberXY_ben), &
         p_d_tot,p_d_tot_2,p_poro(BoxNumberXY_ben),p_p_ae(BoxNumberXY_ben),p_p_an, p_s_ads, &
         p_shift,m_shift,layer_shift,dn, diff, alpha,zuBT,zuD1,zuD2)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration in the first layer.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      cK1p = CalculateSet( KPO4(BoxNumberXY_ben), SET_LAYER_INTEGRAL, LAYER1, &
        LAYER1, dummy, ZERO)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M1p(BoxNumberXY_ben) = K1p(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)/( &
        p_p_ae(BoxNumberXY_ben)+ ONE)/( D1m(BoxNumberXY_ben))
      M11p(BoxNumberXY_ben) =  CalculateFromSet(KPO4(BoxNumberXY_ben),INTEGRAL, &
              STANDARD,D1m(BoxNumberXY_ben),dn)/(dn-D1m(BoxNumberXY_ben));
      M21p(BoxNumberXY_ben) = K21p(BoxNumberXY_ben)/ p_poro(BoxNumberXY_ben)/( p_p_an+ &
        ONE)/( p_d_tot_2- dn )

      if ( InitializeModel== 0) then
         layer_shift=ShiftD2m(BoxNumberXY_ben)*LocalDelta;
         Dnew =   D2m(BoxNumberXY_ben)+ layer_shift

         p_shift =CalculateFromSet( KPO4(BoxNumberXY_ben), SHIFT, &
             LAYER3, D2m(BoxNumberXY_ben), Dnew)
         if ( abs(layer_shift) > 5.0e-7_RLEN) then
           mode=1;
           m_shift=zuD2*IntegralExp(-alpha,layer_shift) 
           call  BenPhosphateEquation(KPO4_2(BoxNumberXY_ben),1, &
            N1p_Ben(BoxNumberXY_ben),K1p(BoxNumberXY_ben),K11p(BoxNumberXY_ben), &
            K21p(BoxNumberXY_ben),D1m(BoxNumberXY_ben),D2m(BoxNumberXY_ben), &
            p_d_tot,p_d_tot_2,p_poro(BoxNumberXY_ben),p_p_ae(BoxNumberXY_ben),p_p_an, p_s_ads, &
            p_shift,m_shift,layer_shift,dn, diff, alpha,zuBT,zuD1,zuD2)
           KSHIFT=KPO4_2(BoxNumberXY_ben)
           cK1p = CalculateSet( KSHIFT, SET_LAYER_INTEGRAL, LAYER1, &
           LAYER1, dummy, ZERO)

           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! Calculate the adaptation time to the steady-state profile
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           Tau = CalculateTau( p_s_ads, diff, p_p_ae(BoxNumberXY_ben), &
             D1m(BoxNumberXY_ben))

           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! Estimate the average value of K1p over the actual time step
           ! (transient value).
           ! This value depends on the adaptation time, the actual time step,
           ! the ''old'' value and the ''equilibrium value''
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           cK1p = cK1p+( K1p(BoxNumberXY_ben)- cK1p)* IntegralExp( - LocalDelta/ &
             Tau, ONE)

           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! Derive the equations for the transient profiles, assuming the same
           ! solution as for the steady-state case and using cK1p as new &
           ! constraint.
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           dummy = CalculateSet( KSHIFT, ADD, 0, 0, dummy, cK1p)

         else

           layer_shift=ZERO
           KSHIFT=KPO4(BoxNumberXY_ben)

         endif
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Start calculation of fluxes:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Vertical fluxes :
        ! There are 2 problems with this model in this version both connected
        ! with the shifting of the layers:
        ! 1. shifting from the oxic+denitrification layer with high
        !  adsorbed fraction phosphate to the lower anoxic layer with a very
        !  low percentage of adsorped phosphate.
        ! 2. Too large changes in spring due to large change of D1.m:
        !  This leads sometimes to a calculated phosphate gradient which
        !  has negative values.
        !
        !  Solution:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate flux at the sediment/water interface:
        ! Check on; to high fluxes from pelagic and on concisteny of gradient
        ! ( only flux of M1p > N1p!)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jK1N1p =  CalculateFromSet( KSHIFT,DERIVATIVE, RFLUX, ZERO, ZERO)

        Call LimitShift(jK1N1p, &
             N1p_Ben(BoxNumberXY_ben)* Depth_Ben(BoxNumberXY_ben), &
             K1p(BoxNumberXY_ben), &
             p_max_state_change);
        call flux(BoxNumberXY_ben, iiBen, ppK1p, ppK1p, -jK1N1p )
        jbotN1p(BoxNumberXY_ben)=jbotN1p(BoxNumberXY_ben) + jK1N1p

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate diffusive flux at the oxic/denitrification interface:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jK11K1p = CalculateFromSet( KSHIFT, DERIVATIVE, RFLUX, &
          D1m(BoxNumberXY_ben), ONE)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate new depth of the oxygen horizon
        ! and the flux of phosphate related to this shifting
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        Dnew  =   D1m(BoxNumberXY_ben)+ shiftD1m(BoxNumberXY_ben)* LocalDelta
        jK11K1p =jK11K1p+  CalculateFromSet( KSHIFT, SHIFT, LAYER1, &
          D1m(BoxNumberXY_ben), Dnew)/ LocalDelta

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! limit flux according to the actual phosphate content in the layer
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call LimitShift(jK11K1p,K1p(BoxNumberXY_ben)-jK1N1p, &
             K11p(BoxNumberXY_ben), &
             p_max_shift_change)

        call flux(BoxNumberXY_ben, iiBen, ppK11p, ppK1p,   jK11K1p* insw( jK11K1p) )
        call flux(BoxNumberXY_ben, iiBen, ppK1p, ppK11p, - jK11K1p* insw(-jK11K1p) )
        

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate diffusive flux at the denitrification/anoxic interface:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jK21K11p=CalculateFromSet( KSHIFT, DERIVATIVE, RFLUX,  dn, ONE)
        Dnew  =   dn+ layer_shift
        jK21K11p =jK21K11p+  CalculateFromSet( KSHIFT, SHIFT, LAYER3+mode, &
          dn, Dnew)/ LocalDelta

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! All the nutrient mineralization source term in the anoxic layer
        ! has been added to K11.p in BenBacDynamics
        ! However in the model this layer is subdivided and hence a partition
        ! flux is here calculated according to the exponential distribution.
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jK21K11p = jK21K11p - zuD2*exp(-alpha*(dn-D2m(BoxNumberXY_ben))) &
           * p_poro(BoxNumberXY_ben)* IntegralExp( -alpha, p_d_tot- dn);

        call LimitShift(jK21K11p, &
             K11p(BoxNumberXY_ben)-jK11K1p, &
             K21p(BoxNumberXY_ben), &
             p_max_shift_change)
        call flux(BoxNumberXY_ben, iiBen, ppK21p, ppK11p, jK21K11p* insw( jK21K11p) )
        call flux(BoxNumberXY_ben, iiBen, ppK11p, ppK21p,-jK21K11p* insw(-jK21K11p) )


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate flux at the lower boundary
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jK31K21p(BoxNumberXY_ben) = CalculateFromSet( KSHIFT, DERIVATIVE, RFLUX, p_d_tot_2, ONE)

        call LimitChange(1, &
             jK31K21p(BoxNumberXY_ben), &
             K21p(BoxNumberXY_ben), &
             p_max_state_change)
        call flux(BoxNumberXY_ben, iiBen, ppK21p, ppK21p, jK31K21p(BoxNumberXY_ben) )

      end if

  end do

  end subroutine BenPhosphateDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenPhosphateEquation
!
! DESCRIPTION
  subroutine BenPhosphateEquation(KPO4,mode,N1p,K1p,K11p,K21p,D1m,D2m, &
              d_tot_1,d_tot_2,p_poro,p_p_ae,p_p_an,p_s_ads, &
              p_shift,m_shift,layer_shift,dn, diff, alpha,zuBT,zuD1,zuD2)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants, ONLY: QUADRATIC_TERM, ZERO_EXPONENTIAL_TERM, LAYERS, &
    LAYER1,LAYER2,LAYER3,LAYER5,LAYER6 ,DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, LAYER4, DEFINE, LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, &
    FLAG, MASS, SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, EXPONENTIAL_TERM, &
    SET_LAYER_INTEGRAL_UNTIL, INPUT_TERM,INPUT_ADD_TERM, PARAMETER, START_ADD_TERM, &
    STANDARD, ADD, DERIVATIVE, RFLUX, SHIFT, ONE_PER_DAY,INTEGRAL,LAYER7

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet, &
    CalculateSet, CalculateTau, CalculateFromSet

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp, insw
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij & M.Vichi
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
  IMPLICIT NONE
  integer,intent(OUT)     :: KPO4
  integer,intent(IN)     :: mode
  real(RLEN),intent(IN)  :: p_shift
  real(RLEN),intent(IN)  :: m_shift
  real(RLEN),intent(IN)  :: layer_shift
  real(RLEN),intent(IN)  :: dn
  real(RLEN),intent(IN)  :: alpha
  real(RLEN),intent(IN)  :: diff
  real(RLEN),intent(IN)  :: zuBT
  real(RLEN),intent(IN)  :: zuD1
  real(RLEN),intent(IN)  :: zuD2
  real(RLEN),intent(IN)  :: N1p
  real(RLEN),intent(IN)  :: K1p
  real(RLEN),intent(IN)  :: K11p
  real(RLEN),intent(IN)  :: K21p
  real(RLEN),intent(IN)  :: D1m
  real(RLEN),intent(IN)  :: D2m
  real(RLEN),intent(IN)  :: d_tot_1
  real(RLEN),intent(IN)  :: d_tot_2
  real(RLEN),intent(IN)  :: p_poro
  real(RLEN),intent(IN)  :: p_p_ae
  real(RLEN),intent(IN)  :: p_p_an 
  real(RLEN),intent(IN)  :: p_s_ads 
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: r
  real(RLEN)  :: pshift
  real(RLEN)  :: lambda
  real(RLEN)  :: dummy=0.0
  integer     :: i
  real(RLEN)  :: dz 
  real(RLEN)  :: ds1,ds2,p_ads_4

      if ( mode <= 0) then
         p_ads_4=p_p_an; ds1=D2m
         ds2=D2m
      else if (layer_shift.lt.ZERO) then
         p_ads_4=p_p_an
         ds1=D2m+min(-1.0e-4_RLEN,layer_shift)
         ds2=D2m
      else
         p_ads_4=p_p_ae
         ds1=D2m
         ds2=ds1 + max(1.0e-4_RLEN,layer_shift)
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Subdivide anoxic layer in two sublayers only for the calculation.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dz  =  ( d_tot_1+ dn )* 0.5_RLEN

      if ( mode .le.0) then
         KPO4 = InitializeSet( KPO4, 5, 14)
         call DefineSet( KPO4, LAYERS, LAYER1,  LAYER2,   D1m, D2m)
         call DefineSet( KPO4, LAYERS, LAYER3,  LAYER4,   dn, dz)
      else
         KPO4 = InitializeSet( KPO4, 6, 18)
         call DefineSet( KPO4, LAYERS,  LAYER1,  LAYER2, D1m, ds1 )
         call DefineSet( KPO4, LAYERS,  LAYER3,  LAYER4, ds2, dn)
         call DefineSet( KPO4, LAYERS,  LAYER5,  0     ,  dz,  dummy)
      endif


      call DefineSet( KPO4, DIFFUSION, FOR_ALL_LAYERS, 0, diff, dummy)

      call DefineSet( KPO4, POROSITY, FOR_ALL_LAYERS, 0, p_poro, dummy)

      call DefineSet( KPO4, ADSORPTION, LAYER1, LAYER2, p_p_ae, p_p_ae)
      call DefineSet( KPO4, ADSORPTION, LAYER3, 0, p_ads_4,dummy)
      if ( mode.le.0) then
         call DefineSet( KPO4, ADSORPTION, LAYER4, LAYER5, p_p_an, p_p_an)
      else
         call DefineSet( KPO4, ADSORPTION, LAYER4, LAYER5, p_p_an, p_p_an)
         call DefineSet( KPO4, ADSORPTION, LAYER6, 0, p_p_an, dummy)
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! P(z) = p13*z^2 + p14*z + p15
      ! 2nd layer:
      ! P(z) = p21*exp[-alpha*(z-D1.m)] + p24*z + p25
      ! 3rd layer:
      ! P(z) = p31*exp[-alpha*(z-D2.m)] + p34*z + p35
      ! 4th layer:
      ! P(z) = p41*exp[-alpha*(z-dn)] + p44*z + p45
      !    p44 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      lambda=sqrt(p_s_ads/diff); 
      call DefineSet( KPO4, DEFINE, 11, EXPONENTIAL_TERM, -lambda, dummy)
      call DefineSet( KPO4, DEFINE, 12, EXPONENTIAL_TERM,  lambda, dummy)
      call DefineSet( KPO4, DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KPO4, DEFINE, 21, ZERO_EXPONENTIAL_TERM, -alpha, dummy)
      call DefineSet( KPO4, DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet( KPO4, DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      if ( mode.le.0) then
        i=30
      else
        call DefineSet( KPO4, DEFINE, 31, EXPONENTIAL_TERM, -lambda,  dummy)
        call DefineSet( KPO4, DEFINE, 32, EXPONENTIAL_TERM,  lambda, dummy)
        call DefineSet( KPO4, DEFINE, 34, LINEAR_TERM, dummy, dummy)
        call DefineSet( KPO4, DEFINE, 35, CONSTANT_TERM, dummy, dummy)
        i=40
      endif

      call DefineSet( KPO4, DEFINE, i+1, ZERO_EXPONENTIAL_TERM, -alpha, dummy)
      call DefineSet( KPO4, DEFINE, i+4, LINEAR_TERM, dummy, dummy)
      call DefineSet( KPO4, DEFINE, i+5, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KPO4, DEFINE, i+11, ZERO_EXPONENTIAL_TERM, -alpha, dummy)
      call DefineSet( KPO4, DEFINE, i+14, LINEAR_TERM, dummy, dummy)
      call DefineSet( KPO4, DEFINE, i+15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KPO4, DEFINE, i+21, ZERO_EXPONENTIAL_TERM, - alpha, dummy)
      call DefineSet( KPO4, DEFINE, i+25, CONSTANT_TERM, dummy, dummy) 

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !1-8/1-10
      call CompleteSet( KPO4, SET_CONTINUITY, FLAG, MASS, dummy)

      !9/11
      call CompleteSet( KPO4, SET_BOUNDARY, LAYER1, EQUATION, ZERO, value=N1p)


      select case ( mode )
        case ( 0 )
          !10-11
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL, LAYER2,LAYER3, dummy, value=K11p)
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL_UNTIL, LAYER4, LAYER5, d_tot_2, value=K21p)

         !12: condition for 30/40 layer....
          call CompleteSet(KPO4,INPUT_TERM,21,PARAMETER,dummy,value=zuD2)
        case ( 1 )
          !12
          r=+p_s_ads * abs(p_shift/layer_shift)/p_ads_4/p_poro + abs(m_shift/layer_shift);
          call CompleteSet( KPO4, INPUT_TERM, 35, STANDARD, dummy, value=r/p_s_ads)

          !13:  next two lines one boundary condition!
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL,  LAYER2, LAYER2,  dummy, value=ZERO)
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL, -LAYER4, LAYER4, dummy, value=K11p-p_shift)

          !14
          call CompleteSet(KPO4,INPUT_TERM,21,PARAMETER,dummy,value=zuD2)
          !15
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL_UNTIL, LAYER5, LAYER6, d_tot_2, value=K21p)

          r  =   exp( - alpha*( dn- ds2))
          call FixProportionCoeff(KPO4,41,51,ONE,r)
        case ( -1 )
          !10-11
          ! The mineralization at D1m equal to the oxic minerlaization at D1m under
          ! assumuption of that the mineralization distribution in oxic layer is distributed
          ! according detritus distribution alpha 
          ! 12-13
          call CompleteSet( KPO4, INPUT_TERM, 21, PARAMETER, dummy, value=zuD1)
          call CompleteSet( KPO4, INPUT_TERM, 41, PARAMETER, dummy,&
                                               value=zuD1* exp(-alpha*(dn-D2m)))
          r  =   exp( - alpha*( D2m- D1m))
          call FixProportionCoeff(KPO4,21,31,ONE,r)
      end select

      !13/!17 : condition for fifth/sixth layer....
      r  =   exp( - alpha*( dz- dn))
      call FixProportionCoeff(KPO4,i+11,i+21,ONE,r)

      !14/18
      r=+p_s_ads * K1p/p_p_ae/p_poro/D1m+zuBT
      call CompleteSet( KPO4, INPUT_TERM, 15, STANDARD, dummy, value=r/p_s_ads)
  end subroutine BenPhosphateEquation
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
