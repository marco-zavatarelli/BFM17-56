#include "cppdefs.h"
#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !MODULE: CO2System
!
! DESCRIPTION
!       Calculate CO2 equilibrium in seawater starting from 
!       total alkalinity and total CO2 (DIC) at
!       a given temperature, salinity and other species
!       Routines are wrapped in a module in order to be used for pelagic
!       and benthic computations
!
!       INPUTS (BFM variables)
!       DIC_in (O3c) = total dissolved inorganic carbon (umol C/kg) 
!       alk    (O3h) = alkalinity (umol eq/kg)
!       N1p     = inorganic phosphate (mmol/m3) 
!       N5s     = inorganic silicate (mmol/m3) 
!       ETW     = temperature (degrees C)
!       ESW     = salinity (PSU)
!       ERHO    = density (kg/m3)
!       EPCO2air        = atmospheric CO2 partial pressure (ppm) 
!
!       diagnostic OUTPUTS
!       co2  = CO2*(aq) (umol/kg)
!       pco2 = oceanic pCO2 (uatm)
!       hco3 = bicarbonate (umol/kg)
!       co3  = carbonate (umol/kg)
!       co2airflux = air-sea CO2 flux (mmol/m^2/d)
!
!       The code has 3 options defined in mod_pelCO2:
!       choice of the equilibrium constants K1K2
!       choice of the pH scale PHSCALE
!       choice of the computation modeCO2
!
! ORIGINAL REFERENCES
!--------------------------------------------------------------------------
!   J. Orr (LODYC) and the OCMIP team 
!   co2calc.f
!   Revision: 1.8  Date: 1999/07/16 11:40:33 
!---------------------------------------------------------------------
!   Zeebe & Wolf-Gladrow, CO2 in Seawater: 
!               Equilibrium, Kinetics, Isotopes. 2001. Elsevier
!               Matlab scripts: csys.m; equic.m
!---------------------------------------------------------------------
!
! AUTHORS
!   M. Vichi and L. Patara (INGV-CMCC, vichi@bo.ingv.it)
!   P. Ruardij and H. Thomas (NIOZ, rua@nioz.nl)
! 
! CHANGE_LOG
!
! COPYING
!   Copyright (C) 2007 M. Vichi and P. Ruardij (vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details. 
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !INTERFACE
  module CO2System

! !USES:
   use global_mem, ONLY:RLEN
! Shared variables
  implicit none
  real(RLEN),public  ::  K0 ! solubility : [Co2]=k0Ac PCO2
  real(RLEN),public  ::  K1 ! carbonate equilibrium I
  real(RLEN),public  ::  K2 ! carbonate equilibrium II
  real(RLEN),public  ::  Kw ! water dissociation
  real(RLEN),public  ::  Kb ! constant for Boron equilibrium
  real(RLEN),public  ::  Ks ! constant for bisulphate equilibrium
  real(RLEN),public  ::  Kf ! constant for hidrogen fluoride equilibirum
  real(RLEN),public  ::  Kp(3) ! constants for phosphate equilibirum
  real(RLEN),public  ::  Ksi ! constant for silicic acid eqilbrium


  real(RLEN),public  :: ldic  ! local DIC in mol/kg 
  real(RLEN),public  :: lpCO2 ! local pCO2
  real(RLEN),public  :: scl   ! chlorinity
  real(RLEN),public  :: ta    ! total alkalinity
  integer,public     :: way

  real(RLEN)         :: bt,ft,st,pt,sit
  real(RLEN),parameter   :: T1=1.0_RLEN,T2=2.0_RLEN,T3=3.0_RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  functions 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  contains

!-------------------------------------------------------------------------!
!BOP
! !IROUTINE:  CalcCO2System
!
! !INTERFACE
  function CalcCO2System(mode,salt,temp,rho,n1p,n5s,alk,CO2,HCO3,CO3,pH, &
                         DIC_in,pCO2_in,DIC_out,pCO2_out)
! !DESCRIPTION
! See module preamble
!
! USES:
  use global_mem, ONLY: ONE,ZERO
  use constants, ONLY: ZERO_KELVIN,MW_C
  use mem_CO2

! !INPUT PARAMETERS:
  IMPLICIT NONE
  integer                          :: mode
  real(RLEN),intent(IN)            :: salt
  real(RLEN),intent(IN)            :: temp
  real(RLEN),intent(IN)            :: rho 
  real(RLEN),intent(IN)            :: n1p
  real(RLEN),intent(IN)            :: n5s
  real(RLEN),intent(IN)            :: alk
  real(RLEN),intent(IN),optional   :: DIC_in
  real(RLEN),intent(IN),optional   :: pCO2_in
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
  real(RLEN),intent(OUT)           :: CO2
  real(RLEN),intent(OUT)           :: HCO3
  real(RLEN),intent(OUT)           :: CO3
  real(RLEN),intent(INOUT)         :: pH
  real(RLEN),intent(OUT),optional  :: DIC_out
  real(RLEN),intent(OUT),optional  :: pCO2_out
  integer                          :: CalcCO2System

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),parameter             :: MEG=1.D6,XACC=1.D-20,  &
                                      PERMIL=ONE/1000.0_RLEN,&
                                      PERMEG=ONE/MEG

  logical     :: small_interval
  integer     :: i, l, error
  real(RLEN)  :: Hplus, Hplus2
  real(RLEN)  :: tmp1, tmp2
  real(RLEN)  :: intercept,lnk,sit,pt
  real(RLEN)  :: h1,h2
  real(RLEN)  :: tk,tk100,tk1002,    &
                 invtk,is,is2,       &
                 dlogtk,dsqrtis,s,s2,&
                 dsqrts,s15
  ! Local variables fro Follows' parameterization
  real(RLEN) :: xx,xx2,xx3,k11,k12,k12p,k123p,k3p,a,c
  real(RLEN) :: cag,dummy,gamm
!EOP
!-------------------------------------------------------------------------!
!BOC

  !---------------------------------------------------------------------
  ! Change units from the input of mmol/m^3 -> mol/kg:
  ! (1 mmol/m^3)  x (1 m^3/1024.5 kg) / 1000
  ! The ocean actual density ERHO is used.
  ! Note: mol/kg are actually what the body of this routine uses 
  ! for all calculations.  
  !---------------------------------------------------------------------
  pt = n1p/rho*PERMIL
  sit = n5s/rho*PERMIL
  if (present(DIC_in)) then
     ! convert from umol/kg to standard DIC units mol/kg 
     ldic = DIC_in*PERMEG
     way = 1 
  elseif (present(pCO2_in)) then
     lpCO2  = pCO2_in
     way = 2 
  endif
  ! convert input variable alkalinity from umol/kg to mol/kg 
  ta = alk*PERMEG

  ! ---------------------------------------------------------------------
  ! Calculate all constants needed to convert between various measured
  ! carbon species. References for each equation are noted in the code. 
  ! The original version of this code was based on the code by Dickson 
  ! in Version 2 of "Handbook of Methods
  ! for the Analysis of the Various Parameters of the Carbon Dioxide System
  ! in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
  !
  ! Derive simple terms used more than once
  ! ---------------------------------------------------------------------
  tk = temp - ZERO_KELVIN
#ifdef DEBUG
  LEVEL2 'Entering CO2System...'
  LEVEL3 'tk',tk,'way',way,'mode',mode
  LEVEL3 'dic',ldic,'ta',ta
  LEVEL3 'pt',pt,'sit',sit
#endif
  tk100 = tk/100.0_RLEN
  tk1002 = tk100*tk100
  invtk = ONE/tk
  dlogtk = dlog(tk)
  ! salinity
  s  = salt
  s2 = salt*salt
  dsqrts = dsqrt(salt)
  s15 = salt**1.5_RLEN
  ! chlorinity
  scl = salt/1.80655_RLEN  
  ! ionic strength
  is = 19.924_RLEN*salt/ (1000._RLEN-1.005_RLEN*salt)
  is2 = is*is
  dsqrtis = dsqrt(is)

  !------------------------------------------------------------------------
  ! Calculate concentrations for borate, sulfate, and fluoride
  ! as a function of chlorinity
  ! ---------------------------------------------------------------------
  ! Uppstrom (1974)
  bt = 0.000232_RLEN * scl/10.811_RLEN
  ! Morris & Riley (1966)
  st = 0.14_RLEN * scl/96.062_RLEN
  ! Riley (1965)
  ft = 0.000067_RLEN * scl/18.9984_RLEN

!MAV this part is never used
  ! ---------------------------------------------------------------------
  ! convert partial pressure to fugacity
  ! ff = k0(1-pH2O)*correction term for non-ideality
  !
  ! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
  ! ---------------------------------------------------------------------
  !ff = exp(-162.8301_RLEN + 218.2968_RLEN/tk100  +      &
  !     90.9241_RLEN*dlog(tk100) - 1.47696_RLEN*tk1002 +  &
  !     s * (.025695_RLEN - .025225_RLEN*tk100 +         &
  !     0.0049867_RLEN*tk1002))

  ! ---------------------------------------------------------------------
  ! K0, solubility of CO2 in the water (K Henry)
  ! from Weiss 1974; K0 = [CO2]/pCO2 [mol kg-1 atm-1]
  ! ---------------------------------------------------------------------
  K0 = CalcK0(salt,temp)

  ! ---------------------------------------------------------------------
  ! Choice of Acidity constants from parameter in CO2 namelist
  ! k1 = [H][HCO3]/[H2CO3]
  ! k2 = [H][CO3]/[HCO3]
  ! If Using Follows et al. parameterization, force K1K2=1
  ! ---------------------------------------------------------------------
  select case (K1K2)
  case (1)
     ! ---------------------------------------------------------------------
     ! Constants according to Roy et al. (1993a). 
     ! Recommended by DOE (1994) and Zeebe / Wolf-Gladrow. 
     ! pH scale: total 
     ! ---------------------------------------------------------------------
     lnK = -2307.1266_RLEN*invtk +2.83655_RLEN-1.5529413_RLEN*dlogtk +  &
          (-4.0484_RLEN*invtk - 0.20760841_RLEN) * dsqrts +            &
          0.08468345_RLEN* s -0.00654208_RLEN * s15+dlog(ONE -0.001005_RLEN* s )
     K1 = dexp(lnK)
     lnK = -3351.6106_RLEN/tk -9.226508_RLEN-0.2005743_RLEN*dlogtk+        &
          (-23.9722_RLEN/tk - 0.10690177_RLEN)*dsqrts + 0.1130822_RLEN* s - &
          0.00846934_RLEN*s15 + dlog(ONE-0.001005_RLEN* s )
     K2 = dexp(lnK)
  case (2)
     ! ---------------------------------------------------------------------
     ! Mehrbach et al. (1973) as refitted by Dickson and Millero (1987) 
     ! pH scale: seawater  (Millero, 1995, p.664)
     ! Standard OCMIP computation. Natural seawater
     ! ---------------------------------------------------------------------
     K1 = 10.0_RLEN**(-ONE*(3670.7_RLEN*invtk - 62.008_RLEN + 9.7944_RLEN*dlogtk - &
          0.0118_RLEN * s + 0.000116_RLEN*s2))
     K2 = 10.0_RLEN**(-1.0*(1394.7_RLEN*invtk + 4.777_RLEN - &
          0.0184_RLEN*s + 0.000118_RLEN*s2))
  case (3)
     ! ---------------------------------------------------------------------
     ! Mehrbach et al (1973) refit by Lueker et al. (2000). 
     ! pH scale: total 
     ! ---------------------------------------------------------------------
     K1 = 10.0_RLEN**(-ONE*(3633.86_RLEN*invtk - 61.2172_RLEN + 9.6777_RLEN*dlogtk - &
          0.011555_RLEN * s + 0.0001152_RLEN * s2))
     K2 = 10.0_RLEN**(-ONE*(471.78_RLEN*invtk + 25.9290_RLEN - &
          3.16967_RLEN*dlogtk - 0.01781_RLEN * s + 0.0001122_RLEN * s2))
  case (4)
     !-----------------------------------------------------------------------
     ! Hansson (1973b) data as refitted by Dickson and Millero (1987).
     ! pH scale: seawater 
     !-----------------------------------------------------------------------
     K1 = 10.0_RLEN**(-ONE*(851.4_RLEN*invtk + 3.237_RLEN - &
          0.0106_RLEN*s + 0.000132_RLEN*s2))
     K2 = 10.0_RLEN**(-ONE*(-3885.4_RLEN*invtk + 125.844_RLEN - 18.141_RLEN*dlogtk -  &
          0.0192_RLEN*s + 0.000132_RLEN* s2))

  end select

  ! ---------------------------------------------------------------------
  ! k1p = [H][H2PO4]/[H3PO4] 
  ! pH scale: total
  !
  ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
  ! Millero p.670 (1995)
  ! ---------------------------------------------------------------------
  lnK = -4576.752_RLEN*invtk + 115.525_RLEN - 18.453_RLEN * dlogtk + &
       (-106.736_RLEN*invtk + 0.69171_RLEN) * dsqrts +               &
       (-0.65643_RLEN*invtk - 0.01844_RLEN) * s
  Kp(1) = dexp(lnK)
  ! ---------------------------------------------------------------------
  ! k2p = [H][HPO4]/[H2PO4] 
  ! pH scale: total
  !
  ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
  ! Millero p.670 (1995)
  ! ---------------------------------------------------------------------
  lnK = -8814.715_RLEN*invtk + 172.0883_RLEN - 27.927_RLEN * dlogtk + &
       (-160.340_RLEN*invtk + 1.3566_RLEN) * dsqrts +                 &
       (0.37335_RLEN*invtk - 0.05778_RLEN) * s
  Kp(2) = dexp(lnK)
  !------------------------------------------------------------------------
  ! k3p = [H][PO4]/[HPO4] 
  ! pH scale: total
  !
  ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  ! ---------------------------------------------------------------------
  lnK = -3070.75_RLEN*invtk - 18.126_RLEN + &
       (17.27039_RLEN*invtk + 2.81197_RLEN) *   &
       dsqrts + (-44.99486_RLEN*invtk - 0.09984_RLEN) * s
  Kp(3) = dexp(lnK)

  !------------------------------------------------------------------------
  ! ksi = [H][SiO(OH)3]/[Si(OH)4] pH on Sea Water Scale
  ! pH scale: total
  !
  ! DOE(1994) as recommended by Zeebe and Wolf-Gladrow
  ! ---------------------------------------------------------------------
  ksi = -8904.2_RLEN*invtk + 117.385_RLEN - 19.334_RLEN * dlogtk + &
       (-458.79_RLEN*invtk + 3.5913_RLEN) * dsqrtis +              &
       (188.74_RLEN*invtk - 1.5998_RLEN) * is +                   &
       (-12.1652_RLEN*invtk + 0.07871_RLEN) * is2 +               &
       dlog(ONE-0.001005_RLEN*s)
  Ksi = dexp(lnK)

  !------------------------------------------------------------------------
  ! kw = [H][OH]  ion product of water
  ! pH scale: see below
  ! ---------------------------------------------------------------------
  select case (phscale)
  case (SWS)
     ! Millero p.670 (1995) using composite data recommended DOE (1994)
     intercept = 148.9802_RLEN
  case (TOTAL)
     ! DOE (1994)
     intercept = 148.96502_RLEN
  end select
  lnK = intercept -13847.26_RLEN*invtk - 23.6521_RLEN * dlogtk + &
       (118.67_RLEN*invtk - 5.977_RLEN + 1.0495_RLEN * dlogtk) *          &
       dsqrts - 0.01615_RLEN * s
  Kw = dexp(lnK)

  !------------------------------------------------------------------------
  ! ks = [H][SO4]/[HSO4] 
  ! pH scale: "free"
  !
  ! Dickson (1990, J. chem. Thermodynamics 22, 113)
  ! ---------------------------------------------------------------------
  lnK = -4276.1_RLEN*invtk + 141.328_RLEN - 23.093_RLEN*dlogtk +      &
       (-13856._RLEN*invtk + 324.57_RLEN - 47.986_RLEN*dlogtk) * dsqrtis + &
       (35474._RLEN*invtk - 771.54_RLEN + 114.723_RLEN*dlogtk) * is -     &
       2698._RLEN*invtk*is**1.5_RLEN + 1776._RLEN*invtk*is2 +              &
       dlog(ONE - 0.001005_RLEN*s)
  Ks = dexp(lnK)

  !------------------------------------------------------------------------
  ! kf = [H][F]/[HF] pH on free scale
  ! pH scale: "free"
  !
  ! Dickson and Riley (1979)  also Dickson and Goyet (1994)
  ! ---------------------------------------------------------------------
  lnK = 1590.2_RLEN*invtk - 12.641_RLEN + 1.525_RLEN*dsqrtis + &
        dlog(ONE - 0.001005_RLEN*s)
  Kf = dexp(lnK)

  !---------------------------------------------------------------------
  ! kb = [H][BO2]/[HBO2] 
  ! pH scale: total
  !
  ! DOE (1994) using data from Dickson p.673 (1990)
  !---------------------------------------------------------------------
  lnK = (-8966.90_RLEN - 2890.53_RLEN*dsqrts - 77.942_RLEN*s + &
       1.728_RLEN*s15 - 0.0996_RLEN*s2)*invtk +                  &
       (148.0248_RLEN + 137.1942_RLEN*dsqrts + 1.62142_RLEN*s) +  &
       (-24.4344_RLEN - 25.085_RLEN*dsqrts - 0.2474_RLEN*s) *     &
       dlogtk + 0.053105_RLEN*dsqrts*tk
  Kb = dexp(lnK)



  select case ( mode)
  case ( DYNAMIC )
     small_interval = ( pH.gt.4.0_RLEN .and. pH.lt.9.0_RLEN) 
     if (small_interval) then
        h1 = 10.0_RLEN**(-(pH+0.5_RLEN))
        h2 = 10.0_RLEN**(-(pH-0.5_RLEN))
        Hplus = drtsafe2(h1, h2, XACC, error)
        !Hplus = drtsafe(h1, h2, XACC)
     end if
     if ((.not.small_interval) .or. error>0) then
        h1 = 10.0_RLEN**(-11.0_RLEN)
        h2 = 10.0_RLEN**(-2.0_RLEN)
        Hplus = drtsafe2(h1, h2, XACC, error)
        !Hplus = drtsafe(h1, h2, XACC)
     end if
     if ( error >0 ) then
        CalcCO2System=error
        return
     endif 
     !---------------------------------------------------------------
     ! Derive [CO2] as defined in DOE Methods Handbook 1994 Ver.2, 
     ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
     ! Compute other diagnostic variables (HCO3, CO3 and pH) and 
     ! pCO2, the CO2 partial pressure in the water (if way=1)
     !---------------------------------------------------------------
     Hplus2  =   Hplus*Hplus
     select case (way)
        case (1)
           CO2 = ldic*Hplus2/( Hplus2+ K1* Hplus+ K1* K2)
           pCO2_out  =   CO2/K0
           lpCO2     = pCO2_out ! defined for consistency with way 2
        case (2)
           CO2 = lpCO2*K0
           ldic = CO2*(Hplus2+K1*Hplus +K1*K2)/Hplus2
     end select
     HCO3  =   K1 * CO2 / Hplus
     CO3   =   K2 * HCO3 / Hplus
     pH    =  -dlog10(Hplus)

  case (FOLLOWS)
        !--------------------------------------------------
        ! Approximate solution by Follows et al., 2006
        ! First guess of H+: previous pH value
        !--------------------------------------------------
         xx= 10.0_RLEN**(-pH)
         xx2 = xx*xx
         xx3 = xx2*xx
         k11 = K1*K1
         k12 = K1*K2
         k12p = Kp(1)*Kp(2)
         k123p = k12p*Kp(3)
         a = xx3 + Kp(1)*xx2 + k12p*xx + k123p
         c = ONE + st/Ks
         cag = - bt/(ONE + xx/kb) - &
              kw/xx + xx -          &
              pt*k12p*xx/a -        &
              T2*pt*k123p/a -       &
              sit/(ONE + xx/ksi) +  &
              st/(ONE + ks/(xx/c))+ &
              ft/(ONE + kf/xx)+     &
              pt*xx3/a +            &
              ta
         gamm = ldic/cag
         dummy= (ONE-gamm)*(ONE-gamm)*k11 - 4.0_RLEN*k12*(ONE-2.0_RLEN*gamm)
         Hplus = 0.5_RLEN*((gamm-ONE)*K1 + dsqrt(dummy))
         !---------------------------------------------------------------
         ! Derive [CO2] as defined in DOE Methods Handbook 1994 Ver.2, 
         ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
         ! Compute other diagnostic variables (HCO3, CO3 and pH) and 
         ! pCO2, the CO2 partial pressure in the water (if way=1)
         !---------------------------------------------------------------
         Hplus2  =   Hplus*Hplus
         select case (way)
            case (1)
               CO2 = ldic*Hplus2/( Hplus2+ K1* Hplus+ K1* K2)
               pCO2_out  =   CO2/K0
               lpCO2     = pCO2_out ! defined for consistency with way 2
            case (2)
               CO2 = lpCO2*K0
               ldic = CO2*(Hplus2+K1*Hplus +K1*K2)/Hplus2
         end select
         HCO3  =   K1 * CO2 / Hplus
         CO3   =   K2 * HCO3 / Hplus
         pH    =  -dlog10(Hplus)

  case ( STATIC )
        !--------------------------------------------------
        ! Static approximate solution (for deep ocean)
        !--------------------------------------------------
        !-p/2(p) = tmp1:
        tmp1 = - 0.5_RLEN*( ta-DIC_in+ 4.0_RLEN* K2/ K1*( &
          2.0_RLEN* DIC_in- ta))/( 1.0_RLEN- 4.0_RLEN* K2/ K1)
        !q(p) = tmp2:
        tmp2 = - K2* (( 2.0_RLEN* DIC_in- ta))**(2.0_RLEN)/( &
          K1*( 1.0_RLEN- 4.0_RLEN* K2/ K1))
        !pCO2:
        pCO2_out = ( tmp1+ sqrt( tmp1* tmp1- tmp2))/ K0
        !CO2:
        CO2  =   K0 * pCO2_out
        !pH:
        pH = - dlog( K1* CO2/(2.0_RLEN* DIC_in &
          - ta- 2.0_RLEN* CO2))/ dlog( 10.0_RLEN)
        !CO3:
        CO3  =   ta - DIC_in + CO2
        !HCO3:
        HCO3  =   DIC_in - CO2 - CO3
  end select

  !---------------------------------------------------------------
  ! NOTE: transformation from mol/kg -----> umol/kg
  !       only for output purposes
  ! pco2 is converted from atm  -----> uatm
  !---------------------------------------------------------------
  if (present(DIC_out)) DIC_out  = ldic * MEG
  ta   = ta * MEG
  CO2  = co2 * MEG
  HCO3 = hco3 * MEG
  CO3  = co3 * MEG
  if (present(pCO2_out)) pCO2_out = lpCO2 * MEG
#ifdef DEBUG
  LEVEL3 'pco2',pco2_out
  LEVEL3 'pH',pH
  LEVEL3 'co2',co2
  LEVEL3 'K0',K0
  LEVEL3 'k1',K1
  LEVEL3 'k2',K2
  LEVEL3 'Kp',Kp(:)
  LEVEL3 'Kw',Kw
  LEVEL3 'Ksi',Ksi
  LEVEL3 'Ks',Ks
  LEVEL3 'Kf',Kf
#endif

  ! return no error
  CalcCO2System = 0

  end function CalcCO2System
!EOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !IROUTINE: CalcHplus
!
! DESCRIPTION
! This routine expresses TA as a function of DIC, hSWS (H+ on
! sea water scale) and constants.
! It also calculates the derivative of this function with respect to
! hSWS. It is used in the iterative solution for hSWS. In the call
! "x" is the input value for hSWS, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhSWS
!
!   INTENT(IN)  lx = pH
!   INTENT(OUT) fn = calculated value for TA
!   INTENT(OUT) df = calculated value for dTA/dHtotal
!
!       fn = hco3(x) + co3(x) + borate(x) + oh(x) + hpo4(x) +
!            2*po4(x) + silicate(x) + hfree(x) + hso4(x) +
!            hf(x) + h3po4(x) - ta
!
!       df = dfn/dx
!
! !INTERFACE
  subroutine CalcHplus(x,fn,df)
!
! !USES:
  use global_mem, ONLY:RLEN
  IMPLICIT NONE
! !INPUT PARAMETERS:
    real(RLEN),intent(IN)  :: x
! !OUTPUT PARAMETERS:
    real(RLEN),intent(OUT) :: fn,df
!
! AUTHORS
!   the OCMIP team
!   Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
!   - by A. Mouchet, 2004:
!   Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
!    1) fixed errors in signs, parenthesis and coefficient c in derivative
!    2) changed from Total to Seawater Scale
!       * c defined for seawater H scale;
!       * fn and df adapted to KF on free H scale
!       * comments have been adapted
! 
!EOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: lfn
  real(RLEN)  :: ldf
  real(RLEN)  :: x2,x3,k12,k12p,k123p,c,a,a2,da,b,b2,db


  !------------------------------------------------------------------------
  ! Derive H+ and other constants
  ! ---------------------------------------------------------------------
  x2    =   x* x
  x3    =   x2* x
  k12   =   K1* K2
  k12p  =   Kp(1)* Kp(2)
  k123p =   k12p* Kp(3)
  c  =   T1 + st/Ks + ft/Kf
  a  =   x3 + Kp(1)*x2 + k12p*x + k123p
  a2 =   a*a
  da =   T3*x2 + T2*Kp(1)*x + k12p
  b  =   x2 + K1*x + k12
  b2 =   b*b
  db =   T2*x+ K1

  !------------------------------------------------------------------------
  ! This computation depends on the input of CO2System above
  ! ---------------------------------------------------------------------
  select case ( way ) 
    case (1)  ! DIC and ALK
      !
      !	fn = hco3+co3
      !	dfn = dhco3/dx+co3/dx
      !
       fn = K1*x*ldic/b +        &
            T2*ldic*k12/b +      &
             bt/(T1 + x/Kb) 
       df = ((K1*ldic*b) - K1*x*ldic*db)/b2  &
            - T2*ldic*K12*db/b2              &
            - bt/Kb/(T1+x/Kb)**T2
    case (2) !pCO2 and ALK
       fn =  K0*(K1 * lpCO2/x  +    k12 * lpCO2/x2)
       df = -K0*(K1 * lpCO2/x2 + T2*k12 * lpCO2/x3)
  end select

  ! 0 = f(HCO3,CO3)+f(B,OH,HPO4,PO4,H+,Si(OH)) - TALK
        fn = fn +  Kw/x +        &
             pt*k12p*x/a +       &
             T2*pt*k123p/a +     &
             sit/(T1 + x/Ksi) -  &
             x/c -               &
             st/(T1 + Ks/(x/c))- &
             ft/(T1 + Kf/(x/c))- &
             pt*x3/a -           &
             ta

         df = df - Kw/x2                            &
             + (pt*k12p*(a - x*da))/a2              &
             - T2*pt*k123p*da/a2                    &
             - sit/Ksi/(T1+x/ksi)**T2               &
             - T1/c                                 &
             - st*(T1 + Ks/(x/c))**(-T2)*(Ks*c/x2)  &
             - ft*(T1 + Kf/(x/c))**(-T2)*(Kf*c/x2)    &
             - pt*x2*(T3*a-x*da)/a2

  end subroutine CalcHplus
!EOC


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcK0
!
! DESCRIPTION
! Computes K0, solubility of CO2 in the water (K Henry) from Weiss 1974 
! K0 = [CO2]/pCO2 [mol kg-1 atm-1]
! Defined as a function because it is called by other routines
!
! !INTERFACE
  function CalcK0(salt,temp)
!
! !USES:
  use constants, ONLY: ZERO_KELVIN
! !INPUT:
!
  real(RLEN),intent(IN)  :: salt
  real(RLEN),intent(IN)  :: temp
! !OUTPUT:
  real(RLEN) :: CalcK0
!EOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOC

! !LOCAL
  real(RLEN) :: tk100,tk1002

  tk100  =  (temp - ZERO_KELVIN)/ 100.0_RLEN
  tk1002 = tk100*tk100

  CalcK0 = exp(93.4517_RLEN/tk100 - 60.2409_RLEN + 23.3585_RLEN * dlog(tk100) +   &
       salt * (.023517_RLEN - 0.023656_RLEN * tk100 + 0.0047036_RLEN * tk1002))

  end function CalcK0
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Drtsafe
!
! !INTERFACE:
function drtsafe2(x1,x2,xacc,error)
!
! !DESCRIPTION:
!   Find roots of the Total Alkalinity function CalcHplus (see this module)
!   by Newton-Raphson and bisection
!   Adapted and optimized from Numerical Recipes drtsafe.f90 
!
! !USES:
   use global_mem, ONLY: RLEN
   use mem_CO2
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(IN)  :: x1,x2,xacc
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   integer,intent(OUT)    :: error
   logical                :: ready
   real(RLEN)             :: drtsafe2
!
! !REVISION HISTORY:
!  Author(s):
!   WH Press, SA Teukolsky, WT Vetterling and BP Flannery
!   Adapted from OCMIP standard files by M. Vichi and P. Ruardij
!
! !LOCAL VARIABLES:

   real(RLEN)             :: swap  
   real(RLEN)             :: df,dx,dxold,f,fh,fl,temp,xh,xl
   integer                :: j
   integer,parameter      :: MAXIT=100
!
!EOP
!-----------------------------------------------------------------------
!BOC

   error=0
   call CalcHplus(x1,fl,df)
   call CalcHplus(x2,fh,df)

   if (fl == ZERO) then
     drtsafe2=x1
     error = 1
     return
   else if (fh == ZERO) then
     drtsafe2=x2
     error = 1
     return
   else if (fl < ZERO) then
     xl=x1
     xh=x2
   else
     xh=x1
     xl=x2
     swap=fl
     fl=fh
     fh=swap
   end if

   drtsafe2=0.5_RLEN*(x1+x2)
   dxold=dabs(x2-x1)
   dx=dxold
   call CalcHplus(drtsafe2,f,df)

   j = 0
   ready = .FALSE.
   do while ( .NOT.ready .AND. j<MAXIT)
      j = j+1
      if (((drtsafe2-xh)*df-f)*((drtsafe2-xl)*df-f) >= ZERO .OR. &
         dabs(2.0_RLEN*f) > dabs(dxold*df) ) then
         dxold=dx
         dx=0.5_RLEN*(xh-xl)
         drtsafe2=xl+dx
         ready = (xl == drtsafe2) 
      else
         dxold=dx
         dx=f/df
         temp=drtsafe2
         drtsafe2=drtsafe2-dx
         ready = (temp == drtsafe2)
      end if
      ready = (dabs(dx) < xacc) 
      if (.NOT.ready) then
         call CalcHplus(drtsafe2,f,df)
         if (f < ZERO) then
            xl=drtsafe2
            fl=f
         else
            xh=drtsafe2
            fh=f
         end if
      end if
   end do
   if ( j.ge.MAXIT) error=2

   return

!-----------------------------------------------------------------------
end function drtsafe2
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Drtsafe
!
! !INTERFACE:
function drtsafe(x1,x2,xacc)
!
! !DESCRIPTION:
!   find roots of the Total Alkalinity function ta_iter_1 
!   by Newton-Raphson and bisection
!   Adapted and optimized from Numerical Recipes rtsafe.f90 
!   (error checking have been removed)
!
! !USES:
   use mem_co2
   use global_mem
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(IN)  :: x1,x2,xacc
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   real(RLEN)             :: drtsafe
!
! !REVISION HISTORY:
!  Author(s):
!   WH Press, SA Teukolsky, WT Vetterling and BP Flannery
!   Adapted from OCMIP standard files by M. Vichi
!
! !LOCAL VARIABLES:
   real(RLEN)             :: swap  
   real(RLEN)             :: df,dx,dxold,f,fh,fl,temp,xh,xl
   integer                :: j
   integer,parameter      :: MAXIT=100
!
!EOP
!-----------------------------------------------------------------------
!BOC


   call Ta_iter_1(x1,fl,df)
   call Ta_iter_1(x2,fh,df)

   if (fl == ZERO) then
     drtsafe=x1
     return
   else if (fh == ZERO) then
     drtsafe=x2
     return
   else if (fl < ZERO) then
     xl=x1
     xh=x2
   else
     xh=x1
     xl=x2
     swap=fl
     fl=fh
     fh=swap
   end if

   drtsafe=0.5_RLEN*(x1+x2)
   dxold=abs(x2-x1)
   dx=dxold
   call Ta_iter_1(drtsafe,f,df)

do j=1,MAXIT
  if (((drtsafe-xh)*df-f)*((drtsafe-xl)*df-f) >= ZERO .OR. &
       abs(2.0_RLEN*f) > abs(dxold*df) ) then
     dxold=dx
     dx=0.5_RLEN*(xh-xl)
     drtsafe=xl+dx
     if (xl == drtsafe) return
  else
     dxold=dx
     dx=f/df
     temp=drtsafe
     drtsafe=drtsafe-dx
     if (temp == drtsafe) return
  end if
  if (abs(dx) < xacc) return
  call Ta_iter_1(drtsafe,f,df)
  if (f < ZERO) then
     xl=drtsafe
     fl=f
  else
     xh=drtsafe
     fh=f
  end if
end do

!-----------------------------------------------------------------------
end function drtsafe
!-----------------------------------------------------------------------
!EOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE:  ta_iter_1
!
! !DESCRIPTION
! This routine expresses TA as a function of DIC, hSWS (H+ on
! sea water scale) and constants.
! It also calculates the derivative of this function with respect to
! hSWS. It is used in the iterative solution for hSWS. In the call
! "x" is the input value for hSWS, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhSWS
!
!   INTENT(IN) x = H+ total on Sea Water Scale, 
!   INTENT(OUT) fn = calculated value for TA
!   INTENT(OUT) df = calculated value for dTA/dHtotal
!
!       fn = hco3(x) + co3(x) + borate(x) + oh(x) + hpo4(x) +
!            2*po4(x) + silicate(x) + hfree(x) + hso4(x) +
!            hf(x) + h3po4(x) - ta
!
!       df = dfn/dx
!
! !INTERFACE
  subroutine ta_iter_1(x,fn,df)
!
! !USES
  use mem_CO2
!
! AUTHORS
!   the OCMIP team
!   Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
!   - by A. Mouchet, 2004:
!   Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
!    1) fixed errors in signs, parenthesis and coefficient c in derivative
!    2) changed from Total to Seawater Scale
!       * c defined for seawater H scale;
!       * fn and df adapted to KF on free H scale
!       * comments have been adapted
! 
!EOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
IMPLICIT NONE

    real(RLEN),intent(IN)  :: x
    real(RLEN),intent(OUT) :: fn,df
    real(RLEN)             :: x2,x3,k12,k12p,k123p,c,a,a2,da,b,b2,db  
    real(RLEN),parameter   :: T1=1.0_RLEN,T2=2.0_RLEN,T3=3.0_RLEN

        x2 = x*x
        x3 = x2*x
  k12   =   K1* K2
  k12p  =   Kp(1)* Kp(2)
  k123p =   k12p* Kp(3)
  c  =   T1 + st/Ks + ft/Kf
  a  =   x3 + Kp(1)*x2 + k12p*x + k123p
  a2 =   a*a
  da =   T3*x2 + T2*Kp(1)*x + k12p
  b  =   x2 + K1*x + k12
  b2 =   b*b
  db =   T2*x+ K1

        fn = k1*x*ldic/b +        &
             T2*ldic*k12/b +        &
             bt/(T1 + x/kb) +    &
             kw/x +                &
             pt*k12p*x/a +         &
             T2*pt*k123p/a +       &
             sit/(T1 + x/ksi) -  &
             x/c -                   &
             st/(T1 + ks/(x/c))- &
             ft/(T1 + kf/(x/c))- &
             pt*x3/a -             &
             ta

        df = ((k1*ldic*b) - k1*x*ldic*db)/b2 -     &
             T2*ldic*k12*db/b2 -                        &
             bt/kb/(T1+x/kb)**T2 -                 &
             kw/x2 +                                   &
             (pt*k12p*(a - x*da))/a2 -                 &
             T2*pt*k123p*da/a2 -                       &
             sit/ksi/(T1+x/ksi)**T2 -              &
             T1/c -                                      &
             st*(T1 + ks/(x/c))**(-T2)*(ks*c/x2) - &
             ft*(T1 + kf/(x/c))**(-T2)*(kf*c/x2) -   &
             pt*x2*(T3*a-x*da)/a2

  return

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end subroutine ta_iter_1
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!EOC
!EOC

end module CO2System
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
