#include "cppdefs.h"
#include "DEBUG.h"
#include "INCLUDE.h"

#if defined INCLUDE_PELCO2 || defined INCLUDE_BENCO2

module CO2System
!/*
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !MODULE: CO2System
!
! DESCRIPTION
!       Calculate co2 equilibrium in seawater starting from 
!       total alkalinity and total co2 (dic) at
!       a given temperature, salinity and other species
!       Routines are wrapped in a module in order to be used for pelagic
!       and benthic computations
!
!       INPUTS (BFM variables)
!       dic_in (O3c) = total dissolved inorganic carbon (umol C/kg) 
!       alk    (O3h) = alkalinity (umol eq/kg)
!       N1p     = inorganic phosphate (mmol/m3) 
!       N5s     = inorganic silicate (mmol/m3) 
!       ETW     = temperature (degrees C)
!       ESW     = salinity (PSU)
!       ERHO    = density (kg/m3)
!       press   = water pressure of the cell (dbar)
!
!       diagnostic OUTPUTS
!       co2  = co2*(aq) (umol/kg)
!       pco2 = oceanic pco2 (uatm)
!       hco3 = bicarbonate (umol/kg)
!       co3  = carbonate (umol/kg)
!       Ocalc = saturation for Calcite
!       Oarag = saturation for Aragonite
!       The code has 3 options defined in mod_pelco2:
!       choice of the equilibrium constants K1K2
!       choice of the ph scale PHSCALE
!       choice of the computation modeco2
!
! ORIGINAL REFERENCES
!--------------------------------------------------------------------------
!   J. Orr (LODYC) and the OCMIP team 
!   co2calc.f
!   Revision: 1.8  Date: 1999/07/16 11:40:33 
!---------------------------------------------------------------------
!   Zeebe & Wolf-Gladrow, co2 in Seawater: 
!               Equilibrium, Kinetics, Isotopes. 2001. Elsevier
!               Matlab scripts: csys.m; equic.m
!---------------------------------------------------------------------
!*/
!
! AUTHORS
!   M. Vichi and L. Patara (INGV-CMCC, vichi@bo.ingv.it)
!   P. Ruardij and H. Thomas (NIOZ, rua@nioz.nl)
! 
! CHANGE_LOG
!
! COPYING
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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


! !USES:
  use global_mem, ONLY: RLEN, LOGUNIT

! Shared variables
  implicit none
  real(RLEN),public  ::  K0 ! solubility : [Co2]=k0Ac Pco2
  real(RLEN),public  ::  K1 ! carbonate equilibrium I
  real(RLEN),public  ::  K2 ! carbonate equilibrium II
  real(RLEN),public  ::  Kw ! water dissociation
  real(RLEN),public  ::  Kb ! constant for Boron equilibrium
  real(RLEN),public  ::  Ks ! constant for bisulphate equilibrium
  real(RLEN),public  ::  Kf ! constant for hidrogen fluoride equilibirum
  real(RLEN),public  ::  Kp(3) ! constants for phosphate equilibirum
  real(RLEN),public  ::  Ksi ! constant for silicic acid eqilbrium


  real(RLEN),public  :: press=0.0_RLEN ! local pressure in bar 
  real(RLEN),public  :: ldic  ! local dic in mol/kg 
  real(RLEN),public  :: lpco2 ! local pco2
  real(RLEN),public  :: scl   ! chlorinity
  real(RLEN),public  :: ta    ! total alkalinity
  integer,public     :: way

  real(RLEN)         :: bt,ft,st,pt,sit
  real(RLEN),parameter   :: T1=1.0_RLEN,T2=2.0_RLEN,T3=3.0_RLEN
  
  public CalcCO2System, CalcHplus, CalcK0, drtsafe2
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  functions 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  contains

!-------------------------------------------------------------------------!
!BOP
! !IROUTINE:  CalcCO2System
!
! !INTERFACE
  integer function CalcCO2System(mode,salt,temp,rho,n1p,n5s,alk,co2,hco3,co3,ph, &
                         pr_in,dic_in,pco2_in,dic_out,pco2_out,omegacal,omegarag)
! !DESCRIPTION
! See module preamble
!
! USES:
  use global_mem, ONLY: ONE,ZERO
  use constants, ONLY: ZERO_KELVIN,MW_C,Rgas
  use mem_co2
  
!  use mem_co2

! !INPUT PARAMETERS:
  IMPLICIT NONE
  integer                          :: mode
  real(RLEN),intent(IN)            :: salt
  real(RLEN),intent(IN)            :: temp
  real(RLEN),intent(IN)            :: rho 
  real(RLEN),intent(IN)            :: n1p
  real(RLEN),intent(IN)            :: n5s
  real(RLEN),intent(IN)            :: alk
  real(RLEN),intent(IN),optional   :: pr_in
  real(RLEN),intent(IN),optional   :: dic_in
  real(RLEN),intent(IN),optional   :: pco2_in
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
  real(RLEN),intent(OUT)           :: co2
  real(RLEN),intent(OUT)           :: hco3
  real(RLEN),intent(OUT)           :: co3
  real(RLEN),intent(INOUT)         :: ph
  real(RLEN),intent(OUT),optional  :: dic_out
  real(RLEN),intent(OUT),optional  :: pco2_out
  real(RLEN),intent(OUT),optional  :: omegacal
  real(RLEN),intent(OUT),optional  :: omegarag

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),parameter             :: MEG=1.E6_RLEN, &
                                      PERMIL=ONE/1000.0_RLEN,&
                                      PERMEG=ONE/MEG

  logical     :: small_interval
  integer     :: i, l, error
  real(RLEN)  :: Hplus, Hplus2
  real(RLEN)  :: tmp1, tmp2
  real(RLEN)  :: intercept,lnk,sit,pt
  real(RLEN)  :: h1,h2
  real(RLEN)  :: tk,tk100,tk1002,    &
                 tc2,pr,pr2,         &
                 invtk,is,is2,       &
                 dlogtk,dsqrtis,s,s2,&
                 dsqrts,s15
  ! Variables for Follows' parameterization
  real(RLEN) :: xx,xx2,xx3,k11,k12,k12p,k123p,k3p,a,c
  real(RLEN) :: cag,dummy,gamm
  ! Variables for Calcite and Aragonite omegas computation
  real(RLEN) :: Caconc, kspcal, kspcal2, ksparag, ksparag2,lomegacal,lomegarag
  ! Pressure correction
  real(RLEN),parameter,dimension(11) :: &
  a0 = (/ 25.5_RLEN,15.82_RLEN,29.48_RLEN,25.60_RLEN,18.03_RLEN, &
          9.78_RLEN,48.76_RLEN,46.0_RLEN,14.51_RLEN,23.12_RLEN,26.57_RLEN /), &
  a1 = (/ 0.1271_RLEN,-0.0219_RLEN,0.1622_RLEN,0.2324_RLEN,0.0466_RLEN, &
          -0.0090_RLEN, 0.5304_RLEN,  0.5304_RLEN,0.1211_RLEN,0.1758_RLEN,0.2020_RLEN/), &
  a2 = (/ 0.0_RLEN,0.0_RLEN,2.608_RLEN,-3.6246_RLEN,0.316_RLEN,-0.942_RLEN, &
          0.0_RLEN,0.0_RLEN,-0.321_RLEN,-2.647_RLEN,-3.042_RLEN /), &
  b0 = (/ 3.08_RLEN,-1.13_RLEN,2.84_RLEN,5.13_RLEN,4.53_RLEN,3.91_RLEN, &
          11.76_RLEN,11.76_RLEN,2.67_RLEN,5.15_RLEN,4.08_RLEN /), &
  b1 = (/ 0.0877_RLEN,-0.1475_RLEN,0.0_RLEN,0.0794_RLEN,0.09_RLEN,0.054_RLEN, &
          0.3692_RLEN,0.3692_RLEN,0.0427_RLEN,0.09_RLEN,0.0714_RLEN /)
  real(RLEN) :: deltav,deltak
  real(RLEN),dimension(11) :: lnkpok0
!EOP
!-------------------------------------------------------------------------!
!BOC

  !---------------------------------------------------------------------
  ! Change units from the input of mmol/m^3 -> mol/kg:
  ! (1 mmol/m^3)  x (1 m^3/1024.5 kg) / 1000
  ! The ocean actual density /*ERHO*/ is used.
  ! Note: mol/kg are actually what the body of this routine uses 
  ! for all calculations.  
  !---------------------------------------------------------------------
  pt = n1p/rho*PERMIL
  sit = n5s/rho*PERMIL
  if (present(dic_in)) then
     ! convert from umol/kg to standard dic units mol/kg 
     ldic = dic_in*PERMEG
     way = 1 
  elseif (present(pco2_in)) then
     lpco2  = pco2_in
     way = 2 
  endif
  ! convert input variable alkalinity from umol/kg to mol/kg 
  ta = alk * PERMEG
  ! convert from dbar to bar
  if (present(pr_in)) press = pr_in * 0.1_RLEN
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
  tk = temp - ZERO_KELVIN  !ZERO_KELVIN=-273.16; ETW is in degC; tk is in degK
  tc2 = temp * temp
#ifdef DEBUG
  LEVEL2 'Entering CO2System...'
  LEVEL3 'temp',temp,'sal',salt
  LEVEL3 'tk',tk,'way',way,'mode',mode
  LEVEL3 'dic',ldic,'ta',ta
  LEVEL3 'pt',pt,'sit',sit
#endif
  tk100 = tk/100.0_RLEN
  tk1002 = tk100*tk100
  invtk = ONE/tk
  dlogtk = log(tk)
  ! salinity
  s  = salt
  s2 = salt*salt
  dsqrts = sqrt(salt)
  s15 = salt**1.5_RLEN
  ! chlorinity
  scl = salt/1.80655_RLEN  
  ! ionic strength
  is = 19.924_RLEN*salt/ (1000._RLEN-1.005_RLEN*salt)
  is2 = is*is
  dsqrtis = sqrt(is)
  ! divide water pressure by the gas constant
  if (present(pr_in)) then 
     pr2 = press * press / Rgas
     pr = press / Rgas
  endif
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
!TOM Should we keep it to compare with SOCAT
  ! ---------------------------------------------------------------------
  ! convert partial pressure to fugacity
  ! ff = k0(1-ph2O)*correction term for non-ideality
  !
  ! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
  ! ---------------------------------------------------------------------
  !ff = exp(-162.8301_RLEN + 218.2968_RLEN/tk100  +      &
  !     90.9241_RLEN*log(tk100) - 1.47696_RLEN*tk1002 +  &
  !     s * (.025695_RLEN - .025225_RLEN*tk100 +         &
  !     0.0049867_RLEN*tk1002))

  ! ---------------------------------------------------------------------
  ! K0, solubility of co2 in the water (K Henry)
  ! from Weiss 1974; K0 = [co2]/pco2 [mol kg-1 atm-1]
  ! ---------------------------------------------------------------------
  K0 = CalcK0(salt,temp)

  ! ---------------------------------------------------------------------
  ! Choice of Acidity constants from parameter in co2 namelist
  ! k1 = [H][hco3]/[H2co3]
  ! k2 = [H][co3]/[hco3]
  ! If Using Follows et al. parameterization, force K1K2=1
  ! ---------------------------------------------------------------------
  select case (K1K2)
  case (1)
     ! ---------------------------------------------------------------------
     ! Constants according to Roy et al. (1993a). 
     ! Recommended by DOE (1994) and Zeebe / Wolf-Gladrow. 
     ! ph scale: total 
     ! ---------------------------------------------------------------------
     lnK = -2307.1266_RLEN*invtk +2.83655_RLEN-1.5529413_RLEN*dlogtk +  &
          (-4.0484_RLEN*invtk - 0.20760841_RLEN) * dsqrts +            &
          0.08468345_RLEN* s -0.00654208_RLEN * s15+log(ONE -0.001005_RLEN* s )
     K1 = exp(lnK)
     lnK = -3351.6106_RLEN/tk -9.226508_RLEN-0.2005743_RLEN*dlogtk+        &
          (-23.9722_RLEN/tk - 0.10690177_RLEN)*dsqrts + 0.1130822_RLEN* s - &
          0.00846934_RLEN*s15 + log(ONE-0.001005_RLEN* s )
     K2 = exp(lnK)
  case (2)
     ! ---------------------------------------------------------------------
     ! Mehrbach et al. (1973) as refitted by Dickson and Millero (1987) 
     ! ph scale: seawater  (Millero, 1995, p.664)
     ! Standard OCMIP computation. Natural seawater
     ! ---------------------------------------------------------------------
     K1 = 10.0_RLEN**(-ONE*(3670.7_RLEN*invtk - 62.008_RLEN + 9.7944_RLEN*dlogtk - &
          0.0118_RLEN * s + 0.000116_RLEN*s2))
     K2 = 10.0_RLEN**(-1.0*(1394.7_RLEN*invtk + 4.777_RLEN - &
          0.0184_RLEN*s + 0.000118_RLEN*s2))
  case (3)
     ! ---------------------------------------------------------------------
     ! Mehrbach et al (1973) refit by Lueker et al. (2000). 
     ! ph scale: total 
     ! ---------------------------------------------------------------------
     K1 = 10.0_RLEN**(-ONE*(3633.86_RLEN*invtk - 61.2172_RLEN + 9.6777_RLEN*dlogtk - &
          0.011555_RLEN * s + 0.0001152_RLEN * s2))
     K2 = 10.0_RLEN**(-ONE*(471.78_RLEN*invtk + 25.9290_RLEN - &
          3.16967_RLEN*dlogtk - 0.01781_RLEN * s + 0.0001122_RLEN * s2))
  case (4)
     !-----------------------------------------------------------------------
     ! Hansson (1973b) data as refitted by Dickson and Millero (1987).
     ! ph scale: seawater 
     !-----------------------------------------------------------------------
     K1 = 10.0_RLEN**(-ONE*(851.4_RLEN*invtk + 3.237_RLEN - &
          0.0106_RLEN*s + 0.000132_RLEN*s2))
     K2 = 10.0_RLEN**(-ONE*(-3885.4_RLEN*invtk + 125.844_RLEN - 18.141_RLEN*dlogtk -  &
          0.0192_RLEN*s + 0.000132_RLEN* s2))

  end select

  ! ---------------------------------------------------------------------
  ! k1p = [H][H2PO4]/[H3PO4] 
  ! ph scale: total
  !
  ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
  ! Millero p.670 (1995)
  ! ---------------------------------------------------------------------
  lnK = -4576.752_RLEN*invtk + 115.525_RLEN - 18.453_RLEN * dlogtk + &
       (-106.736_RLEN*invtk + 0.69171_RLEN) * dsqrts +               &
       (-0.65643_RLEN*invtk - 0.01844_RLEN) * s
  Kp(1) = exp(lnK)
  ! ---------------------------------------------------------------------
  ! k2p = [H][HPO4]/[H2PO4] 
  ! ph scale: total
  !
  ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
  ! Millero p.670 (1995)
  ! ---------------------------------------------------------------------
  lnK = -8814.715_RLEN*invtk + 172.0883_RLEN - 27.927_RLEN * dlogtk + &
       (-160.340_RLEN*invtk + 1.3566_RLEN) * dsqrts +                 &
       (0.37335_RLEN*invtk - 0.05778_RLEN) * s
  Kp(2) = exp(lnK)
  !------------------------------------------------------------------------
  ! k3p = [H][PO4]/[HPO4] 
  ! ph scale: total
  !
  ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  ! ---------------------------------------------------------------------
  lnK = -3070.75_RLEN*invtk - 18.126_RLEN + &
       (17.27039_RLEN*invtk + 2.81197_RLEN) *   &
       dsqrts + (-44.99486_RLEN*invtk - 0.09984_RLEN) * s
  Kp(3) = exp(lnK)

  !------------------------------------------------------------------------
  ! ksi = [H][SiO(OH)3]/[Si(OH)4] ph on Sea Water Scale
  ! ph scale: total
  !
  ! DOE(1994) as recommended by Zeebe and Wolf-Gladrow
  ! ---------------------------------------------------------------------
  ksi = -8904.2_RLEN*invtk + 117.385_RLEN - 19.334_RLEN * dlogtk + &
       (-458.79_RLEN*invtk + 3.5913_RLEN) * dsqrtis +              &
       (188.74_RLEN*invtk - 1.5998_RLEN) * is +                   &
       (-12.1652_RLEN*invtk + 0.07871_RLEN) * is2 +               &
       log(ONE-0.001005_RLEN*s)
  Ksi = exp(lnK)

  !------------------------------------------------------------------------
  ! kw = [H][OH]  ion product of water
  ! ph scale: see below
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
  Kw = exp(lnK)

  !------------------------------------------------------------------------
  ! ks = [H][SO4]/[HSO4] 
  ! ph scale: "free"
  !
  ! Dickson (1990, J. chem. Thermodynamics 22, 113)
  ! ---------------------------------------------------------------------
  lnK = -4276.1_RLEN*invtk + 141.328_RLEN - 23.093_RLEN*dlogtk +      &
       (-13856._RLEN*invtk + 324.57_RLEN - 47.986_RLEN*dlogtk) * dsqrtis + &
       (35474._RLEN*invtk - 771.54_RLEN + 114.723_RLEN*dlogtk) * is -     &
       2698._RLEN*invtk*is**1.5_RLEN + 1776._RLEN*invtk*is2 +              &
       log(ONE - 0.001005_RLEN*s)
  Ks = exp(lnK)

  !------------------------------------------------------------------------
  ! kf = [H][F]/[HF] ph on free scale
  ! ph scale: "free"
  !
  ! Dickson and Riley (1979)  also Dickson and Goyet (1994)
  ! ---------------------------------------------------------------------
  lnK = 1590.2_RLEN*invtk - 12.641_RLEN + 1.525_RLEN*dsqrtis + &
        log(ONE - 0.001005_RLEN*s)
  Kf = exp(lnK)

  !---------------------------------------------------------------------
  ! kb = [H][BO2]/[HBO2] 
  ! ph scale: total
  !
  ! DOE (1994) using data from Dickson p.673 (1990)
  !---------------------------------------------------------------------
  lnK = (-8966.90_RLEN - 2890.53_RLEN*dsqrts - 77.942_RLEN*s + &
       1.728_RLEN*s15 - 0.0996_RLEN*s2)*invtk +                  &
       (148.0248_RLEN + 137.1942_RLEN*dsqrts + 1.62142_RLEN*s) +  &
       (-24.4344_RLEN - 25.085_RLEN*dsqrts - 0.2474_RLEN*s) *     &
       dlogtk + 0.053105_RLEN*dsqrts*tk
  Kb = exp(lnK)
  !------------------------------------------------------------------------
  ! Effect of pressure on equilibrium constant (Millero 1995)
  ! Zeebe and Wolf-Gladrow, CO2 in Seawater, Elsevier, 2001
  ! Pag. 267
  !------------------------------------------------------------------------
  if (present(pr_in) .AND. press .GT. 0) then
     do l=1,11
        deltav  =  -a0(l) + a1(l)*temp + 1.0e-3_RLEN*a2(l)*tc2
        deltak  = 1.0e-3_RLEN*(-b0(l) + b1(l)*temp)
        lnkpok0(l) = -deltav*invtk*pr + 0.5_RLEN*deltak*invtk*pr2
     end do
     k1 = k1*exp(lnkpok0(1))
     k2 = k2*exp(lnkpok0(2))
     kb = kb*exp(lnkpok0(3))
     kw = kw*exp(lnkpok0(4))
     ks = ks*exp(lnkpok0(5))
     kf = kf*exp(lnkpok0(6))
     kp(1) = kp(1)*exp(lnkpok0(9))
     kp(2) = kp(2)*exp(lnkpok0(10))
     kp(3) = kp(3)*exp(lnkpok0(11))
  endif 
  !------------------------------------------------------------------------
  !
  ! Calculate [H+] total when DIC and TA are known.
  ! Available Methods: 1 - STATIC
  !                    
  !                    2 - FOLLOWS
  !
  !------------------------------------------------------------------------

  select case ( mode)
  case ( DYNAMIC )
     small_interval = ( ph.gt.4.0_RLEN .and. ph.lt.9.0_RLEN) 
     if (small_interval) then
        h1 = 10.0_RLEN**(-(ph+M2PHDELT))
        h2 = 10.0_RLEN**(-(ph-M2PHDELT))
        Hplus = drtsafe2(h1, h2, M2XACC, M2MAXIT, error)
     end if
     if ((.not.small_interval) .or. error>0) then
        h1 = 10.0_RLEN**(-11.0_RLEN)
        h2 = 10.0_RLEN**(-2.0_RLEN)
        Hplus = drtsafe2(h1, h2, M2XACC, M2MAXIT, error)
     end if
     if ( error >0 ) then
        CalcCO2System=error
        return
     endif 
     !---------------------------------------------------------------
     ! Derive [co2] as defined in DOE Methods Handbook 1994 Ver.2, 
     ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
     ! Compute other diagnostic variables (hco3, co3 and ph) and 
     ! pco2, the co2 partial pressure in the water (if way=1)
     !---------------------------------------------------------------
     Hplus2  =   Hplus*Hplus
     select case (way)
        case (1)
           co2 = ldic*Hplus2/( Hplus2+ K1* Hplus+ K1* K2)
           pco2_out  =   co2/K0
           lpco2     = pco2_out ! defined for consistency with way 2
        case (2)
           co2 = lpco2*K0
           ldic = co2*(Hplus2+K1*Hplus +K1*K2)/Hplus2
     end select
     hco3  =   K1 * co2 / Hplus
     co3   =   K2 * hco3 / Hplus
     ph    =  -log10(Hplus)

  case (FOLLOWS)
        !--------------------------------------------------
        ! Approximate solution by Follows et al., 2006
        ! First guess of H+: previous ph value
        !--------------------------------------------------
         xx= 10.0_RLEN**(-ph)
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
         Hplus = 0.5_RLEN*((gamm-ONE)*K1 + sqrt(dummy))
         !---------------------------------------------------------------
         ! Derive [co2] as defined in DOE Methods Handbook 1994 Ver.2, 
         ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
         ! Compute other diagnostic variables (hco3, co3 and ph) and 
         ! pco2, the co2 partial pressure in the water (if way=1)
         !---------------------------------------------------------------
         Hplus2  =   Hplus*Hplus
         select case (way)
            case (1)
               co2 = ldic*Hplus2/( Hplus2+ K1* Hplus+ K1* K2)
               pco2_out  =   co2/K0
               lpco2     = pco2_out ! defined for consistency with way 2
            case (2)
               co2 = lpco2*K0
               ldic = co2*(Hplus2+K1*Hplus +K1*K2)/Hplus2
         end select
         hco3  =   K1 * co2 / Hplus
         co3   =   K2 * hco3 / Hplus
         ph    =  -log10(Hplus)

  case ( STATIC )
        !--------------------------------------------------
        ! Static approximate solution (for deep ocean)
        !--------------------------------------------------
        !-p/2(p) = tmp1:
        tmp1 = - 0.5_RLEN*( ta-dic_in+ 4.0_RLEN* K2/ K1*( &
          2.0_RLEN* dic_in- ta))/( 1.0_RLEN- 4.0_RLEN* K2/ K1)
        !q(p) = tmp2:
        tmp2 = - K2* (( 2.0_RLEN* dic_in- ta))**(2.0_RLEN)/( &
          K1*( 1.0_RLEN- 4.0_RLEN* K2/ K1))
        !pco2:
        pco2_out = ( tmp1+ sqrt( tmp1* tmp1- tmp2))/ K0
        !co2:
        co2  =   K0 * pco2_out
        !ph:
        ph = - log( K1* co2/(2.0_RLEN* dic_in &
          - ta- 2.0_RLEN* co2))/ log( 10.0_RLEN)
        !co3:
        co3  =   ta - dic_in + co2
        !hco3:
        hco3  =   dic_in - co2 - co3
  end select

! ======================================================================
!  Calculate CaCO3 saturation state parameters
! ======================================================================
!
! Omega for calcite and aragonite (< 1 = undersaturation; > 1 = 
! supersaturation)
! [derived from Zeebe & Wolf-Gladrow Matlab routine; but apparently
!  originally from Mucci, 1983 (?)]
!
!     convert Ca concentration from mol/m3 to mol/kg
      Caconc  = Caconc0 / rho          ! concentration in mol/kg
!
!     Ca concentration normalised to salinity ! specified in namelist
      if (Canorm) Caconc = Caconc * (salt / 35.0_RLEN)
!
!     Ksp_calcite
      kspcal   = -171.9065_RLEN - 0.077993_RLEN * tk + 2839.319_RLEN * invtk &
           + 71.595_RLEN * log10(tk) + dsqrts * (-0.77712_RLEN +  &
           0.0028426_RLEN * tk + 178.34_RLEN * invtk)            &
           - 0.07711_RLEN * salt + 0.0041249_RLEN *s15
      kspcal2  = 10.0_RLEN**kspcal      ! in mol/kg
!
!     Ksp_aragonite
      ksparag  = -171.945_RLEN - 0.077993_RLEN * tk + 2903.293_RLEN * invtk   &
           + 71.595_RLEN * log10(tk) + dsqrts * (-0.068393_RLEN +   &
           0.0017276_RLEN * tk + 88.135_RLEN * invtk)   &
           - 0.10018_RLEN * salt + 0.0059415_RLEN * s15
      ksparag2 = 10.0_RLEN**ksparag     ! in mol/kg

!     Pressure effect (see above)
      if (present(pr_in) .AND. press .GT. 0) then
         kspcal2 = kspcal2*exp(lnkpok0(7));
         ksparag2 = ksparag2*exp(lnkpok0(8));
      endif 
!
!     Omega_calcite
      lomegacal = Caconc * co3 / kspcal2
!
!     Omega_aragonite
      lomegarag = Caconc * co3 / ksparag2
  !---------------------------------------------------------------
  ! NOTE: transformation from mol/kg -----> umol/kg
  !       only for output purposes
  ! pco2 is converted from atm  -----> uatm
  !---------------------------------------------------------------
  if (present(dic_out)) dic_out  = ldic * MEG
  ta   = ta * MEG
  co2  = co2 * MEG
  hco3 = hco3 * MEG
  co3  = co3 * MEG
  if (present(pco2_out)) pco2_out = lpco2 * MEG
  if (present(omegacal)) omegacal = lomegacal
  if (present(omegarag)) omegarag = lomegarag 
#ifdef DEBUG
  LEVEL3 'pco2',pco2_out
  LEVEL3 'ph',ph
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
! This routine expresses TA as a function of dic, hSWS (H+ on
! sea water scale) and constants.
! It also calculates the derivative of this function with respect to
! hSWS. It is used in the iterative solution for hSWS. In the call
! "x" is the input value for hSWS, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhSWS
!
!   INTENT(IN)  lx = ph
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
    case (1)  ! dic and ALK
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
    case (2) !pco2 and ALK
       fn =  K0*(K1 * lpco2/x  +    k12 * lpco2/x2)
       df = -K0*(K1 * lpco2/x2 + T2*k12 * lpco2/x3)
  end select

  ! 0 = f(hco3,co3)+f(B,OH,HPO4,PO4,H+,Si(OH)) - TALK
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
! Computes K0, solubility of co2 in the water (K Henry) from Weiss 1974 
! K0 = [co2]/pco2 [mol kg-1 atm-1]
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

  CalcK0 = exp(93.4517_RLEN/tk100 - 60.2409_RLEN + 23.3585_RLEN * log(tk100) +   &
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
function drtsafe2(x1,x2,xacc,maxit,error)
!
! !DESCRIPTION:
!   Find roots of the Total Alkalinity function CalcHplus (see this module)
!   by Newton-Raphson and bisection
!   Adapted and optimized from Numerical Recipes drtsafe.f90 
!
! !USES:
   use global_mem, ONLY: RLEN, ZERO
!   use mem_co2
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(IN)  :: x1,x2,xacc
   integer,intent(IN)     :: maxit
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
   dxold=abs(x2-x1)
   dx=dxold
   call CalcHplus(drtsafe2,f,df)

   j = 0
   ready = .FALSE.
   do while ( .NOT.ready .AND. j<maxit)
      j = j+1
      if (((drtsafe2-xh)*df-f)*((drtsafe2-xl)*df-f) >= ZERO .OR. &
         abs(2.0_RLEN*f) > abs(dxold*df) ) then
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
      ready = (abs(dx) < xacc) 
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
   if ( j.ge.maxit) error=2

   return

!-----------------------------------------------------------------------
end function drtsafe2
!-----------------------------------------------------------------------

end module CO2System
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
