#include "cppdefs.h"
#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE:  CO2Calc
!
! !DESCRIPTION
!	Calculate pCO2 in the water from total alkalinity and total CO2 at
!       a given model temperature and salinity 
!
!       INPUTS 
!	O3c 	= total dissolved inorganic carbon (mgC/m3) 
!	O3h  	= alkalinity (umol/kg)
!	N1p  	= inorganic phosphate (mmol/m^3) 
!	N5s 	= inorganic silicate (mmol/m^3) 
!	ETW     = temperature (degrees C)
!	ESW     = salinity (PSU)
!	ERHO    = density (kg/m3)
!	EPCO2air	= atmospheric CO2 partial pressure (ppm) 
!
!	phlo   	= lower limit of pH range
!	phhi   	= upper limit of pH range
!
!       diagnostic OUTPUTS
!	co2  = CO2*(aq) (umol/kg)
!       pco2 = oceanic pCO2 (uatm)
!       hco3 = bicarbonate (umol/kg)
!       co3  = carbonate (umol/kg)
!       co2airflux = air-sea CO2 flux (mmol/m^2/d)
!
!       The code has to options defined in mod_CO2:
!       choice of the equilibrium constants K1K2
!       choice of the computation CalcCO2ITER
!
! ORIGINAL REFERENCES
!--------------------------------------------------------------------------
!   J. Orr (LODYC) and the OCMIP team 
!   co2calc.f
!   Revision: 1.8  Date: 1999/07/16 11:40:33 
!---------------------------------------------------------------------
!   Zeebe & Wolf-Gladrow, CO2 in Seawater: 
!		Equilibrium, Kinetics, Isotopes. 2001. Elsevier
!               Matlab scripts: csys.m; equic.m
!---------------------------------------------------------------------
!
! AUTHOR
!   M. Vichi and L. Patara (INGV, vichi@bo.ingv.it)
! 
! CHANGE_LOG
!
! COPYING
!   Copyright (C) 2006 Marcello Vichi (vichi@bo.ingv.it)
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

subroutine CO2Calc

  use global_mem, ONLY:RLEN,ONE,ZERO_KELVIN
  use constants, ONLY: MW_C
  use mem, ONLY: NO_BOXES,NO_BOXES_XY,N1p,N5s,O3c,O3h,    &
                 ETW,ESW,EWIND,ERHO,EICE,EPCO2air,Depth, &
                 flux_vector,Erho,P1c,P2c,P3c,B1c,R6c,R1c,  &
                 Z4c,Z5c,Z6c,ppN1p,ppN5s,ppO3c,ppO3h, co2, co3,  &
                 hco3, pco2, ppco2, pppco2, pphco3, pH, pppH, &
                 dic, ppdic, CO2airflux,iiPel 
  use mem, ONLY: flux_vector
  use mem_CO2
  use api_bfm, ONLY: SRFindices               


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  implicit none

  interface
    function drtsafe(x1,x2,xacc) ! finds roots of the f(dic,ta,htotal,...)
      use constants, ONLY: RLEN
      IMPLICIT NONE
      real(RLEN) :: x1,x2,xacc
      real(RLEN) :: drtsafe
    end function
  end interface

  real(RLEN),parameter             :: MEG=1.D6,XACC=1.D-25,  &
                                      PERMIL=ONE/1000.0_RLEN,&
                                      PERMEG=ONE/MEG
  real(RLEN),parameter   :: T1=1.0_RLEN,T2=2.0_RLEN,T3=3.0_RLEN

  ! Local variables
  integer                          :: l
  real(RLEN),dimension(NO_BOXES) :: k0,k12,k1,k2,kw,kb,ks,kf, &
                                    k1p,k2p,k3p,k12p,k123p,ksi,ff
  real(RLEN),dimension(NO_BOXES) :: bt,st,ft,sit,pt,ta

  real(RLEN)                       :: x1,x2
  real(RLEN),dimension(NO_BOXES)   :: tk,tk100,tk1002,    &
                                      invtk,is,is2,       &
                                      dlogtk,dsqrtis,s,s2, &
                                      dsqrts,s15,scl,htotal2 
  real(RLEN)             :: tmpflux(NO_BOXES)

  real(RLEN),dimension(NO_BOXES) :: xx,xx2,xx3,k11,a,c
  real(RLEN),dimension(NO_BOXES) :: cag,dummy,gamm
  real(RLEN)                     :: phscale


  !---------------------------------------------------------------------
  ! Change units from the input of mmol/m^3 -> mol/kg:
  ! (1 mmol/m^3)  x (1 m^3/1024.5 kg) / 1000
  ! The ocean actual surface density ERHO is used.
  ! Note: mol/kg are actually what the body of this routine uses 
  ! for calculations.  
  !---------------------------------------------------------------------
  pt = N1p/ERHO*PERMIL
  sit= N5s/ERHO*PERMIL
  ! convert O3c from mg/m3 to standard DIC units mol/kg 
  dic = O3c/MW_C/ERHO*PERMIL  
  ! convert alkalinity O3h from umol/kg to mol/kg 
  ta = O3h*PERMEG
  !ta = O3h + 2210.0 ! Alternative, if ta is changed by NO3 uptake

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
  tk = ETW - ZERO_KELVIN
#ifdef DEBUG
  LEVEL2 'Entering CO2Calc...'
  LEVEL3 'tk',tk,'O3c',O3c
  LEVEL3 'dic',dic,'ta',ta
  LEVEL3 'pt',pt,'sit',sit
#endif
  tk100 = tk/100.0_RLEN
  tk1002 = tk100*tk100
  invtk = ONE/tk
  dlogtk = dlog(tk)
  is = 19.924_RLEN*ESW / (1000._RLEN-1.005_RLEN*ESW)
  is2 = is*is
  dsqrtis = dsqrt(is)
  s  = ESW
  s2 = ESW*ESW
  dsqrts = dsqrt(ESW)
  s15 = ESW**1.5_RLEN
  scl = ESW/1.80655_RLEN

  !------------------------------------------------------------------------
  ! Calculate concentrations for borate, sulfate, and fluoride
  ! ---------------------------------------------------------------------
  ! Uppstrom (1974)
  bt = 0.000232_RLEN * scl/10.811_RLEN
  ! Morris & Riley (1966)
  st = 0.14_RLEN * scl/96.062_RLEN
  ! Riley (1965)
  ft = 0.000067_RLEN * scl/18.9984_RLEN

  ! ---------------------------------------------------------------------
  ! ff = k0(1-pH2O)*correction term for non-ideality
  !
  ! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
  ! ---------------------------------------------------------------------
  ff = exp(-162.8301_RLEN + 218.2968_RLEN/tk100  +      &
       90.9241_RLEN*dlog(tk100) - 1.47696_RLEN*tk1002 +  &
       s * (.025695_RLEN - .025225_RLEN*tk100 +         &
       0.0049867_RLEN*tk1002))
   
  ! ---------------------------------------------------------------------
  ! K0, solubility of CO2 in the water (K Henry)
  ! from Weiss 1974; K0 = [CO2]/pCO2 [mol kg-1 atm-1]
  ! ---------------------------------------------------------------------
  K0 = exp(93.4517_RLEN/tk100 - 60.2409_RLEN + 23.3585_RLEN * dlog(tk100) +   &
       s * (.023517_RLEN - 0.023656_RLEN * tk100 + 0.0047036_RLEN * tk1002))

  ! ---------------------------------------------------------------------
  ! Choice of Acidity constants from parameter in CO2 namelist
  ! k1 = [H][HCO3]/[H2CO3]
  ! k2 = [H][CO3]/[HCO3]
  ! If Using Follows et al. parameterization, force K1K2=1
  ! ---------------------------------------------------------------------
  if (.NOT.CalcCO2ITER) K1K2=1
  select case (K1K2)
  case (1)
     ! ---------------------------------------------------------------------
     ! Constants according to Roy et al. (1993a). Recommended by DOE (1994). 
     ! pH total scale
     ! ---------------------------------------------------------------------
     k1 = exp(-2307.1266_RLEN*invtk +2.83655_RLEN-1.5529413_RLEN*dlogtk +  &
          (-4.0484_RLEN*invtk - 0.20760841_RLEN) * dsqrts +            &
          0.08468345_RLEN* s -0.00654208_RLEN * s15+dlog(ONE -0.001005_RLEN* s ))
     k2 = exp(-3351.6106_RLEN/tk -9.226508_RLEN-0.2005743_RLEN*dlogtk+        &
          (-23.9722_RLEN/tk - 0.10690177_RLEN)*dsqrts + 0.1130822_RLEN* s - &
          0.00846934_RLEN*s15 + dlog(ONE-0.001005_RLEN* s ))
  case (2)
     ! ---------------------------------------------------------------------
     ! Mehrbach et al. (1973) as refitted by Dickson and Millero (1987) 
     ! pH on seawater scale  (Millero, 1995, p.664)
     ! Standard OCMIP computation
     ! ---------------------------------------------------------------------
     k1 = 10.0_RLEN**(-ONE*(3670.7_RLEN*invtk - 62.008_RLEN + 9.7944_RLEN*dlogtk - &
          0.0118_RLEN * s + 0.000116_RLEN*s2))
     k2 = 10.0_RLEN**(-1.0*(1394.7_RLEN*invtk + 4.777_RLEN - &
          0.0184_RLEN*s + 0.000118_RLEN*s2))
  case (3)
     ! ---------------------------------------------------------------------
     ! Mehrbach et al (1973) refit by Lueker et al. (2000). pH total scale 
     ! ---------------------------------------------------------------------
     k1 = 10.0_RLEN**(-ONE*(3633.86_RLEN*invtk - 61.2172_RLEN + 9.6777_RLEN*dlogtk - &
          0.011555_RLEN * s + 0.0001152_RLEN * s2))
     k2 = 10.0_RLEN**(-ONE*(471.78_RLEN*invtk + 25.9290_RLEN - &
          3.16967_RLEN*dlogtk - 0.01781_RLEN * s + 0.0001122_RLEN * s2))
  case (4) 
     !-----------------------------------------------------------------------
     ! Hansson (1973b) data as refitted by Dickson and Millero (1987).
     ! pH seawater scale
     !-----------------------------------------------------------------------
     k1 = 10.0_RLEN**(-ONE*(851.4_RLEN*invtk + 3.237_RLEN - 0.0106_RLEN*s + 0.000132_RLEN*s2))
     k2 = 10.0_RLEN**(-ONE*(-3885.4_RLEN*invtk + 125.844_RLEN - 18.141_RLEN*dlogtk -  &
          0.0192_RLEN*s + 0.000132_RLEN* s2))    

  end select

   
  ! ---------------------------------------------------------------------
  ! k1p = [H][H2PO4]/[H3PO4] on hSWS
  !
  ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
  ! Millero p.670 (1995)
  ! ---------------------------------------------------------------------
  k1p = exp(-4576.752_RLEN*invtk + 115.525_RLEN - 18.453_RLEN * dlogtk + &
       (-106.736_RLEN*invtk + 0.69171_RLEN) * dsqrts +               &
       (-0.65643_RLEN*invtk - 0.01844_RLEN) * s)
   
  ! ---------------------------------------------------------------------
  ! k2p = [H][HPO4]/[H2PO4] 
  !
  ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
  ! Millero p.670 (1995)
  ! ---------------------------------------------------------------------
  if (CalcCO2ITER) then
     ! pH on Sea Water Scale
     phscale = 172.0883_RLEN
  else
     ! pH on total Scale
     phscale = 172.1033_RLEN
  end if
  k2p = exp(-8814.715_RLEN*invtk + phscale - 27.927_RLEN * dlogtk + &
       (-160.340_RLEN*invtk + 1.3566_RLEN) * dsqrts +                 &
       (0.37335_RLEN*invtk - 0.05778_RLEN) * s)
   
  !------------------------------------------------------------------------
  ! k3p = [H][PO4]/[HPO4] pH on Sea Water Scale
  !
  ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  ! ---------------------------------------------------------------------
  k3p = exp(-3070.75_RLEN*invtk - 18.126_RLEN + &
       (17.27039_RLEN*invtk + 2.81197_RLEN) *   &
       dsqrts + (-44.99486_RLEN*invtk - 0.09984_RLEN) * s)

  !------------------------------------------------------------------------
  ! ksi = [H][SiO(OH)3]/[Si(OH)4] pH on Sea Water Scale
  !
  ! Millero p.671 (1995) using data from Yao and Millero (1995)
  ! ---------------------------------------------------------------------
  ksi = exp(-8904.2_RLEN*invtk + 117.400_RLEN - 19.334_RLEN * dlogtk + &
       (-458.79_RLEN*invtk + 3.5913_RLEN) * dsqrtis +              &
       (188.74_RLEN*invtk - 1.5998_RLEN) * is +                   &
       (-12.1652_RLEN*invtk + 0.07871_RLEN) * is2 +               &
       log(ONE-0.001005_RLEN*s))
   
  !------------------------------------------------------------------------
  ! kw = [H][OH]  ion product of water; pH on Sea Water Scale
  !
  ! Millero p.670 (1995) using composite data recommended DOE (1994)
  ! ---------------------------------------------------------------------
  kw = exp(-13847.26_RLEN*invtk + 148.9802_RLEN - 23.6521_RLEN * dlogtk + &
       (118.67_RLEN*invtk - 5.977_RLEN + 1.0495_RLEN * dlogtk) *          &
       dsqrts - 0.01615_RLEN * s)
   
  !------------------------------------------------------------------------
  ! ks = [H][SO4]/[HSO4] pH on free scale
  ! (converted to hSWS in ta_iter)
  ! Dickson (1990, J. chem. Thermodynamics 22, 113)
  ! ---------------------------------------------------------------------
  ks = exp(-4276.1_RLEN*invtk + 141.328_RLEN - 23.093_RLEN*dlogtk +      &
       (-13856._RLEN*invtk + 324.57_RLEN - 47.986_RLEN*dlogtk) * dsqrtis + &
       (35474._RLEN*invtk - 771.54_RLEN + 114.723_RLEN*dlogtk) * is -     &
       2698._RLEN*invtk*is**1.5_RLEN + 1776._RLEN*invtk*is2 +              &
       dlog(ONE - 0.001005_RLEN*s))
   
  if (CalcCO2ITER) then
     !------------------------------------------------------------------------
     ! kf = [H][F]/[HF] pH on free scale
     ! (converted to hSWS in ta_iter)
     ! Dickson and Riley (1979)  also Dickson and Goyet (1994)
     ! ---------------------------------------------------------------------
     kf = exp(1590.2_RLEN*invtk - 12.641_RLEN + 1.525_RLEN*dsqrtis + &
          dlog(ONE - 0.001005_RLEN*s)) 
  else
     ! for Follows et al use pH total scale
     kf = exp(1590.2_RLEN*invtk - 12.641_RLEN + 1.525_RLEN*dsqrtis + &
          dlog(ONE - 0.001005_RLEN*s) + &
          dlog(ONE + (0.1400_RLEN/96.062_RLEN)*(scl)/ks))
  end if

   
  ! ---------------------------------------------------------------------
  ! kb = [H][BO2]/[HBO2] 
  !
  ! Millero p.669 (1995) using data from Dickson p.673 (1990)
  ! ---------------------------------------------------------------------
  if (CalcCO2ITER) then
  ! pH on Sea Water Scale
  kb = exp((-8966.90_RLEN - 2890.53_RLEN*dsqrts - 77.942_RLEN*s + &
       1.728_RLEN*s15 - 0.0996_RLEN*s2)*invtk +                  &
       (148.0248_RLEN + 137.1942_RLEN*dsqrts + 1.62142_RLEN*s) +  &
       (-24.4344_RLEN - 25.085_RLEN*dsqrts - 0.2474_RLEN*s) *     &
       dlogtk + 0.053105_RLEN*dsqrts*tk  +                        &
       dlog((ONE+(st/ks)+(ft/kf))/(ONE+(st/ks))) ) 
  else
  ! pH on Total Scale
  kb = exp((-8966.90_RLEN - 2890.53_RLEN*dsqrts - 77.942_RLEN*s + &
       1.728_RLEN*s15 - 0.0996_RLEN*s2)*invtk +                  &
       (148.0248_RLEN + 137.1942_RLEN*dsqrts + 1.62142_RLEN*s) +  &
       (-24.4344_RLEN - 25.085_RLEN*dsqrts - 0.2474_RLEN*s) *     &
       dlogtk + 0.053105_RLEN*dsqrts*tk)
  end if
   

  if (CalcCO2ITER) then
  !------------------------------------------------------------------------
  !
  ! Calculate [H+] at SWS when DIC and TA are known at T, S and 1 atm.
  ! The solution converges to err of XACC. The solution must be within
  ! the range x1 to x2.
  !
  ! If DIC and TA are known then either a root finding or iterative method
  ! must be used to calculate htotal. In this case we use the Newton-Raphson
  ! "safe" method taken from "Numerical Recipes" (function "rtsafe.f90" with
  ! error trapping removed).
  !
  ! As currently set, this procedure iterates about 12 times. The x1 and x2
  ! values set in mod_co2 will accomodate ANY oceanographic values. 
  ! The first time step is run with x1 and x2 set as default values. 
  ! After that, x1 and x2 are set to the previous value of the pH +/- 0.5. 
  ! The current setting of XACC will result in output accurate to 3 
  ! significant figures (xx.y). Making XACC bigger will result in faster 
  ! convergence also, but this is not recommended 
  ! (XACC of 10**-9 drops precision to 2 significant figures)
  ! ---------------------------------------------------------------------
  ! IMPORTANT NOTE:
  ! This part is optimized for scalar looping on a vector machine.
  ! May not perform well on single cache architectures.
  ! The input arrays are scanned and each element assigned to a scalar 
  ! which is defined in the module mem_CO2
  ! ---------------------------------------------------------------------
     do l = 1,NO_BOXES
        x1 = 10.0_RLEN**(-phhi(l))
        x2 = 10.0_RLEN**(-phlo(l))
        s_bt = bt(l)
        s_st = st(l)
        s_ft = ft(l)
        s_sit= sit(l)
        s_pt = pt(l)
        s_dic= dic(l)
        s_ta = ta(l)
        s_k1 = k1(l)
        s_k2 = k2(l)
        s_kw = kw(l)
        s_kb = kb(l)
        s_ks = ks(l)
        s_kf = kf(l)
        s_k1p= k1p(l)
        s_k2p= k2p(l)
        s_k3p= k3p(l)
        s_ksi= ksi(l)
        htotal(l) = drtsafe(x1,x2,XACC)
        ! define the new range of pH for next calculations
        phlo(l) = -dlog10(htotal(l))-0.5_RLEN
        phhi(l) = -dlog10(htotal(l))+0.5_RLEN
     end do
  else
  !---------------------------------------------------------------
  ! Approximate solution by Follows et al., 2006
  !---------------------------------------------------------------
#ifdef DEBUG
     LEVEL3 'htotal before',minval(htotal),maxval(htotal)
#endif
     xx=htotal
     xx2 = xx*xx
     xx3 = xx2*xx
     k11 = k1*k1
     k12 = k1*k2
     k12p = k1p*k2p
     k123p = k12p*k3p
     a = xx3 + k1p*xx2 + k12p*xx + k123p
     c = T1 + st/ks
     cag = - bt/(T1 + xx/kb) - &
          kw/xx + xx -         &
          pt*k12p*xx/a -       &
          T2*pt*k123p/a -      &
          sit/(T1 + xx/ksi) +  &
          st/(T1 + ks/(xx/c))+ &
          ft/(T1 + kf/xx)+     &
          pt*xx3/a +           &
          ta
     gamm = dic/cag
     dummy= (ONE-gamm)*(ONE-gamm)*k11 - 4.0_RLEN*k12*(ONE-2.0_RLEN*gamm)
     htotal = 0.5_RLEN*((gamm-ONE)*k1 + dsqrt(dummy))
#ifdef DEBUG
     LEVEL3 'talk',minval(ta),maxval(ta)
     LEVEL3 'dsqrt(dummy)',maxval(dsqrt(dummy))
     LEVEL3 'dic after',minval(dic),maxval(dic)
     LEVEL3 'cag',minval(cag),maxval(cag)
     LEVEL3 'gamm',minval(gamm),maxval(gamm)
     LEVEL3 'htotal after',minval(htotal),maxval(htotal)
#endif
  end if !CalcCO2ITER
  

  !---------------------------------------------------------------
  ! Calculate [CO2] as defined in DOE Methods Handbook 1994 Ver.2, 
  ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
  ! Compute other diagnostic variables (HCO3, CO3 and pH) and 
  ! pCO2, the CO2 partial pressure in the water
  !---------------------------------------------------------------
  htotal2 = htotal*htotal
  co2 = dic*htotal2 / (htotal2 + k1*htotal + k1*k2)
  hco3 = k1 * co2 / htotal
  co3  = k2 * hco3/ htotal
  pH = -dlog10(htotal)
  pco2 = co2 / K0

! check values (Zeebe and Wolf-Gladrow)
! T = 25 degC
! S = 35 
! pco2(atm)=365 uatm
! pH = 8.1
! DIC = 2.1 mmol/kg
! co2 = 10.4 umol/kg
! hco3 = 1818 umol/kg
! co3 = 272 umol/kg


  !---------------------------------------------------------------
  ! NOTE: transformation from mol/kg -----> umol/kg
  !       only for output purposes
  ! pco2 is converted from atm  -----> uatm
  !---------------------------------------------------------------
  dic  = dic * MEG
  ta   = ta * MEG
  co2  = co2 * MEG
  hco3 = hco3 * MEG
  co3  = co3 * MEG
  pco2 = pco2 * MEG
#ifdef DEBUG
  LEVEL3 'pco2',minval(pco2),maxval(pco2)
  LEVEL3 'pH',minval(pH),maxval(pH)  
  LEVEL3 'co2',minval(co2),maxval(co2)
  LEVEL3 'K0',minval(K0),maxval(K0)
  LEVEL3 'k1',minval(k1),maxval(k1)
  LEVEL3 'k2',minval(k2),maxval(k2)
#endif

  !---------------------------------------------------------------
  ! Prepares input for air-sea flux computation
  ! Output is in mmol/m2/d
  !---------------------------------------------------------------
  call co2flux(EWIND,ETW(SRFindices),ERHO(SRFindices),EPCO2air, &
               pco2(SRFindices),K0(SRFindices),                &
               NO_BOXES_XY,CO2airflux)
#ifdef DEBUG
  LEVEL3 'pco2air',EPCO2air
  LEVEL3 'co2airflux',CO2airflux
#endif

  !---------------------------------------------------------------
  ! flux is positive downward. 
  ! Conversion from mmolC/m2/d to mgC/m3/d.
  ! In the water, the flux is subtracted from
  ! (or added to) the diagonal element of O3c (i.e. infinite source)
  ! The fraction of ice-free water is also considered
  !---------------------------------------------------------------
  tmpflux(SRFindices) = (ONE-EIce)*CO2airflux * MW_C / Depth(SRFindices)
  call flux_vector( iiPel, ppO3c,ppO3c, tmpflux )

#ifdef DEBUG
  write(*,*) 'tmpflux',minval(tmpflux),maxval(tmpflux) 
  write(*,"(4A11)") "dic","ta","pH","K Henry"
  write(*,"(4G12.6)") dic(1),ta(1),pH(1),K0(1)
  write(*,"(4A11)") "pco2","co2","co3","hco3"
  write(*,"(4G12.6)") pco2(1),co2(1),co3(1),hco3(1)
#endif

#ifdef PIPPO
! ======================================================================
!  Calculate CaCO3 saturation state parameters
! ======================================================================
!
! Omega for calcite and aragonite (< 1 = undersaturation; > 1 = 
! supersaturation)
! [derived from Zeebe & Wolf-Gladrow Matlab routine; but apparently
!  originally from Mucci, 1983 (?)]
!
!
!     conversion factor for mol/m3 --> mol/kg
      m3_to_kg  = 1. / (sigma0(1) * 1000.)
!
!     convert CO3 concentration to mol/kg       
      co3conc  = co3 * m3_to_kg     ! concentration in mol/kg
!
!     calcium ion concentration
!     [from : "Seawater : Its composition, properties and behaviour"
!      (2nd Edition), Open University Course Team, 1995]
!     seawater concentration    = 412 mg / l
!     atomic weight             = 40.078 g
!     therefore, concentration  = 10.279 mmol / l
!                               = 10.279 mol / m3
!
!     Ca concentration normalised to salinity
      calconc  = 10.279 * (sal(1) / 35.)
!     calconc  = 10.279            ! or not if one prefers
!
!     convert Ca concentration to mol/kg
      calconc  = calconc * m3_to_kg ! concentration in mol/kg
!
!     Ksp_calcite
!     tmp1 = -171.9065-0.077993.*T+2839.319./T+71.595.*dlog10(T);
!     tmp2 = +(-0.77712+0.0028426.*T+178.34./T).*dsqrt(S);
!     tmp3 = -0.07711.*S+0.0041249.*S.^1.5;
!     dlog10Kspc = tmp1 + tmp2 + tmp3;
!
      kspcal   = -171.9065 - (0.077993 * tk) + (2839.319 / tk) &
           + (71.595 * dlog10(tk)) + (dsqrt(sal(1)) * (-0.77712 +  &
           (0.0028426 * tk) + (178.34 / tk))) &
           - (0.07711 * sal(1)) + (0.0041249 * (sal(1)**1.5))
      kspcal2  = 10.0**kspcal      ! in mol/kg
!
!     Ksp_aragonite
!     tmp1 = -171.945-0.077993.*T+2903.293./T+71.595.*dlog10(T);
!     tmp2 = +(-0.068393+0.0017276.*T+88.135./T).*dsqrt(S);
!     tmp3 = -0.10018.*S+0.0059415.*S.^1.5;
!     dlog10Kspa = tmp1 + tmp2 + tmp3;
!
      ksparag  = -171.945 - (0.077993 * tk) + (2903.293 / tk)   &
           + (71.595 * dlog10(tk)) + (dsqrt(sal(1)) * (-0.068393 +   &
           (0.0017276 * tk) + (88.135 / tk)))   &
           - (0.10018 * sal(1)) + (0.0059415 * (sal(1)**1.5))   
      ksparag2 = 10.0**ksparag     ! in mol/kg
!
!     Omega_calcite
      omegacal = calconc * co3conc / kspcal2
!
!     Omega_aragonite
      omegarag = calconc * co3conc / ksparag2
#endif

  return
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end subroutine CO2Calc
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


