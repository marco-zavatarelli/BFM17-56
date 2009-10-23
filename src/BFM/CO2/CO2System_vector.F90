#include "cppdefs.h"
#include "DEBUG.h"
#include "INCLUDE.h"
!/*
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !SUBROUTINE: CO2System_vector
!
! DESCRIPTION
!       Calculate co2 equilibrium in seawater starting from 
!       total Alkalinity and total co2 (dic) at
!       a given temperature, salinity and other species
!       This is a specialized vectorized versione for NEC-SX architecture
!
!       INPUTS (BFM variables)
!       dic (O3c) = total dissolved inorganic carbon (umol C/kg) 
!       ac  (O3h) = Alkalinity (umol eq/kg)
!       N1p     = inorganic phosphate (mmol/m3) 
!       N5s     = inorganic silicate (mmol/m3) 
!       ETW     = temperature (degrees C)
!       ESW     = salinity (PSU)
!       ERHO    = density (kg/m3)
!       EPco2air        = atmospheric co2 partial pressure (ppm) 
!
!       diagnostic OUTPUTS
!       co2  = co2*(aq) (umol/kg)
!       pco2 = oceanic pco2 (uatm)
!       hco3 = bicarbonate (umol/kg)
!       co3  = carbonate (umol/kg)
!       co2airflux = air-sea co2 flux (mmol/m^2/d)
!
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
  subroutine CalcCO2System_vector(mode)

! !USES:
  use global_mem, ONLY: RLEN,ONE,ZERO
  use constants, ONLY: ZERO_KELVIN,MW_C
  use mem
  use mem_co2
  implicit none
! !INPUT PARAMETERS:
  integer,intent(IN)               :: mode
! !LOCAL VARIABLES:
  integer,save :: first=0
  integer :: AllocStatus, DeallocStatus
  real(RLEN),allocatable,save,dimension(:)  ::  K0 ! solubility : [Co2]=k0Ac Pco2
  real(RLEN),allocatable,save,dimension(:)  ::  K1 ! carbonate equilibrium I
  real(RLEN),allocatable,save,dimension(:)  ::  K2 ! carbonate equilibrium II
  real(RLEN),allocatable,save,dimension(:)  ::  Kw ! water dissociation
  real(RLEN),allocatable,save,dimension(:)  ::  Kb ! constant for Boron equilibrium
  real(RLEN),allocatable,save,dimension(:)  ::  Ks ! constant for bisulphate equilibrium
  real(RLEN),allocatable,save,dimension(:)  ::  Kf ! constant for hidrogen fluoride equilibirum
  real(RLEN),allocatable,save,dimension(:,:)::  Kp ! constants for phosphate equilibirum
  real(RLEN),allocatable,save,dimension(:)  ::  Ksi ! constant for silicic acid eqilbrium
  real(RLEN),allocatable,save,dimension(:)  :: ldic  ! local dic in mol/kg 
  real(RLEN),allocatable,save,dimension(:)  :: lpco2 ! local pco2
  real(RLEN),allocatable,save,dimension(:)  :: scl   ! chlorinity
  real(RLEN),allocatable,save,dimension(:)  :: ta    ! total Alkalinity
  real(RLEN),allocatable,save,dimension(:)  :: lnK,bt,ft,st,pt,sit
  real(RLEN),parameter             :: MEG=1.D6,XACC=1.D-20,  &
                                      PERMIL=ONE/1000.0_RLEN,&
                                      PERMEG=ONE/MEG
  integer     :: i, l, error
  real(RLEN),allocatable,save,dimension(:)  :: Hplus, Hplus2
  real(RLEN),allocatable,save,dimension(:)  :: tmp1, tmp2
  real(RLEN)                      :: intercept
  real(RLEN),allocatable,save,dimension(:)  :: tk,tk100,tk1002,    &
                                     invtk,is,is2,       &
                                     dlogtk,dsqrtis,s,s2,&
                                     dsqrts,s15
  ! Local variables for Follows' parameterization
  real(RLEN),allocatable,save,dimension(:) :: xx,xx2,xx3,k11,k12,k12p,k123p,a,c
  real(RLEN),allocatable,save,dimension(:) :: cag,ldummy,gamm
!EOP
!-------------------------------------------------------------------------!
!BOC
  if (first==0) then 
     first=1
     allocate(K0(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating K0"
     allocate(K1(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating K1"
     allocate(Kw(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Kw"
     allocate(Kb(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Kb"
     allocate(Ks(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Ks"
     allocate(Kf(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Kf"
     allocate(Kp(3,NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Kp"
     allocate(Ksi(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Ksi"
     allocate(ldic(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ldic"
     allocate(lpco2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating lpco2"
     allocate(scl(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating scl"
     allocate(ta(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ta"
     allocate(lnK(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating lNK"
     allocate(bt(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating bt"
     allocate(st(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating st"
     allocate(pt(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pt"
     allocate(sit(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sit"
     allocate(Hplus(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Hplus"
     allocate(Hplus2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Hplus2"
     allocate(tmp1(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tmp1"
     allocate(tmp2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tmp2"
     allocate(tk(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tk"
     allocate(tk100(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tk100"
     allocate(tk1002(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tk1002"
     allocate(invtk(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating invtk"
     allocate(is(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating is"
     allocate(is2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating is2"
     allocate(dlogtk(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating dlogtk"
     allocate(dsqrtis(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating dsqrtis"
     allocate(s(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating s"
     allocate(s2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating s2"
     allocate(dsqrts(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating dsqrts"
     allocate(s15(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating s15"
     allocate(xx(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating xx"
     allocate(xx2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating xx2"
     allocate(xx3(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating xx3"
     allocate(k11(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating k11"
     allocate(k12(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating k12"
     allocate(k12p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating k12p"
     allocate(k123p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating k123p"
     allocate(a(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating a"
     allocate(c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating c"
     allocate(cag(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating cag"
     allocate(ldummy(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ldummy"
     allocate(gamm(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating gamm"
   end  if

  !---------------------------------------------------------------------
  ! convert DIC and alkalinity from BFM units to diagnostic output
  ! mg C/m3 --> umol/kg
  ! mmol eq/m3 --> umol/kg
  !---------------------------------------------------------------------
  DIC(:) = O3c(:)/MW_C/ERHO(:)*1000.0_RLEN
  Ac(:) = O3h(:)/ERHO(:)*1000.0_RLEN
  !---------------------------------------------------------------------
  ! Change units from the input of mmol/m^3 -> mol/kg:
  ! (1 mmol/m^3)  x (1 m^3/1024.5 kg) / 1000
  ! The ocean actual density /*ERHO*/ is used.
  ! Note: mol/kg are actually what the body of this routine uses 
  ! for all calculations.  
  !---------------------------------------------------------------------
  pt = N1p(:)/ERHO(:)*PERMIL
  sit = N5s(:)/ERHO(:)*PERMIL
  ! convert from umol/kg to standard dic units mol/kg 
  ldic = DIC(:)*PERMEG
  ! convert input variable Alkalinity from umol/kg to mol/kg 
  ta = Ac(:)*PERMEG

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
  tk = ETW(:) - ZERO_KELVIN
  tk100 = tk/100.0_RLEN
  tk1002 = tk100*tk100
  invtk = ONE/tk
  dlogtk = dlog(tk)
  ! salinity
  s  = ESW(:)
  s2 = ESW(:)*ESW(:)
  dsqrts = dsqrt(ESW(:))
  s15 = ESW(:)**1.5_RLEN
  ! chlorinity
  scl = ESW(:)/1.80655_RLEN  
  ! ionic strength
  is = 19.924_RLEN*ESW(:)/ (1000._RLEN-1.005_RLEN*ESW(:))
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
  ! ff = k0(1-ph2O)*correction term for non-ideality
  !
  ! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
  ! ---------------------------------------------------------------------
  !ff = exp(-162.8301_RLEN + 218.2968_RLEN/tk100  +      &
  !     90.9241_RLEN*dlog(tk100) - 1.47696_RLEN*tk1002 +  &
  !     s * (.025695_RLEN - .025225_RLEN*tk100 +         &
  !     0.0049867_RLEN*tk1002))

  ! ---------------------------------------------------------------------
  ! K0, solubility of co2 in the water (K Henry)
  ! from Weiss 1974; K0 = [co2]/pco2 [mol kg-1 atm-1]
  ! ---------------------------------------------------------------------
  K0 = exp(93.4517_RLEN/tk100 - 60.2409_RLEN + 23.3585_RLEN * dlog(tk100) +   &
       ESW(:) * (.023517_RLEN - 0.023656_RLEN * tk100 + 0.0047036_RLEN * tk1002))

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
          0.08468345_RLEN* s -0.00654208_RLEN * s15+dlog(ONE -0.001005_RLEN* s )
     K1 = dexp(lnK)
     lnK = -3351.6106_RLEN/tk -9.226508_RLEN-0.2005743_RLEN*dlogtk+        &
          (-23.9722_RLEN/tk - 0.10690177_RLEN)*dsqrts + 0.1130822_RLEN* s - &
          0.00846934_RLEN*s15 + dlog(ONE-0.001005_RLEN* s )
     K2 = dexp(lnK)
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
  Kp(1,:) = dexp(lnK)
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
  Kp(2,:) = dexp(lnK)
  !------------------------------------------------------------------------
  ! k3p = [H][PO4]/[HPO4] 
  ! ph scale: total
  !
  ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
  ! ---------------------------------------------------------------------
  lnK = -3070.75_RLEN*invtk - 18.126_RLEN + &
       (17.27039_RLEN*invtk + 2.81197_RLEN) *   &
       dsqrts + (-44.99486_RLEN*invtk - 0.09984_RLEN) * s
  Kp(3,:) = dexp(lnK)

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
       dlog(ONE-0.001005_RLEN*s)
  Ksi = dexp(lnK)

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
  Kw = dexp(lnK)

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
       dlog(ONE - 0.001005_RLEN*s)
  Ks = dexp(lnK)

  !------------------------------------------------------------------------
  ! kf = [H][F]/[HF] ph on free scale
  ! ph scale: "free"
  !
  ! Dickson and Riley (1979)  also Dickson and Goyet (1994)
  ! ---------------------------------------------------------------------
  lnK = 1590.2_RLEN*invtk - 12.641_RLEN + 1.525_RLEN*dsqrtis + &
        dlog(ONE - 0.001005_RLEN*s)
  Kf = dexp(lnK)

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
  Kb = dexp(lnK)



  select case ( mode)
  case ( DYNAMIC )
        stop "BFM: Dynamic CO2 Not implemented for NECSX"
  case (FOLLOWS)
        !--------------------------------------------------
        ! Approximate solution by Follows et al., 2006
        ! First guess of H+: previous ph value
        !--------------------------------------------------
         xx= 10.0_RLEN**(-pH(:))
         xx2 = xx*xx
         xx3 = xx2*xx
         k11 = K1*K1
         k12 = K1*K2
         k12p = Kp(1,:)*Kp(2,:)
         k123p = k12p*Kp(3,:)
         a = xx3 + Kp(1,:)*xx2 + k12p*xx + k123p
         c = ONE + st/Ks
         cag = - bt/(ONE + xx/kb) - &
              kw/xx + xx -          &
              pt*k12p*xx/a -        &
              2.0_RLEN*pt*k123p/a -       &
              sit/(ONE + xx/ksi) +  &
              st/(ONE + ks/(xx/c))+ &
              ft/(ONE + kf/xx)+     &
              pt*xx3/a +            &
              ta
         gamm = ldic/cag
         ldummy= (ONE-gamm)*(ONE-gamm)*k11 - 4.0_RLEN*k12*(ONE-2.0_RLEN*gamm)
         Hplus = 0.5_RLEN*((gamm-ONE)*K1 + dsqrt(ldummy))
         !---------------------------------------------------------------
         ! Derive [co2] as defined in DOE Methods Handbook 1994 Ver.2, 
         ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
         ! Compute other diagnostic variables (hco3, co3 and ph) and 
         ! pco2, the co2 partial pressure in the water 
         !---------------------------------------------------------------
         Hplus2  =   Hplus*Hplus
         CO2(:)   = ldic*Hplus2/( Hplus2+ K1* Hplus+ K1* K2)
         pCO2(:)  =   CO2(:)/K0
         HCO3(:)  =   K1 * CO2(:) / Hplus
         CO3(:)   =   K2 * HCO3(:) / Hplus
         pH(:)    =  -dlog10(Hplus)

  case ( STATIC )
        !--------------------------------------------------
        ! Static approximate solution (for deep ocean)
        !--------------------------------------------------
        !-p/2(p) = tmp1:
        tmp1 = - 0.5_RLEN*( ta-ldic+ 4.0_RLEN* K2/ K1*( &
          2.0_RLEN* ldic- ta))/( 1.0_RLEN- 4.0_RLEN* K2/ K1)
        !q(p) = tmp2:
        tmp2 = - K2* (( 2.0_RLEN* ldic- ta))**(2.0_RLEN)/( &
          K1*( 1.0_RLEN- 4.0_RLEN* K2/ K1))
        !pco2:
        pCO2(:) = ( tmp1+ sqrt( tmp1* tmp1- tmp2))/ K0
        !co2:
        CO2(:)  =   K0 * pCO2(:)
        !ph:
        pH(:) = - dlog( K1* CO2(:)/(2.0_RLEN* ldic &
          - ta- 2.0_RLEN* CO2(:)))/ dlog( 10.0_RLEN)
        !co3:
        CO3(:)  =   ta - ldic + CO2(:)
        !hco3:
        HCO3(:)  =   ldic - CO2(:) - CO3(:)
  end select

  !---------------------------------------------------------------
  ! NOTE: transformation from mol/kg -----> umol/kg
  !       only for output purposes
  ! pco2 is converted from atm  -----> uatm
  !---------------------------------------------------------------
  CO2(:)  = CO2(:) * MEG
  HCO3(:) = HCO3(:) * MEG
  CO3(:)  = CO3(:) * MEG
  pCO2(:)  = pCO2(:) * MEG

end subroutine CalcCO2System_vector
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
