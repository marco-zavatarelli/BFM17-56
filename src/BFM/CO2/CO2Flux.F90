#include "DEBUG.h"
#include "INCLUDE.h"

#ifdef INCLUDE_PELCO2
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!  BFM - Biogeochemical Flux Model 
!
! SUBROUTINE
!
! FILE
!
! /*
! DESCRIPTION
!   This function computes the CO2 flux at the air-sea interface
!   as forced by temperature, wind and chemical enhancement
!   (Schneider et al., 1999).
!   Uses equation from:
!   R. Wanninkhof (1992), Relationship between windspeed and gas
!   exchange over the oecean
!   J. GeoPhys. Res. 97, 7373-7382
!
!   notes: 
!   - K0 = co2/pco2 [1.e-6 mol / (l * 1.e-6 atm)]
!   - exchange coefficient: deltapCO2 * bt * K0 
!   [1.e-6atm * cm/hr * 1.e-6mol/(l * 1.e-6atm)] = [cm/hr * 1.e-6mol / l]
!   - Temperature in degrees C	
!   test parameter  (DIC=2133, AC=2260, pco2=341), O7.c = AC-2210,
!*/
! AUTHORS
!   H. Thomas (NIOZ) adapted from the OCMIP standard files
!   T. Lovato (CMCC) added atmospheric pCO2 and chemical enhancement
! 
! CHANGE_LOG
!
! COPYING
!   Copyright (C) 2013, the BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2004  H. Thomas and the ERSEM team (vichi@ingv.it)
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
  subroutine CO2Flux()
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  use constants, ONLY: RLEN,HOURS_PER_DAY,ZERO_KELVIN,MW_C
  use global_mem, ONLY: RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: EWIND,ETW,ESW,ERHO,EPCO2air,pCO2,  &
                 NO_BOXES,NO_BOXES_XY,CO2airflux,EICE
  use mem, ONLY: iiPel, ppO3c, jsurO3c, CO2airflux, &
                 Depth, flux_vector
#endif
   use mem_Param, ONLY:  AssignAirPelFluxesInBFMFlag
#ifdef BFM_GOTM
  use bio_var, ONLY: SRFindices
#else
  use api_bfm, ONLY: SRFindices
#endif

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
IMPLICIT NONE

    real(RLEN) :: k0(NO_BOXES_XY)    ! 1.e-6 mol / (l * 1.e-6 atm)
    real(RLEN) :: wind(NO_BOXES_XY)  ! m/s
    real(RLEN) :: ice(NO_BOXES_XY)   ! fraction
    real(RLEN) :: temp(NO_BOXES_XY)  ! deg C
    real(RLEN) :: salt(NO_BOXES_XY)  ! -
    real(RLEN) :: rho(NO_BOXES_XY)   ! kg/m3
    real(RLEN) :: pco2air(NO_BOXES_XY) ! uatm
    real(RLEN) :: pco2sea(NO_BOXES_XY) ! uatm
    real(RLEN) :: tmpflux(NO_BOXES)

    real(RLEN),parameter  :: C1=2073.1_RLEN
    real(RLEN),parameter  :: C2=125.62_RLEN
    real(RLEN),parameter  :: C3=3.6276_RLEN
    real(RLEN),parameter  :: C4=0.043219_RLEN
    real(RLEN),parameter  :: CO2SCHMIDT=660.0_RLEN
    real(RLEN),parameter  :: d=0.31_RLEN
    real(RLEN),parameter  :: CM2M=0.01_RLEN
    integer, save :: first=0
    real(RLEN),allocatable,save,dimension(:) :: pschmidt,temp2,bt, &
                                         ken,tk,tk100,tk1002,ScRatio
    integer :: AllocStatus
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BEGIN compute
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   if (first==0) then
      first=1
      allocate(pschmidt(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating pschmidt"
      allocate(temp2(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating temp2"
      allocate(bt(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating bt"
      allocate(ken(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating ken"
      allocate(tk(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating tk"
      allocate(tk100(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating tk100"
      allocate(tk1002(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating tk1002"
      allocate(ScRatio(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating ScRatio"
   end if

    temp = ETW(SRFindices)
    salt = ESW(SRFindices)
    wind = EWIND(:)
    ice = EICE(:)
    rho = ERHO(SRFindices)
    pco2air = EPCO2air(:)
    pco2sea = pCO2(SRFindices)
    tmpflux(:) = ZERO
    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate Schmidt number,
    ! ratio between the kinematic viscosity and the molecular 
    ! diffusivity of carbon dioxide.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    temp2 = temp*temp
    pschmidt = (C1 - C2*temp + C3*temp2 - C4*temp2*temp)

    ScRatio = CO2SCHMIDT/pschmidt
    !
    ! ScRatio is limited to 0 when T > 40 °C  
    WHERE(ScRatio .le. 0.0_RLEN); ScRatio=0.0_RLEN ; END WHERE    
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Compute Chemical enhancement the Temperature dependent
    ! gas transfer 
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    bt = 2.5_RLEN*(0.5246_RLEN + 1.6256E-02_RLEN*temp + 4.9946E-04_RLEN*temp2) 

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Calculate wind dependency + Chemical enhancement
    ! including conversion cm/hr => m/day :
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ken = (bt + d*wind*wind)*sqrt(ScRatio)* CM2M*HOURS_PER_DAY
    !
    !Alternative way, without enhancement
    !kun = (d*wind*wind)*sqrt(ScRatio)*  CM2M*HOURS_PER_DAY
    !
    ! ---------------------------------------------------------------------
    ! K0, solubility of co2 in the water (K Henry)
    ! from Weiss 1974; K0 = [co2]/pco2 [mol kg-1 atm-1]
    ! ---------------------------------------------------------------------
    tk = temp - ZERO_KELVIN
    tk100 = tk/100.0_RLEN
    tk1002 = tk100*tk100
    k0 = exp(93.4517_RLEN/tk100 - 60.2409_RLEN + 23.3585_RLEN * log(tk100) +   &
       salt * (.023517_RLEN - 0.023656_RLEN * tk100 + 0.0047036_RLEN * tk1002))

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! flux co2 in mmol/m2/day   
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! (m * d-1) * uatm * (mol * kg-1 * atm-1) * (kg * m-3)
    !     d-1   1.e-6      mol   m-2
    !     umol m-2 d-1 / 1000 = mmol/m2/d
    CO2airflux(:) = ken * (pco2air - pco2sea) * k0 * rho / 1000.0_RLEN

  !---------------------------------------------------------------
  ! flux is positive downward. 
  ! Conversion from mmolC/m2/d to mgC/m3/d.
  ! The fraction of ice-free water is also considered
  ! Boundary variable first assigned, then the source term is 
  ! added to the Source/Sink arrays if the Flag is TRUE
  ! In the water, the flux is subtracted from
  ! (or added to) the diagonal element of O3c (i.e. infinite source)
  !---------------------------------------------------------------
  jsurO3c(:) = jsurO3c(:) + (ONE-ice(:)) * CO2airflux(:) * MW_C
  tmpflux(SRFindices) = jsurO3c(:) / Depth(SRFindices)
  if ( AssignAirPelFluxesInBFMFlag) then
     call flux_vector( iiPel, ppO3c,ppO3c, tmpflux )
  end if
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  return

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end subroutine CO2Flux

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Compute the pCO2 in the air at sea level
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine CalcPCO2Air()
  !
  ! Covert the atmospheric CO2 conctration into pCO2 
  ! 
  ! VARIABLES:
  ! AtmCO2   CO2 air Mixing Ratio                    ppmv
  ! EPCO2air Partial Pressure of atmospheric CO2     uatm
  ! EAPR     Atmospheric Sea Level Pressure          hPa
  ! ETDP     Air Dew Point temperature               °
  !  
  ! NOTES:
  ! The p(H2O vapor) is computed wiht August-Roche-Magnus formula,
  ! coefficients for water (aw,bw, cw) from Lawrence(2005)
  ! and ice (ai,bi,ci) surfaces from WMO 2000 Technical Regulations. 
  ! These parameters are to be used with temperature in  °C.
  ! Alternatively it is possible to use the formulations of Buck (1996) or
  ! Goff (1957) for water and Goff and Gratch (1946) for ice.
  !
  use global_mem,  ONLY: ONE,ZERO,RLEN,LOGUNIT
#ifdef NOPOINTERS
  use mem
#else
  use mem,         ONLY: NO_BOXES_XY, EPCO2air
#endif
  use mem_CO2,     ONLY: AtmCO2, AtmSLP, AtmTDP, pCO2Method
  use constants,   ONLY: ZERO_KELVIN
  !
  ! LOCAL variables
  ! 
  IMPLICIT NONE

  real(RLEN),dimension(NO_BOXES_XY) :: EATD, EAPR, e
  real(RLEN),parameter              :: aw=17.625_RLEN, bw=243.04_RLEN, cw=6.1094_RLEN
  real(RLEN),parameter              :: ai=21.8745584_RLEN, bi=265.5_RLEN, ci=6.1078_RLEN
  real(RLEN)                        :: atm2pa=9.8692326671601E-06_RLEN

! BOC

  EAPR = AtmSLP%fnow

  select case (pCO2Method) 
    ! 
    ! Approximate pCO2, pCO2 = Mixing ratio * p(atm)
  case(1)
       EPCO2air = AtmCO2%fnow * (EAPR * 100.0_RLEN) * atm2pa
#ifdef DEBUG
    write(LOGUNIT,*)
    write(LOGUNIT,*) " Control on pCO2 calculation, with method : ", pCO2Method
    write(LOGUNIT,*) " Atm Press: ",EAPR(1)," CO2 ppm: ",AtmCO2%fnow(1)," pCO2 : ", EPCO2air(1)
    write(LOGUNIT,*)
#endif
    ! 
    ! pCO2 = Mixing ratio * (p(air) - p(water vapor))
  case(2)
       EATD = AtmTDP%fnow
       ! Convert Temperature at Dew Point to Celsius degrees
       if (EATD(1) > 200.0) EATD = EATD + ZERO_KELVIN 
       ! August-Roche-Magnus formula, with coefficients from Lawrence(2005)
       ! Input : ETDP and EAPR
       ! Partial pressure of water vapor 
       e = cw * exp((aw * EATD) / (bw + EATD))  
       ! Partial pressure of CO2 
       EPCO2air = AtmCO2%fnow * (EAPR - e) * 100.0_RLEN * atm2pa
#ifdef DEBUG
    write(LOGUNIT,*)
    write(LOGUNIT,*) " Control on pCO2 calculation, with method : ", pCO2Method
    write(LOGUNIT,*) " Atm Press: ",EAPR(1)," CO2 ppm: ",AtmCO2%fnow(1)," pCO2 : ", EPCO2air(1)
    write(LOGUNIT,*) " e (pH2Ovapor): ",e(1)," T dew point :", EATD(1)
#endif
  end select
!EOC
  end subroutine CalcPCO2Air

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
