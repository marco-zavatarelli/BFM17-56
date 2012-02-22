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
!   - exchange coefficient: deltapCO2 * k660 * K0 
!   [1.e-6atm * cm/hr * 1.e-6mol/(l * 1.e-6atm)] = [cm/hr * 1.e-6mol / l]
!   - Temperature in degrees C	
!   test parameter  (DIC=2133, AC=2260, pco2=341), O7.c = AC-2210,
!*/
! AUTHORS
!   H. Thomas (NIOZ) adapted from the OCMIP standard files
! 
! CHANGE_LOG
!
! COPYING
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
  use mem, ONLY: iiPel, ppO3c, D3STATE, jsurO3c, CO2airflux, &
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
    real(RLEN),parameter  :: CO2SCHMIDT=660._RLEN
    real(RLEN),parameter  :: CM2M=0.01_RLEN
    integer, save :: first=0
    real(RLEN),allocatable,save,dimension(:) :: pschmidt,reacon,temp2, &
                                         k660,kex,tk,tk100,tk1002
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BEGIN compute
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer :: AllocStatus, DeallocStatus
   if (first==0) then
      first=1
      allocate(pschmidt(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating pschmidt"
      allocate(reacon(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating reacon"
      allocate(temp2(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating temp2"
      allocate(k660(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating k660"
      allocate(kex(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating kex"
      allocate(tk(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating tk"
      allocate(tk100(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating tk100"
      allocate(tk1002(NO_BOXES_XY),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating tk1002"
   end if

    temp = ETW(SRFindices)
    salt = ESW(SRFindices)
    wind = EWIND(:)
    ice = EICE(:)
    rho = ERHO(SRFindices)
    pco2air = EPCO2air(:)
    pco2sea = pCO2(SRFindices)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate Schmidt number,
    ! ratio between the kinematic viscosity and the molecular 
    ! diffusivity of carbon dioxide.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    temp2 = temp*temp
    pschmidt = (C1 - C2*temp + C3*temp2 - C4*temp2*temp)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! gas transfer velocity at a Schmidt number of 660
    ! Temperature dependent
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    k660 = 2.5_RLEN*(0.5246_RLEN + 1.6256D-02*temp + 4.9946D-04*temp2) 

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Calculate wind dependency
    ! including conversion cm/hr => m/day :
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    kex = (k660 + 0.3_RLEN*wind*wind)*sqrt(CO2SCHMIDT/pschmidt)* &
          CM2M*HOURS_PER_DAY
    !Alternative way 
    !kex = (0.3_RLEN*wind*wind)*sqrt(CO2SCHMIDT/pschmidt)* &
    !      CM2M*HOURS_PER_DAY

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
    CO2airflux(:) = kex * (pco2air - pco2sea) * k0 * rho / 1000.0_RLEN

  !---------------------------------------------------------------
  ! flux is positive downward. 
  ! Conversion from mmolC/m2/d to mgC/m3/d.
  ! The fraction of ice-free water is also considered
  ! Boundary variable first assigned, then the source term is 
  ! added to the Source/Sink arrays if the Flag is TRUE
  ! In the water, the flux is subtracted from
  ! (or added to) the diagonal element of O3c (i.e. infinite source)
  !---------------------------------------------------------------
  jsurO3c =  (ONE-ice(:))*CO2airflux(:) * MW_C
  tmpflux(:) = ZERO
  tmpflux(SRFindices) = jsurO3c / Depth(SRFindices)
  if ( AssignAirPelFluxesInBFMFlag) then
     call flux_vector( iiPel, ppO3c,ppO3c, tmpflux )
  end if
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  return

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end subroutine CO2Flux
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
