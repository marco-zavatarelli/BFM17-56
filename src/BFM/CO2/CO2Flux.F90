#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!  BFM - Biogeochemical Flux Model 
!
! SUBROUTINE
!
! FILE
!
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
!
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
  subroutine CO2Flux(wind,temp,rho,pco2air,pco2sea,k0,ncomp,flux)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
use constants, ONLY: RLEN,HOURS_PER_DAY

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
IMPLICIT NONE

    integer   ,intent(IN) :: ncomp
    real(RLEN),intent(IN) :: wind(ncomp)  ! m/s
    real(RLEN),intent(IN) :: temp(ncomp)  ! deg C
    real(RLEN),intent(IN) :: rho(ncomp)   ! kg/m3
    real(RLEN),intent(IN) :: pco2air(ncomp) ! uatm
    real(RLEN),intent(IN) :: pco2sea(ncomp) ! uatm
    real(RLEN),intent(IN) :: k0(ncomp)    ! 1.e-6 mol / (l * 1.e-6 atm)
    real(RLEN),intent(OUT):: flux(ncomp)

    real(RLEN),parameter  :: C1=2073.1_RLEN
    real(RLEN),parameter  :: C2=125.62_RLEN
    real(RLEN),parameter  :: C3=3.6276_RLEN
    real(RLEN),parameter  :: C4=0.043219_RLEN
    real(RLEN),parameter  :: CO2SCHMIDT=660._RLEN
    real(RLEN),parameter  :: CM2M=0.01_RLEN

    real(RLEN),dimension(size(wind,1)) :: pschmidt,reacon,temp2, &
                                          k660,kex

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BEGIN compute
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
    !k660 = 2.5_RLEN*(0.5246_RLEN + 1.6256D-02*temp + 4.9946D-04*temp2) 

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Calculate wind dependency
    ! including conversion cm/hr => m/day :
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !kex = (k660 + 0.3_RLEN*wind*wind)*sqrt(pschmidt/CO2SCHMIDT)* &
    !      CM2M*HOURS_PER_DAY
    !Alternative way 
     kex = (0.3_RLEN*wind*wind)*sqrt(pschmidt/CO2SCHMIDT)* &
          CM2M*HOURS_PER_DAY

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! flux co2 in mmol/m2/day   
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! (m * d-1) * uatm * (mol * kg-1 * atm-1) * (kg * m-3)
    !     d-1   1.e-6      mol   m-2
    !     umol m-2 d-1 / 1000 = mmol/m2/d
    flux = kex * (pco2air - pco2sea) * k0 * rho / 1000.0_RLEN

#ifdef DEBUG
     write(0,*) 'co2air-co2sea',pco2air - pco2sea
#endif

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  return

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end subroutine CO2Flux
