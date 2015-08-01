#include "cppdefs.h"
#include "DEBUG.h"
#include "INCLUDE.h"

#ifdef INCLUDE_PELCO2
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelCO2Dynamics
!
! DESCRIPTION
!   !
!
! !INTERFACE
  subroutine PelCO2Dynamics()
!
! !USES:

  use global_mem, ONLY:RLEN,ONE,ZERO
  use constants, ONLY:MW_C
  use mem_Param, ONLY:  AssignAirPelFluxesInBFMFlag
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: iiPel, O3h, O3c, D3STATE, jsurO3c, CO2airflux,    &
                 Depth, flux_vector, DIC, EPCO2air, ALK, DIC
  use mem, ONLY: ppO3h, ppO3c, NO_BOXES, NO_BOXES_XY, BoxNumber,   &
    N1p,N5s,CO2, HCO3, CO3, pCO2, pH, ETW, ESW, ERHO, EWIND, EICE, &
    OCalc, OArag, EPR
#endif
  use CO2System, ONLY: CalcCO2System,CalcK0
  use mem_CO2    
  use BFM_ERROR_MSG, ONLY: BFM_ERROR
#ifdef BFM_GOTM
  use bio_var, ONLY: SRFindices
#else
  use api_bfm, ONLY: SRFindices
#endif
  IMPLICIT NONE

!  
!
! !AUTHORS
!   M. Vichi, H. Thomas and P. Ruardij
!
! !REVISION_HISTORY

! !LOCAL VARIABLES:
  integer            ::error=0
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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
  ! To use the Pressure correction of CSYS here the pr_in=EPS value
  do BoxNumber=1,NO_BOXES
     ! convert DIC and alkalinity from model units to diagnostic output
     ! mg C/m3 --> umol/kg
     ! mmol eq/m3 --> umol/kg
     DIC(BoxNumber) = O3c(BoxNumber)/MW_C/ERHO(BoxNumber)*1000.0_RLEN
     ALK(BoxNumber) = O3h(BoxNumber)/ERHO(BoxNumber)*1000.0_RLEN
     error= CalcCO2System(MethodCalcCO2,ESW(BoxNumber),    &
              ETW(BoxNumber),ERHO(BoxNumber),  &
              N1p(BoxNumber),N5s(BoxNumber),ALK(BoxNumber),&
              CO2(BoxNumber),HCO3(BoxNumber),CO3(BoxNumber),pH(BoxNumber),&
              pr_in=EPR(BoxNumber), DIC_in=DIC(BoxNumber),pCO2_out=pCO2(BoxNumber),& 
              omegacal=OCalc(BoxNumber),omegarag=OArag(BoxNumber))
#ifdef DEBUG
     write(LOGUNIT,*) "in PelCO2:"
     write(LOGUNIT,'(A,'' ='',f12.6)') 'ERHO',ERHO(BoxNumber)
     write(LOGUNIT,'(A,'' ='',f12.6)') 'ESW',ESW(BoxNumber)
     write(LOGUNIT,'(A,'' ='',f12.6)') 'N1p',N1p(BoxNumber)
     write(LOGUNIT,'(A,'' ='',f12.6)') 'N5s',N5s(BoxNumber)
     write(LOGUNIT,'(A,'' ='',f12.6)') 'DIC',DIC(BoxNumber)
     write(LOGUNIT,'(A,'' ='',f12.6)') 'ALK',ALK(BoxNumber)
     write(LOGUNIT,'(A,'' ='',f12.6)') 'OCalc',OCalc(BoxNumber)
     write(LOGUNIT,'(A,'' ='',f12.6)') 'OArag',OArag(BoxNumber)
     write(LOGUNIT,'(''layer:'',I4,'' pH='',f12.6)') BoxNumber,pH(BoxNumber)
#endif
     if ( error > 0 ) then
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',f12.6)') 'ERHO',ERHO(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'ETW',ETW(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'ESW',ESW(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'EPR',EPR(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'N1p',N1p(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'N5s',N5s(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'DIC',DIC(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'OCalc',OCalc(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'OArag',OArag(BoxNumber)
            write(LOGUNIT,'(A,'' ='',f12.6)') 'ALK',O3h(BoxNumber)
            write(LOGUNIT,'(''layer:'',I4,'' pH='',f12.6)') BoxNumber,pH(BoxNumber)
            call BFM_ERROR("PelCO2Dynamics","pH outside range 2-11")
     endif
  end do

  ! Rough approximation: pCO2 is assumed equal to the mixing ratio of CO2
  if (.not. calcAtmpCO2) EPCO2air = AtmCO2%fnow

  !---------------------------------------------------------------
  ! Computes Atmospheric pCO2 value
  !---------------------------------------------------------------
  if (calcAtmpCO2) call CalcPCO2Air()

  !---------------------------------------------------------------
  ! Computes air-sea flux (only at surface points)
  !---------------------------------------------------------------
  call CO2Flux()

  end subroutine PelCO2Dynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
