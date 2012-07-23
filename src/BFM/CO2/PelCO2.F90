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
                 Depth, flux_vector, DIC, EPCO2air, Ac, DIC
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
#ifdef NECSX
  ! use a vectorized version of the code
  call CalcCO2System_vector(MethodCalcCO2)
#else
  ! use a scalar version of the code
  ! To use the Presure correction of CSYS here the pr_in=EPS value
  do BoxNumber=1,NO_BOXES
     ! convert DIC and alkalinity from model units to diagnostic output
     ! mg C/m3 --> umol/kg
     ! mmol eq/m3 --> umol/kg
     DIC(BoxNumber) = O3c(BoxNumber)/MW_C/ERHO(BoxNumber)*1000.0_RLEN
     Ac(BoxNumber) = O3h(BoxNumber)/ERHO(BoxNumber)*1000.0_RLEN
     error= CalcCO2System(MethodCalcCO2,ESW(BoxNumber),    &
              ETW(BoxNumber),ERHO(BoxNumber),  &
              N1p(BoxNumber),N5s(BoxNumber),Ac(BoxNumber),&
              CO2(BoxNumber),HCO3(BoxNumber),CO3(BoxNumber),pH(BoxNumber),&
              pr_in=EPR(BoxNumber), DIC_in=DIC(BoxNumber),pCO2_out=pCO2(BoxNumber),& 
              omegacal=OCalc(BoxNumber),omegarag=OArag(BoxNumber))
#ifdef DEBUG
            write(LOGUNIT,*) "in PelCO2:"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ERHO',ERHO(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ESW',ESW(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'N1p',N1p(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'N5s',N5s(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DIC',DIC(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Ac',Ac(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'OCalc',OCalc(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'OArag',OArag(BoxNumber)
            write(LOGUNIT,'(''layer:'',I4,'' pH='',G12.6)') BoxNumber,pH(BoxNumber)
#endif
     if ( error > 0 ) then
#ifdef DEBUG
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ERHO',ERHO(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ESW',ESW(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'EPR',EPR(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'N1p',N1p(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'N5s',N5s(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DIC',DIC(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'OCalc',OCalc(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'OArag',OArag(BoxNumber)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Ac',O3h(BoxNumber)
            write(LOGUNIT,'(''layer:'',I4,'' pH='',G12.6)') BoxNumber,pH(BoxNumber)
            call BFM_ERROR("PelCO2Dynamics","pH outside range 2-11")
#endif
            pH(BoxNUmber)=-1.0;
            CO2(BoxNUmber)=-1.0;
            HCO3(BoxNUmber)=-1.0;
            CO3(BoxNUmber)=-1.0;
            pCO2(BoxNUmber)=-1.0;
     endif
  end do
#endif

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

  !---------------------------------------------------------------
  ! flux is positive downward. 
  ! Conversion from mmolC/m2/d to mgC/m3/d.
  ! The fraction of ice-free water is also considered
  ! Boundary variable first assigned, then the source term is 
  ! added to the Source/Sink arrays if the Flag is TRUE
  ! In the water, the flux is subtracted from
  ! (or added to) the diagonal element of O3c (i.e. infinite source)
  !---------------------------------------------------------------
  jsurO3c(:) =  jsurO3c(:) + (ONE-EICE(:))*CO2airflux(:) * MW_C
#ifdef DEBUG
!  write(*,"(4A11)") "dic","ta","pH"
!  write(*,"(3G12.6)") DIC(1),Ac(1),pH(1)
!  write(*,"(4A11)") "pco2","co2","co3","hco3"
!  write(*,"(4G12.6)") pco2(1),co2(1),co3(1),hco3(1)
#endif

  end subroutine PelCO2Dynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
