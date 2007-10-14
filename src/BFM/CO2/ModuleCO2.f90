!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CO2
!
! DESCRIPTION
! Module for CO2 dynamics in Sea Water
!
! !INTERFACE
  module mem_CO2
!
! !USES:
  use global_mem
  use constants, ONLY: RLEN,ZERO
  use mem, ONLY: NO_BOXES, NO_BOXES_XY,EPCO2air

!  
!
! !AUTHORS
!   Lavinia Patara and Marcello Vichi
!
!
!
! !REVISION_HISTORY

!
! COPYING
!   Copyright (C) 2006  Marcello Vichi and the BFM team (vichi@ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details. 
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
IMPLICIT NONE

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Default all is public
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   public
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! SHARED GLOBAL VARIABLES (explicit public declaration)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! scalar variables for the iterative computation of the [H+]
   real(RLEN)                          :: s_bt,s_st,s_ft,s_sit, &
                                          s_pt,s_dic,s_ta
   real(RLEN)                          :: s_k1,s_k2,s_kw,s_kb,  &
                                          s_ks,s_kf,s_k1p,s_k2p,&
                                          s_k3p,s_ksi
   ! pH ranges
   real(RLEN),allocatable,dimension(:) :: phlo,phhi,htotal
   ! Initial Partial pressure in the air
   real(RLEN)                          :: pco2air=365.0_RLEN
   ! Choice of the acidity constants parameterization
   ! K1K2==1 Roy et al. (1993); DOE (1994); pH on total scale
   ! K1K2==2 Mehrbach et al (1973) refit by Dickson & Millero (1987)
   !         OCMIP STANDARD; pH on Sea Water Scale
   ! K1K2==3 Mehrbach et al (1973) refit by Lueker et al. (2000)
   !         pH on total scale
   ! K1K2==4 Hansson (1973b) data as refitted by Dickson and 
   !         Millero (1987); pH on Sea Water Scale
   integer                :: K1K2=2

   ! Choice of [H+] numerical computation
   ! CalcCO2ITER=.true. Standard OCMIP iteration
   ! CalcCO2ITER=.false. Follows et al., Ocean Modelling 2006
   logical                :: CalcCO2ITER=.false.

contains

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine InitCo2()
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    namelist /CO2_parameters/pco2air,K1K2,CalcCO2ITER  
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BEGIN compute
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  Open the namelist file(s)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading CO2 Chemistry parameters.."
   open(NMLUNIT,file='CO2.nml',status='old',action='read',err=100)
   read(NMLUNIT,nml=CO2_parameters,err=101)
   close(NMLUNIT)
   write(LOGUNIT,*) "#  Namelist is:"
   write(LOGUNIT,nml=CO2_parameters)

   ! assign initial atmospheric pCO2 and water pH
   EPCO2air(:) = pco2air
   allocate(htotal(NO_BOXES)); htotal(:) = 6.0E-9_RLEN
   allocate(phlo(NO_BOXES)); phlo(:) = 6.0_RLEN
   allocate(phhi(NO_BOXES)); phhi(:) = 9.0_RLEN
   
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!END compute
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  return

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Local Error Messages
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
100 call error_msg_prn(NML_OPEN,"ModuleCO2.f90","CO2.nml")
101 call error_msg_prn(NML_READ,"ModuleCO2.f90","CO2_parameters")

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end subroutine InitCO2

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end module mem_CO2
