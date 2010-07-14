!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModulePelCO2
!
! DESCRIPTION
! Module for CO2 dynamics in Sea Water
!
! !INTERFACE
  module mem_CO2
!
! !USES:
  use global_mem
  IMPLICIT NONE
!  
!
! !AUTHORS
!   L. Patara, M. Vichi (INGV-CMCC)
!   P. Ruardij, H. Thomas (NIOZ)
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2007 P. Ruardij and M. Vichi
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  integer,parameter   :: STATIC=1
  integer,parameter   :: DYNAMIC=2
  integer,parameter   :: FOLLOWS=3
  integer,parameter   :: TOTAL=1
  integer,parameter   :: SWS=2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! PelCO2 PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   ! Initial Partial pressure in the air
   real(RLEN)   :: pco2air=365.0_RLEN ! uatm

   ! Choice of the acidity constants parameterization
   ! K1K2==1 Roy et al. (1993); DOE (1994); pH on total scale
   ! K1K2==2 Default
   !         Mehrbach et al (1973) refit by Dickson & Millero (1987)
   !         OCMIP STANDARD; pH on Sea Water Scale
   ! K1K2==3 Mehrbach et al (1973) refit by Lueker et al. (2000)
   !         pH on total scale
   ! K1K2==4 Hansson (1973b) data as refitted by Dickson and 
   !         Millero (1987); pH on Sea Water Scale
   integer      :: K1K2=2

   ! Choice of [H+] numerical computation
   ! MethodCalcCO2=1 Approximate static solution
   ! MethodCalcCO2=2 Default. Standard OCMIP iteration
   ! MethodCalcCO2=2 Follows et al., Ocean Modelling 2006
   integer      :: MethodCalcCO2=2 

   ! Initial pH value (needed for Follows)
   real(RLEN)   :: phstart=8.0_RLEN ! [-]

   ! Choice of pH scale (forced automatically according to K1K2)
   ! phscale = 1: Total
   ! phscale = 2: SeaWater Scale (OCMIP default)
   integer      :: phscale = SWS

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitCO2

  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitCO2()

  use mem, ONLY: EPCO2air,pH
  use global_mem
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    namelist /CO2_parameters/ pco2air,K1K2,MethodCalcCO2,phscale,phstart
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading PelCO2 parameters.."
    open(NMLUNIT,file='CO2.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=CO2_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=CO2_parameters)

    ! assign initial atmospheric pCO2 
    EPCO2air(:) = pco2air
    ! assign initial pH 
    pH(:) = phstart

    ! check consistency of input parameters
    select case (K1K2)
    case (1) 
        phscale = TOTAL
    case (2) 
        phscale = SWS
    end select

    return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"ModuleCO2.f90","CO2.nml")
101 call error_msg_prn(NML_READ,"ModuleCO2.f90","CO2_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitCO2

  end module mem_CO2

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
