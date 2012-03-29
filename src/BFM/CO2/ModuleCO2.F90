#ifdef INCLUDE_PELCO2
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
  use SystemForcing, only :ForcingName, ForcingField, FieldInit
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
   ! pco2air0  : initial constant value 
   ! AtmCO2    : structure of data from file time series
   real(RLEN)   :: pco2air0=365.0_RLEN ! uatm
   
   type(ForcingName)    :: ATMpCO2_N
   type(ForcingField)   :: ATMpCO2

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
   ! MethodCalcCO2=3 Follows et al., Ocean Modelling 2006
   !
   ! Parameters for MethodCalcCO2=2
   ! M2XACC    :  accuracy of the iterative scheme for OCMIP (default 1.E-10)
   ! M2PHDELT  :  delta of pH for the root search (realized pH+/-DELT)
   !              in the OCMIP scheme (default 0.5)
   ! M2MAXIT   :  maximum number of iterations for OCMIP (default 100 )
   !
   integer      :: MethodCalcCO2=2 
   real(RLEN)   :: M2XACC=1.E-20_RLEN
   real(RLEN)   :: M2PHDELT=0.5_RLEN
   integer      :: M2MAXIT=100

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

#ifdef NOPOINTERS
  use mem
#else
  use mem,          ONLY: EPCO2air,pH,NO_BOXES_XY
#endif
!  use SystemForcing, ONLY : FieldInit

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    namelist /CO2_parameters/ pco2air0,K1K2,MethodCalcCO2,phscale,phstart,  &
                              M2XACC,M2PHDELT,M2MAXIT,ATMpCO2_N
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !--------------------------------------------------------------------------
  ! Initialize the structured array that defines if a variable is initialized  
  ! with external data.
  !---------------------------------------------------------------------------
                        ! Read  !   File     ! Netcdf  !  Var   ! File    ! Input      !   Time   !
                        ! Input !   name     ! Logical !  name  ! RefTime ! Frequency  !  interp  !
    ATMpCO2_N = ForcingName( 0  , "dummy.nc" , .TRUE.  ,"dummy" , "dummy" ,  "dummy"   ,  .TRUE.  )

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
    
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set initial conditions
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ATMpCO2%init = ATMpCO2_N%init 

    ! assign initial atmospheric pCO2
    if (AtmpCO2%init == 0) then
       ! Use constant atmospheric pCO2 
       EPCO2air(:) = pco2air0
       write(*,*) 'Use constant atmCO2', EPCO2air
    else
       ! read external AtmpCO2 
       CALL FieldInit(ATMpCO2_N, ATMpCO2)
       EPCO2air(:) = ATMpCO2%fnow
       write(*,*) 'EPCO2air Initial input', EPCO2air
    endif

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
#endif
