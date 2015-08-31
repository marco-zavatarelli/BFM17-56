#include "cppdefs.h"

#if defined INCLUDE_PELCO2 || defined INCLUDE_BENCO2
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
  use SystemForcing, only :ForcingName, ForcingField, FieldInit, FieldClose
  IMPLICIT NONE
!  
!
! !AUTHORS
!   L. Patara, M. Vichi (INGV-CMCC)
!   P. Ruardij, H. Thomas (NIOZ)
!
! !REVISION_HISTORY
!   T. Lovato (CMCC) 2012
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
  !NAMELIST CO2_parameters
  !-------------------------------------------------------------------------!
  ! CARBONATE SYSYEM SETTING
  ! NAME           [UNIT]/KIND             DESCRIPTION
  ! AtmCO20        [ppmv]           Initial atmospheric concentration of CO2
  ! calcAtmpCO2    logical          Compute the partial pressure of Atmospheric CO2
  ! pCO2Method     integer          pCO2 computation method: 1=MixRatio*slp0, 2=Magnus formula
  ! phstart        [pH]             Initial pH value
  ! K1K2           integer          Switch for the acidity constants parameterization
  !                                 1 : Roy et al. (1993); DOE (1994); pH on total scale
  !                                 2 : Default. OCMIP STANDARD; pH on Sea Water Scale
  !                                     Mehrbach et al (1973) refit by Dickson & Millero (1987)
  !                                 3 : Mehrbach et al (1973) refit by Lueker et al. (2000)
  !                                     pH on total scale
  !                                 4 : Hansson (1973b) data as refitted by Dickson and 
  !                                     Millero (1987);  pH on Sea Water Scale 
  ! MethodCalcCO2  numeric          Switch for the choice of [H+] numerical computation 
  !                                 1 : Approximate static solution 
  !                                 2 : Default. Standard OCMIP iteration 
  !                                 3 : Follows et al., Ocean Modelling 2006 
  ! CalcBioAlkFlag logical          Compute biological processes corrections on total alkalinity
  !              ---------  Parameters for MethodCalcCO2=2 -----------
  ! M2XACC         real             Accuracy of the iterative scheme for OCMIP (default 1.E-10) 
  ! M2PHDELT       [pH]             Delta of pH for the root search (realized pH+/-DELT)
  !                                 in the OCMIP scheme (default 0.5)
  ! M2MAXIT        integer          Maximum number of iterations for OCMIP (default 100 )
  !              ----------------------------------------------------- 
  ! Caconc0        [mol/m3]         Calcium ion concentration 
  !                                 ["Seawater : Its composition, properties and behaviour" 
  !                                 (2nd Edition), Open University Course Team, 1995]
  !                                 Seawater concentration   = 412 mg / l 
  !                                                        -> atomic weight = 40.078 g / mol
  !                                 therefore, concentration = 10.279 mmol / l = 10.279 mol / m3
  ! Canorm          logical         Normalize Calcium ion concentration by sea water salinity
  !              ---------  EXTERNAL DATA INPUT STRUCTURES -----------
  ! AtmCO2_N       structure        Read external data for atmospheric CO2 values
  ! AtmSLP_N       structure        Read external data for atmospheric sea level pressure
  ! AtmTDP_N       structure        Read external data for atmospheric dew-point temperature
  ! Example of general input structure for the data structure:
  !          ! Read  !   File                               ! NetCDF  !  Var    !
  !          ! Input !   name                               ! Logical !  name   !
  !AtmCO2_N  =    0  , 'CMIP5_Historical_GHG_1765_2005.dat' , .FALSE.  , 'CO2'  ,
  !          !  RefTime          ! Input      !   Time   !
  !          !  yyyymmdd         ! Frequency  !  interp  !
  !           '1764-07-01 00:00' ,  'yearly'  ,  .TRUE.
  !
  ! Convention for Input reading : 0 = use constant value (default if struct is not initialized)
  !                               2 = read timeseries file ( e.g. CO2 mixing ratios)
  !                               4 = field from a coupled model (e.g. atmospheric SLP from OGCM)
  ! NOTE: The file "CMIP5_Historical_GHG_1765_2005.dat" is located in "$BFMDIR/tools" folder
  !-----------------------------------------------------------------------------------!
   real(RLEN)   :: AtmCO20=365.0_RLEN ! ppm 
   logical      :: calcAtmpCO2=.FALSE. 
   integer      :: pCO2Method=1
   type(ForcingName)    :: AtmCO2_N, AtmSLP_N, AtmTDP_N
   type(ForcingField)   :: AtmCO2, AtmSLP, AtmTDP
   integer      :: K1K2=2
   integer      :: MethodCalcCO2=2 
   real(RLEN)   :: M2XACC=1.E-20_RLEN
   real(RLEN)   :: M2PHDELT=0.5_RLEN
   integer      :: M2MAXIT=100
   real(RLEN)   :: phstart=8.0_RLEN ! [-]
   integer      :: phscale = SWS
   real(RLEN)   :: Caconc0 = 10.279E0
   logical      :: Canorm = .TRUE.
   logical      :: CalcBioAlkFlag = .FALSE.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitCO2, CloseCO2

  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitCO2()

#ifdef NOPOINTERS
  use mem
#else
  use mem,          ONLY:EPCO2air,pH
  use mem,          ONLY:NO_BOXES_XY
#endif
  use api_bfm, ONLY: bfm_init
  use mem_Param, ONLY: slp0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    namelist /CO2_parameters/ AtmCO20,calcAtmpCO2,pCO2Method,K1K2,MethodCalcCO2,     &
                              phscale,phstart,M2XACC,M2PHDELT,M2MAXIT,         &
                              Caconc0,Canorm,AtmCO2_N,AtmSLP_N,AtmTDP_N, &
                              CalcBioAlkFlag
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer            ::error=0
  !---------------------------------------------------------------------------
  ! Initialize the structured array that defines if a variable is initialized  
  ! with external data.
  !---------------------------------------------------------------------------
                        ! Read  !   File     ! Netcdf  !  Var   ! File    ! Input      !   Time   !
                        ! Input !   name     ! Logical !  name  ! RefTime ! Frequency  !  interp  !
    AtmCO2_N = ForcingName( 0  , "dummy.nc" , .TRUE.  ,"AtmCO2" , "dummy" ,  "dummy"   ,  .TRUE.  )
    AtmSLP_N = ForcingName( 0  , "dummy.nc" , .TRUE.  ,"AtmSLP" , "dummy" ,  "dummy"   ,  .TRUE.  )
    AtmTDP_N = ForcingName( 0  , 'dummy.nc' , .TRUE.  ,'AtmTDP' , 'dummy' ,  'dummy'   ,  .TRUE.  )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
    LEVEL1 ' '
    LEVEL1 '     INITIALIZE PELAGIC CARBONATE SYSTEM       ' 
    LEVEL1 ' '
    LEVEL2 'Namelist content:'
    open(NMLUNIT,file='Carbonate_Dynamics.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=CO2_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,nml=CO2_parameters)
    LEVEL1 ' '
 
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set initial conditions
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Atmospheric CO2 concentration
    AtmCO2%init = AtmCO2_N%init 
    if (AtmCO2%init == 0) then
       ! Use constant  
       CALL FieldInit(AtmCO2_N, AtmCO2)
       ! the following check is needed to avoid allocation of empty arrays
       ! with MPI and land domains
       if (NO_BOXES_XY > 0) then
          AtmCO2%fnow = AtmCO20
          write(LOGUNIT,*) 'Using constant atmospheric CO2 concentration:', AtmCO2%fnow(1)
          write(LOGUNIT,*) ' '
       end if
    elseif (AtmCO2%init == 2) then
       ! read external 0-D timeseries
       CALL FieldInit(AtmCO2_N, AtmCO2)
       ! the following check is needed to avoid allocation of empty arrays
       ! with MPI and land domains
       if (NO_BOXES_XY > 0) then
          write(LOGUNIT,*) 'Using variable atmospheric CO2 concentration. Initial value:', AtmCO2%fnow(1)
          write(LOGUNIT,*) ' '
       end if
    elseif (AtmCO2%init == 4) then
       ! CO2 concentration is provided by external model
       CALL FieldInit(AtmCO2_N, AtmCO2)
       ! the following check is needed to avoid allocation of empty arrays
       ! with MPI and land domains
       if (NO_BOXES_XY > 0) then
          AtmCO2%fnow = AtmCO20
          write(LOGUNIT,*) 'CO2 conc. provided by external model. Initialize with default uniform value', AtmCO2%fnow(1)
          write(LOGUNIT,*) ' '
       end if
    endif
    ! Rough approximation: pCO2 is assumed equal to the mixing ratio of CO2
    if (.not. calcAtmpCO2) EPCO2air = AtmCO2%fnow

    ! COMPUTATION OF pCO2
    ! Sea Level Pressure
    AtmSLP%init = AtmSLP_N%init
    AtmTDP%init = AtmTDP_N%init
    if (calcAtmpCO2) then
       if (AtmSLP%init == 0) then
         ! Use constant
          CALL FieldInit(AtmSLP_N, AtmSLP)
          ! the following check is needed to avoid allocation of empty arrays
          ! with MPI and land domains
          if (NO_BOXES_XY > 0) then
             AtmSLP%fnow = slp0
             write(LOGUNIT,*) 'Using constant atmospheric SLP (see slp0 in BFM_General.nml): ', AtmSLP%fnow(1)
             write(LOGUNIT,*) ' '
          end if
       else
         CALL FieldInit(AtmSLP_N, AtmSLP)
         if (AtmSLP%init .eq. 2 ) &
            write(LOGUNIT,*) 'BFM reads atmospheric SLP from file: ', AtmSLP_N%filename
         if (AtmSLP%init .eq. 4 ) &
            write(LOGUNIT,*) 'BFM receives atmospheric SLP from coupled model, using sbc forcing for O3h (TA)'
         write(LOGUNIT,*) ' '
       endif

       ! Atmospheric Dew Point Temperature
       if (pCO2Method == 2 .AND. AtmTDP%init .ne. 0 ) then
          CALL FieldInit(AtmTDP_N, AtmTDP)
         if (AtmTDP%init .eq. 2 ) &
            write(LOGUNIT,*) 'BFM reads Dew Point Temperature from file: ', AtmTDP_N%filename
         if (AtmTDP%init .eq. 4 ) &
            write(LOGUNIT,*) 'BFM receives Dew Point Temperature from coupled model, using sbc forcing for N6r (Red. Equival.)'
         write(LOGUNIT,*) ' '
       else
          pCO2Method = 1
          write(LOGUNIT,*) 'pCO2Method is forced to 1 because AtmTDP%init is set to 0.'
       endif
    endif

    ! Assign initial pH 
    if (bfm_init /= 1 ) pH(:) = phstart

    ! Check consistency of input parameters
    select case (K1K2)
    case (1) 
        phscale = TOTAL
    case (2) 
        phscale = SWS
    end select
    
    if (calcAtmpCO2) write(LOGUNIT,*) 'BFM computes pCO2 with method: ', pCO2Method
    if (AtmSLP%init == 4 ) write(LOGUNIT,*) 'SLP is provided by external model.' 
    if (AtmTDP%init == 4 ) write(LOGUNIT,*) 'TDP is provided by external model.' 
    LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
    LEVEL1 ' '

    FLUSH(LOGUNIT)
    return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"ModuleCO2.f90","Carbonate_Dynamics.nml")
101 call error_msg_prn(NML_READ,"ModuleCO2.f90","CO2_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitCO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine CloseCO2()
    implicit none
    
    if (AtmCO2%init == 2) then
       ! close external 0-D timeseries
       CALL FieldClose(AtmCO2_N, AtmCO2)
    end if
  end subroutine CloseCO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  end module mem_CO2

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#endif
