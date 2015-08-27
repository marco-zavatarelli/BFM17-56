#include "cppdefs.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Param
!
! DESCRIPTION
!   List of global model parameters 
!   (global variables that can be changed during the model initialization
!
! !INTERFACE
  MODULE mem_Param

!
! !USES:

  USE global_mem
  USE constants
  USE mem

!  
!
! !AUTHORS
!   Piet Ruardij and Marcello Vichi
!
! !REVISION_HISTORY
!   --------
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  PUBLIC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Global Switches : turn on/off or choose model components
  ! NAME                          KIND    DESCRIPTION
  ! CalcPelagicFlag               logical Pelagic System
  ! CalcBenthicFlag               numeric Benthic system
  !                                       0 = No Benthic System
  !                                       The following are Not Yet Activated
  !                                       1 = Simple Benthic Return
  !                                       2 = Benthic organisms and intermediate
  !                                           complexity nutrient regeneration
  !                                       3 = Benthic organisms and full nutrient
  !                                           regeneration (early diagenesis)
  ! CalcTransportFlag             logical Compute Transport Term (when coupled
  !                                       with a OGCM)
  ! CalcConservationFlag          logical Mass Conservation Check
  ! CalcPhytoPlankton             logical Pelagic Phytoplankton (vector)
  ! CalcPelBacteria               logical Pelagic Bacteria (vector)
  ! CalcMesoZooPlankton           logical Mesozooplankton (vector)
  ! CalcMicroZooPlankton          logical Microzooplankton (vector)
  ! CalcPelChemistry              logical Pelagic Hydrochemical Processes
  ! AssignPelBenFluxesInBFMFlag   logical Benthic-pelagic fluxes are added to the
  !                                       time integration
  ! AssignAirPelFluxesInBFMFlag   logical Air-sea fluxes are added to the
  !                                       time integration
  ! ChlDynamicsFlag               numeric Choose the dynamics of Chl-a
  !                                       1 = diagnostic, optimal light property
  !                                           in phytoplankton
  !                                           (Ebenhoeh et al 1995, ERSEM-II)
  !                                       2 = state variable, constituent of
  !                                           phytoplankton
  ! check_fixed_quota             numeric Check whether zooplankton have fixed quota
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  logical   :: CalcPelagicFlag=.TRUE.  
  integer   :: CalcBenthicFlag=0      ! Switch for Benthic system
  logical   :: CalcSeaiceFlag=.TRUE.  ! Switch for Seaice system

  logical   :: CalcTransportFlag=.FALSE.  
  logical   :: CalcConservationFlag=.TRUE.  
  logical   :: CalcPelChemistry=.TRUE.
  logical   :: AssignPelBenFluxesInBFMFlag=.TRUE.
  logical   :: AssignAirPelFluxesInBFMFlag=.TRUE.
  integer   :: ChlDynamicsFlag=2
  integer   :: check_fixed_quota=0

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Global Parameters : used throughout the model and not related 
  !                     to a specific component
  ! NAME          UNIT          DESCRIPTION
  ! p_small      [-]           Smallest numeric value (the model "zero")
  ! slp0         [mbar]        Reference sea level pressure
  ! p_pe_R1c     [-]           Fractional content of C in cytoplasm 
  ! p_pe_R1n     [-]           Fractional content of N in cytoplasm
  ! p_pe_R1p     [-]           Fractional content of P in cytoplasm
  ! p_qro        [mmolHS-/     Stoichiometric coefficient for
  !               mmolO2]      anaerobic reactions
  ! p_qon_dentri [mmolO2/      Stoichiometric coefficient for 
  !               mmolN]       denitrification 
  ! p_qon_nitri  [mmolO2/      Stoichiometric coefficient for 
  !               mmolN]       nitrification (3/2)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic model parameters
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  real(RLEN)   :: &
      p_small=1.0E-20_RLEN,  &
      slp0=1013.25_RLEN
  real(RLEN)   :: &
      p_pe_R1c=0.60_RLEN  ,     &
      p_pe_R1n=0.72_RLEN  ,     &
      p_pe_R1p=0.832_RLEN  ,    &
      p_pe_R1s=0.06_RLEN  ,  &
      p_qro=0.5_RLEN,  &  
      p_qon_dentri=1.25_RLEN, &  
      p_qon_nitri=1.5_RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Benthic model parameters
  ! NAME          UNIT          DESCRIPTION
  ! p_sedlevels   [-]           Number of sigma levels for benthic nutrient
  ! p_poro0       [-]           Constant porosity for 0D and 1D runs
  ! p_InitSink    Logical       parameter to Initialize BenthicSInk var.
  ! p_q10diff     [-]           Temperature-dependency porewater diffusion
  ! p_clDxm       [m]           minimal value of D?.m for calculation of the alpha
  ! p_d_tot       [m]           Thickness of modelled benthic sediment layers
  ! p_clD1D2m     [m]           minimum distance between D1m and D2m
  ! p_d_tot_2     [m]           maximal Thickness of D2m
  ! p_sedsigma    [-]           Parameter for sigma level distribution
  ! p_poro        [-]           Sediment porosity
  ! p_p_ae        [-]           Adsorption coefficient
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! 0d-parameters 
  integer      :: p_sedlevels=20 
  real(RLEN)   :: &
      p_sedsigma=2.0_RLEN, &
      p_d_tot=0.30_RLEN , &
      p_poro0=0.4
  ! 1d-parameters
  real(RLEN),public,dimension(:),allocatable   ::  p_p_ae, p_poro      
#ifdef INCLUDE_BEN
      integer   :: calc_init_bennut_states
  real(RLEN)   :: &
      p_InitSink=100.0_RLEN,  &  
      p_q10diff=1.49_RLEN,  &  
      p_clDxm=0.001_RLEN, &  
      p_clD1D2m=0.01_RLEN,    &  
      p_d_tot_2=0.35_RLEN,  &
      p_qnQIc, &
      p_qpQIc, &
      p_qsQIc
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Seaice model parameters
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#ifdef INCLUDE_SEAICE
#endif


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitParam
  
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitParam()
  use mem
  use constants
  use global_mem, ONLY: bfm_lwp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Param_parameters/ p_small, p_qro, p_qon_dentri,                   &
    p_qon_nitri,                                                              &
    CalcPelagicFlag, CalcBenthicFlag,CalcSeaiceFlag, CalcTransportFlag,       &
    CalcConservationFlag,CalcPhytoPlankton,CalcMicroZooPlankton,              &
    CalcPelChemistry,CalcMesoZooPlankton, CalcPelBacteria,                    &
    AssignPelBenFluxesInBFMFlag, AssignAirPelFluxesInBFMFlag,                 &
    slp0,                                                                     &
    p_pe_R1c, p_pe_R1n, p_pe_R1p, p_pe_R1s,                                   &
    ChlDynamicsFlag, check_fixed_quota
#ifdef INCLUDE_BEN
  namelist /Param_parameters_ben/                                             &
    p_sedlevels, p_sedsigma,p_d_tot, p_poro0,                                 &
    CalcBenOrganisms,CalcBenBacteria,                                         &
    calc_init_bennut_states, p_qnQIc, p_qpQIc, p_qsQIc,                       &
    p_InitSink, p_d_tot_2, p_clD1D2m, p_clDxm, p_q10diff
#endif
#ifdef INCLUDE_SEAICE
  namelist /Param_parameters_ice/                                             &
    CalcSeaiceAlgae,CalcSeaiceZoo,CalcSeaiceBacteria
#endif
   integer :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   LEVEL1 "#  Reading BFM parameters .."
   open(NMLUNIT,file='BFM_General.nml',status='old',action='read',err=100)
   read(NMLUNIT,nml=Param_parameters,err=101)
#ifdef INCLUDE_BEN
   read(NMLUNIT,nml=Param_parameters_ben,err=102)
#endif
#ifdef INCLUDE_SEAICE
   read(NMLUNIT,nml=Param_parameters_ice,err=103)
#endif
   close(NMLUNIT)
   if (bfm_lwp) then 
    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"      
    write(LOGUNIT,*) "#  Reading Param parameters.. "
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Param_parameters)
#ifdef INCLUDE_BEN
    write(LOGUNIT,nml=Param_parameters_ben)
#endif
#ifdef INCLUDE_SEAICE
    write(LOGUNIT,nml=Param_parameters_ice)
#endif
   endif 
   ! These initializations are done here because some compilers do not
   ! allow the initialization of constants with intrinsic functions
   MIN_VAL_EXPFUN=log(DBL_MIN)  
   DAY_PER_SEC=ONE/SEC_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitParam.f90","BFM_General.nml")
101 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters")
102 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters_ben")
103 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters_ice")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitParam

  end module
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
