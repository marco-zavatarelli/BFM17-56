!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: mem_Phyto
!
! DESCRIPTION
!   Parameter values for the phytoplankton groups
!
!
! !INTERFACE
  module mem_Phyto
!
! !USES:

  use global_mem
  use mem,  ONLY: iiPhytoPlankton
  use bfm_error_msg
!  
!
! !AUTHORS
!   ERSEMII version by J.W. Baretta, H. Baretta-Bekker and W. Ebenhoeh
!   Additional parametrizations by P. Ruardij and M. Vichi 
!   Dynamical allocation by G. Mattia 
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij and M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)!
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
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Phyto PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  ! NAME         [UNIT]/KIND            DESCRIPTION
  !        :     --------- Physiological parameters -----------------
  !  p_q10       [-]            Characteristic Q10 coefficient
  !  p_qtemp     [-]            Cut-off threshold for temperature factor
  !  p_sum       [1/d]          Maximal productivity at 10 degrees C
  !  p_srs       [1/d]          Respiration rate at 10 degrees C
  !  p_sdmo      [1/d]          Max.specific nutrient-stress lysis rate
  !  p_thdo      [-]            Half saturation constant for nutrient stress lysis
  !  p_seo       [1/d]          Extra lysis rate (biomass density-dependent)
  !  p_sheo      [mg C/3]       Half saturation constant for extra lysis
  !  p_pu_ea     [-]            Excreted fraction of primary production
  !  p_pu_ra     [-]            Activity respiration fraction
  !  p_switchDOC [1-3]          Switch for the type of DOC excretion
  !                             This choice must be consistent with bacteria
  !                             1. All DOC is released as R1c (Vichi et al., 2007)
  !                             2. Activity DOC is released as R2c (Vichi et al., 2004)
  !                                (there is no nutrient-stress excretion)
  !                             3. All DOC is released as R2c (Polimene et al., 2006)
  !                             
  real(RLEN)  :: p_q10(iiPhytoPlankton)
  real(RLEN)  :: p_temp(iiPhytoPlankton)=ZERO
  real(RLEN)  :: p_sum(iiPhytoPlankton)
  real(RLEN)  :: p_srs(iiPhytoPlankton)
  real(RLEN)  :: p_sdmo(iiPhytoPlankton)
  real(RLEN)  :: p_thdo(iiPhytoPlankton)
  real(RLEN)  :: p_seo(iiPhytoPlankton)
  real(RLEN)  :: p_sheo(iiPhytoPlankton)
  real(RLEN)  :: p_pu_ea(iiPhytoPlankton)
  real(RLEN)  :: p_pu_ra(iiPhytoPlankton)
  integer(RLEN)  :: p_switchDOC(iiPhytoPlankton)
  !
  !              --------- Nutrient parameters in phytoplankton -----------------
  !  p_netgrowth [T or F]       Logical switch for nutrient-limited growth
  !                             .T. nutrient-balanced growth (Vichi et al.2004)
  !                             .F. nutrient-stress carbon excretion
  !                               (Baretta-Bekker et al.1995 and Vichi et al.2007)
  !  p_limnut    [1-3]          Switch for N-P co-limitation
  !                             0. Geometric mean
  !                             1. Threshold (Liebig-like)
  !                             2. Combined
  !                   ---- N limitation control ----
  !  p_qun       [m3/mgC/d]     Membrane affinity for N
  !  p_lN4       [mmolN/m3]     Half saturation constant for NH4 uptake preference over NO3
  !  p_qnlc      [mmolN/mgC]    Minimum quotum Si:C 
  !  p_qncPPY    [mmolN/mgC]    Reference quotum Si:C
  !  p_xqn       [-]            Multiplication factor for luxury storage
  !                   ---- P limitation control ----
  !  p_qup       [m3/mgC/d]     Membrane affinity for P
  !  p_qplc      [mmolP/mgC]    Minimum quotum Si:C 
  !  p_qpcPPY      [mmolP/mgC]    Reference quotum Si:C
  !  p_xqp       [-]            Multiplication factor for luxury storage
  !                   ---- Si limitation control ----
  !  p_switchSi  [1-2]          Switch for Silica limitation
  !                             1. Si limitation is controlled by external Si conc.
  !                             2. Si limitation is controlled by internal quota
  !  p_chPs      [mmolSi/m3]    Half saturation conc. for dissolved Si limitation
  !  p_Contois   [>=0]          If >0, use Contois formulation
  !  p_qus       [m3/mgC/d]     Membrane affinity for Si
  !  p_qslc      [mmolSi/mgC]   Minimum quotum Si:C 
  !  p_qscPPY      [mmolSi/mgC]   Reference quotum Si:C
  !                   ---- nutrient stressed sinking ----
  !  p_esNI      [-]            Nutrient stress threshold for sinking
  !  p_res       [m/d]          Maximum Sinking velocity (m/d)
  logical     :: p_netgrowth(iiPhytoPlankton)=.FALSE.
  integer     :: p_limnut(iiPhytoPlankton) 
  real(RLEN)  :: p_qun(iiPhytoPlankton)
  real(RLEN)  :: p_lN4(iiPhytoPlankton)
  real(RLEN)  :: p_qnlc(iiPhytoPlankton)
  real(RLEN)  :: p_qncPPY(iiPhytoPlankton)
  real(RLEN)  :: p_xqn(iiPhytoPlankton)
  real(RLEN)  :: p_qup(iiPhytoPlankton)
  real(RLEN)  :: p_qplc(iiPhytoPlankton)
  real(RLEN)  :: p_qpcPPY(iiPhytoPlankton)
  real(RLEN)  :: p_xqp(iiPhytoPlankton)
  integer     :: p_switchSi(iiPhytoPlankton) 
  real(RLEN)  :: p_chPs(iiPhytoPlankton)
  real(RLEN)  :: p_Contois(iiPhytoplankton)
  real(RLEN)  :: p_qus(iiPhytoPlankton)
  real(RLEN)  :: p_qslc(iiPhytoPlankton)
  real(RLEN)  :: p_qscPPY(iiPhytoPlankton)
  real(RLEN)  :: p_esNI(iiPhytoPlankton)
  real(RLEN)  :: p_res(iiPhytoPlankton)
  !
  !              --------- Chlorophyll parameters -----------
  !  p_switchChl [1-4]          Switch for Chla-a synthesis
  !  p_sdchl     [1/d]          Specific turnover rate for Chla 
  !  p_alpha_chl [mgC s m2/     Initial slope of the P-E curve
  !               mgChl/uE] 
  !  p_qlcPPY    [mgChla/mgC]   Reference quotum Chla:C 
  !  p_epsChla   [m2/mgChla]    Chla-specific extinction coefficient
  !  p_tochl_relt  [1/d]        Relaxation rate towards maximum Chla:C
  !  p_EpEk_or   [-]            Optimal value of E_PAR/E_K
  integer     :: p_switchChl(iiPhytoPlankton) 
  real(RLEN)  :: p_sdchl(iiPhytoPlankton)
  real(RLEN)  :: p_alpha_chl(iiPhytoPlankton)
  real(RLEN)  :: p_qlcPPY(iiPhytoPlankton)
  real(RLEN)  :: p_epsChla(iiPhytoPlankton)
  real(RLEN)  :: p_tochl_relt(iiPhytoplankton)
  real(RLEN)  :: p_EpEk_or(iiPhytoplankton)
  !
  !              --------- Light Adaptation parameters -----------
  !  p_iswLtyp   [0-6]    Shape of the productivity function ChlDynamicsFlag=1
  !                         0 : Steele (old ERSEM)  y*exp(1-y)
  !                         1 : Steele (Simpson)    y*exp(1-y)
  !                         2 : Ebenhoeh            2y/(1+y^2)
  !                         3 : ramp                min(1,y)
  !                         4 : step                1 if y>1 , 0 elsewhere
  !                         5 : Smith_average
  !                         6 : Smith II (actual_Irr)
  !  p_chELiPPY  [W/m2]   Maximum Iopt
  !  p_clELiPPY  [W/m2]   Minimum Iopt
  !  p_ruELiPPY  [1/d]    Maximum daily shift in Iopt (1/d)
  !  p_addepth   [m]      Adaptation depth. Meaningless with high-res models
  integer     :: p_iswLtyp(iiPhytoPlankton)
  real(RLEN)  :: p_chELiPPY(iiPhytoPlankton)
  real(RLEN)  :: p_clELiPPY(iiPhytoPlankton)
  real(RLEN)  :: p_ruELiPPY(iiPhytoPlankton)
  real(RLEN)  :: p_addepth(iiPhytoPlankton)
  !
  !              --------- Sinking parameters -----------
  !  p_rPIm      [m/d]    Phytoplankton background sinking rate
  real(RLEN)  :: p_rPIm(4)
#ifdef INCLUDE_PELFE
  !
  !              --------- Iron parameters -----------
  !  p_quf       [m3/mgC/d]     Membrane affinity for Fe
  !  p_qflc      [umolFe/mgC]   Minimum quotum Fe:C derived from 3 umol Fe/mol C
  !                             Sunda & Huntsman (1997), Nature, 390, p 389-392
  !  p_qfcPPY      [umolFe/mgC]   Reference quotum Fe:C 
  !  p_xqf       [-]            Multiplication factor for luxury storage
  real(RLEN)  :: p_quf(iiPhytoPlankton)
  real(RLEN)  :: p_qflc(iiPhytoPlankton)
  real(RLEN)  :: p_qfcPPY(iiPhytoPlankton)
  real(RLEN)  :: p_xqf(iiPhytoPlankton)
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitPhyto

  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPhyto()
  integer :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Phyto_parameters/ p_q10, p_sum, p_srs, p_sdmo, p_seo, p_pu_ea, &
                              p_temp, p_netgrowth,p_limnut, &
                              p_pu_ra, p_qnlc, p_qplc, p_qslc, &
                              p_qncPPY, p_qpcPPY, p_qscPPY, p_qlcPPY, &
                              p_qun, p_qup, p_qus, &
                              p_xqn, p_xqp, p_sheo, &
                              p_esNI, p_thdo, p_res, p_lN4, p_chPs, &
                              p_Contois, p_EpEk_or, p_tochl_relt,   &
                              p_switchDOC,p_switchSi,p_switchChl,  &
                              p_alpha_chl, p_sdchl, p_epsChla, p_iswLtyp, &
                              p_addepth, p_chELiPPY, p_clELiPPY, &
                              p_ruELiPPY, p_rPIm

#ifdef INCLUDE_PELFE
  namelist /Phyto_parameters_iron/ p_qflc, p_qfcPPY, p_xqf, p_quf
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading Phyto parameters.."
  open(NMLUNIT,file='Pelagic_Ecology.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=Phyto_parameters,err=101)
#ifdef INCLUDE_PELFE
  read(NMLUNIT,nml=Phyto_parameters_iron,err=101)
#endif
  close(NMLUNIT)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Check consistency of parameters according to the parametrization
  ! (loop over the phytoplankton groups and stop
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do i=1,iiPhytoPlankton
     write(LOGUNIT,*) "#  Checking Phyto parameters for group:",i
     select case ( p_switchSi(i) )
       case ( 1 ) ! external limitation
          if (p_chps(i)==ZERO)  &
          call bfm_error("Pelagic_Ecology.nml","Phyto_parameters p_switchS1=1: "//&
                         "External silica control is "//&
                         "selected but half saturation constant p_chps=0")
       case ( 2 ) ! internal limitation
          if (p_qus(i)==ZERO)  &
          call bfm_error("Pelagic_Ecology.nml","Phyto_parameters p_switchS1=2: "//&
                         "Internal silica control is "//&
                         "selected but membrane affinity p_qus=0")
     end select
     if (p_netgrowth(i)) then
        p_switchDOC(i) = 2
        write(LOGUNIT,*) "#  Balanced growth is activated: p_netgrowth=",p_netgrowth(i)
        write(LOGUNIT,*) "#  forcing p_switchDOC = 2"
     else if (p_switchDOC(i)==2) then
        write(LOGUNIT,*) "#  Balanced growth is not activated: p_netgrowth=",p_netgrowth(i)
        write(LOGUNIT,*) "#  do you really want p_switchDOC = 2?"
     end if
     write(LOGUNIT,*) "#  OK"
  end do

  write(LOGUNIT,*) "#  Namelist is:"
  write(LOGUNIT,nml=Phyto_parameters)
#ifdef INCLUDE_PELFE
  write(LOGUNIT,nml=Phyto_parameters_iron)
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPhyto.f90","Pelagic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitPhyto.f90","Phyto_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPhyto
  end module mem_Phyto
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
