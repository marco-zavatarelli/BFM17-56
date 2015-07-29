!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Seaicealgae
!
! DESCRIPTION
!   Parameter values for the sea ice algae group
!
! !INTERFACE
  module mem_SeaiceAlgae
!
! !USES:

  use global_mem
  use mem,  ONLY: iiSeaiceAlgae, iiS1, iiS2
!
!
! !AUTHORS
!   Letizia Tedesco and Marcello Vichi
!
!
!
! !REVISION_HISTORY
!   !
!
!
! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Sea ice algae PARAMETERS (read from nml)
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
  !  p_sheo      [mg C/3]       Half saturation constant for extra lysis
  !  p_pu_ea     [-]            Excreted fraction of primary production
  !  p_pu_ra     [-]            Activity respiration fraction
  real(RLEN)  :: p_q10(iiSeaiceAlgae)
  real(RLEN)  :: p_sum(iiSeaiceAlgae)
  real(RLEN)  :: p_srs(iiSeaiceAlgae)
  real(RLEN)  :: p_sdmo(iiSeaiceAlgae)
  real(RLEN)  :: p_thdo(iiSeaiceAlgae)
  real(RLEN)  :: p_pu_ea(iiSeaiceAlgae)
  real(RLEN)  :: p_pu_ra(iiSeaiceAlgae)
  !
  !  ---------------- Nutrient parameters in sea ice algae -----------------
  !  p_limnut    [1-3]          Switch for N-P co-limitation
  !                             0. Geometric mean
  !                             1. Threshold (Liebig-like)
  !                             2. Combined
  !                   ---- N limitation control ----
  !  p_qun       [m3/mgC/d]     Membrane affinity for N
  !  p_lN4       [mmolN/m3]     Half saturation constant for NH4 uptake preference over NO3
  !  p_qnlc      [mmolN/mgC]    Minimum quotum Si:C
  !  p_qnRc      [mmolN/mgC]    Reference quotum Si:C
  !  p_xqn       [-]            Multiplication factor for luxury storage
  !                   ---- P limitation control ----
  !  p_qup       [m3/mgC/d]     Membrane affinity for P
  !  p_qplc      [mmolP/mgC]    Minimum quotum Si:C
  !  p_qpRc      [mmolP/mgC]    Reference quotum Si:C
  !  p_xqp       [-]            Multiplication factor for luxury storage
  !                   ---- Si limitation control ----
  !  p_chsSAL    [mmolSi/m3]    Half saturation conc. for dissolved Si limitation
  !  p_qscSAL    [mmolSi/mgC]   Reference quotum Si:C
  integer  :: p_limnut(iiSeaiceAlgae)
  real(RLEN)  :: p_qun(iiSeaiceAlgae)
  real(RLEN)  :: p_lN4(iiSeaiceAlgae)
  real(RLEN)  :: p_qnlc(iiSeaiceAlgae)
  real(RLEN)  :: p_qncSAL(iiSeaiceAlgae)
  real(RLEN)  :: p_xqn(iiSeaiceAlgae)
  real(RLEN)  :: p_qup(iiSeaiceAlgae)
  real(RLEN)  :: p_qplc(iiSeaiceAlgae)
  real(RLEN)  :: p_qpcSAL(iiSeaiceAlgae)
  real(RLEN)  :: p_xqp(iiSeaiceAlgae)
  real(RLEN)  :: p_chsSAL(iiSeaiceAlgae)
  real(RLEN)  :: p_qscSAL(iiSeaiceAlgae)
  !
  !  ------------- Chlorophyll parameters -----------
  !  p_alpha_chl [mgC s m2/     Initial slope of the P-E curve
  !               mgChl/uE]
  !  p_qlcSAL    [mgChla/mgC]   Maximum quotum Chla:C
  !  p_epsSAL    [m2/mgChla]    Chla-specific extinction coefficient
  real(RLEN)  :: p_alpha_chl(iiSeaiceAlgae)
  real(RLEN)  :: p_qlcSAL(iiSeaiceAlgae)
  real(RLEN)  :: p_epsSAL(iiSeaiceAlgae)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitSeaiceAlgae
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSeaiceAlgae()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Seaicealgae_parameters/ p_q10, p_sum, p_srs, p_sdmo, p_pu_ea, &
    p_pu_ra, p_qnlc, p_qplc, p_qncSAL, p_qpcSAL, p_qscSAL, p_qun, p_qup, &
    p_xqn, p_xqp, p_thdo, p_lN4, p_chsSAL, &
    p_limnut, p_alpha_chl, p_qlcSAL, p_epsSAL
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading Sea ice algae parameters.."
    open(NMLUNIT,file='Seaice_Ecology.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=Seaicealgae_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Seaicealgae_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"ModuleSeaiceAlgae.f90","Seaice_Ecology.nml")
101 call error_msg_prn(NML_READ,"ModuleSeaiceAlgae.f90","Seaicealgae_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSeaiceAlgae
  end module mem_SeaiceAlgae
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
