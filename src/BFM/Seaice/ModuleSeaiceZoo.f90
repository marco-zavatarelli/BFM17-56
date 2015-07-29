!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: SeaiceZoo
!
! DESCRIPTION
!
! !INTERFACE
  module mem_SeaiceZoo
!
! !USES:

  use global_mem
  use mem,  ONLY: iiSeaiceZoo, iiX1, iiSeaiceAlgae, iiSeaiceBacteria

!  
!
! !AUTHORS
!   Letizia Tedesco and Marcello Vichi (CMCC)
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
!   Copyright (C) 2007 the BFM team
!   (vichi@bo.ingv.it)
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
  !NAMELIST SeaiceZoo_parameters
  !-------------------------------------------------------------------------!
  ! NAME         [UNIT]/KIND           DESCRIPTION
  ! p_q10        [-]             Q10 value for physiological rates
  ! p_srs        [1/d]           Respiration rate at 10 degrees Celsius
  ! p_sum        [1/d]           Potential growth rate
  ! p_sdo        [1/d]           Mortality rate due to oxygen limitation
  ! p_sd         [1/d]           Temperature independent mortality rate
  ! p_pu         [-]             Assimilation efficiency
  ! p_pu_ea      [-]             Fraction of activity excretion
  ! p_chro       [mmolO2/m3]     Half-saturation oxygen concentration
  ! p_chuc       [mgC/m3]        Half-saturation Food concentration for Type II
  ! p_minfood    [mgC/m3]        Half-saturation food concentration for
  !                              preference factor
  ! p_qncSZO     [mmolN/mgC]     Maximum quotum P:C
  ! p_qpcSZO     [mmolN/mgC]     Maximum quotum N:C
  ! p_paSBA(z,b) [-]             Availability of sea ice Bacteria group b
  !                              to Zooplankton group z
  ! p_paSAL(z,p) [-]             Availability of Algae group p
  !                              to Zooplankton group z
  ! p_paSZO(z,m) [-]             Availability of Zooplankton group m
  !                              to Zooplankton group z
  !-------------------------------------------------------------------------!
  !
  real(RLEN)  :: p_q10(iiSeaiceZoo)
  real(RLEN)  :: p_srs(iiSeaiceZoo)
  real(RLEN)  :: p_sum(iiSeaiceZoo)
  real(RLEN)  :: p_sdo(iiSeaiceZoo)
  real(RLEN)  :: p_sd(iiSeaiceZoo)
  real(RLEN)  :: p_pu(iiSeaiceZoo)
  real(RLEN)  :: p_pu_ea(iiSeaiceZoo)
  real(RLEN)  :: p_chro(iiSeaiceZoo)
  real(RLEN)  :: p_chuc(iiSeaiceZoo) 
  real(RLEN)  :: p_minfood(iiSeaiceZoo)
  real(RLEN)  :: p_qncSZO(iiSeaiceZoo)
  real(RLEN)  :: p_qpcSZO(iiSeaiceZoo)
  real(RLEN)  :: p_paSBA(iiSeaiceZoo,iiSeaiceBacteria) 
  real(RLEN)  :: p_paSAL(iiSeaiceZoo,iiSeaiceAlgae)
  real(RLEN)  :: p_paSZO(iiSeaiceZoo,iiSeaiceZoo)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitSeaiceZoo
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSeaiceZoo()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /SeaiceZoo_parameters/ p_q10, p_srs, p_sum, p_sdo, p_sd, p_pu, &
    p_pu_ea, p_chro, p_chuc, p_minfood, p_paSAL, p_paSZO, p_paSBA, p_qpcSZO, p_qncSZO
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading SeaiceZoo parameters.."
  open(NMLUNIT,file='Seaice_Ecology.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=SeaiceZoo_parameters,err=101)
  close(NMLUNIT)
  write(LOGUNIT,*) "#  Namelist is:"
  write(LOGUNIT,nml=SeaiceZoo_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"ModuleSeaiceZoo.f90","Seaice_Ecology.nml")
101 call error_msg_prn(NML_READ,"ModuleSeaiceZoo.f90","SeaiceZoo_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSeaiceZoo
  end module mem_SeaiceZoo
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
