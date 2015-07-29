!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MicroZoo
!
! DESCRIPTION
!   Definition and reading of the parameters for the microzooplankton submodel
!
!
! !INTERFACE
  module mem_MicroZoo
!
! !USES:
  use global_mem
  use mem,  ONLY: iiMicroZooPlankton, iiPelBacteria, &
                  iiMesoZooPlankton, iiPhytoPlankton
!  
!
! !AUTHORS
!   Original: Hanneke Baretta-Bekker and Job Baretta
!   Additional parameterizations: P. Ruardij and M. Vichi
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-------------------------------------------------------------------------!
  !NAMELIST MicroZoo_parameters
  !-------------------------------------------------------------------------!
  !  MICRO-ZOOPLANKTON
  !
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
  ! p_qncMIZ     [mmolN/mgC]     Maximum quotum P:C
  ! p_qpcMIZ     [mmolN/mgC]     Maximum quotum N:C
  ! p_paPBA(z,b) [-]             Availability of pelagic Bacteria group b 
  !                              to Zooplankton group z
  ! p_paPPY(z,p) [-]             Availability of PhytoPlankton group p
  !                              to Zooplankton group z
  ! p_paMIZ(z,m) [-]             Availability of MicroZooplankton group m 
  !                              to Zooplankton group z
  !-------------------------------------------------------------------------!
  real(RLEN)  :: p_q10(iiMicroZooPlankton)  
  real(RLEN)  :: p_srs(iiMicroZooPlankton)
  real(RLEN)  :: p_sum(iiMicroZooPlankton)
  real(RLEN)  :: p_sdo(iiMicroZooPlankton)
  real(RLEN)  :: p_sd(iiMicroZooPlankton)
  real(RLEN)  :: p_pu(iiMicroZooPlankton)
  real(RLEN)  :: p_pu_ea(iiMicroZooPlankton)
  real(RLEN)  :: p_chro(iiMicroZooPlankton)
  real(RLEN)  :: p_chuc(iiMicroZooPlankton)
  real(RLEN)  :: p_minfood(iiMicroZooPlankton)
  real(RLEN)  :: p_qpcMIZ(iiMicroZooPlankton)
  real(RLEN)  :: p_qncMIZ(iiMicroZooPlankton)
  real(RLEN)  :: p_paPBA(iiMicroZooPlankton,iiPelBacteria)
  real(RLEN)  :: p_paPPY(iiMicroZooPlankton,iiPhytoPlankton)
  real(RLEN)  :: p_paMIZ(iiMicroZooPlankton,iiMicroZooPlankton)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitMicroZoo
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitMicroZoo()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /MicroZoo_parameters/ p_q10, p_srs, p_sum, p_sdo, p_sd, p_pu, &
    p_pu_ea, p_chro, p_chuc, p_minfood, p_qncMIZ, p_qpcMIZ,p_paPPY, p_paMIZ, p_paPBA
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading MicroZoo parameters.."
  open(NMLUNIT,file='Pelagic_Ecology.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=MicroZoo_parameters,err=101)
  close(NMLUNIT)
  write(LOGUNIT,*) "#  Namelist is:"
  write(LOGUNIT,nml=MicroZoo_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitMicroZoo.f90","Pelagic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitMicroZoo.f90","MicroZoo_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitMicroZoo
  end module mem_MicroZoo
!EOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
