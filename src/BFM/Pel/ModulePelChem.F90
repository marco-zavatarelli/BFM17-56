!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelChem
!
! DESCRIPTION
!   This process descibes the additional dynamics of dissolved
!   compounds in the water. Parameterized processes are:
!       - nitrification
!       - denitrification
!       - reoxidation of reduction equivalents
!       - iron scavenging and remineralization (INCLUDE_PELFE)
!
! !INTERFACE
  module mem_PelChem
!
! !USES:
  use global_mem
!  
!
! !AUTHORS
!   Original version by P. Ruardij and M. Vichi
!
!
!
! !REVISION_HISTORY
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
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! PelChem PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !NAMELIST PelChem_parameters, PelChem_parameters_iron
  !-------------------------------------------------------------------------!
  !  Pelagic Chemistry parameters
  ! NAME        [UNIT]/KIND      DESCRIPTION
  ! p_q10N4N3   [-]              Q10 factor for nitrification/denit
  ! p_sN4N3     [1/d]            Specific nitrification rate at 10 degC
  ! p_clO2o     [mmolO2/m3]      Half-saturation O2 concentration for
  !                              nitrification and reoxidation
  ! p_rOS       [1/d]            Specific reoxidation rate of reduction
  !                              equivalents
  ! p_sN3O4n    [1/d]            Specific denitrification rate                                           
  ! p_clN6r     [mmolHS/m3]      Half-saturation concentration of
  !                              reduction equivalents for denitrification
  ! p_rPAo      [mmolO2/m3/d]    Reference anoxic mineralization rate
  ! p_q10R6N5   [-]              Q10 factor for biogenic silica
  ! p_sR6N5     [1/d]            Specific remineralization rate of
  !                              biogenic silica
  real(RLEN)  :: p_q10N4N3
  real(RLEN)  :: p_sN4N3
  real(RLEN)  :: p_clO2o
  real(RLEN)  :: p_rOS
  real(RLEN)  :: p_sN3O4n
  real(RLEN)  :: p_clN6r
  real(RLEN)  :: p_rPAo
  real(RLEN)  :: p_q10R6N5
  real(RLEN)  :: p_sR6N5
#ifdef INCLUDE_PELFE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !NAMELIST PelChem_parameters, PelChem_parameters_iron
  !-------------------------------------------------------------------------!
  !              --------- Iron parameters -----------
  ! p_q10R6N7   [-]              Q10 temperature dependence
  ! p_sR6N7     [1/d]            Specific remineralization rate of particulate
  ! p_sR1N7     [1/d]            Specific remineralization rate of dissolved
  ! p_scavN7f   [1/d]            Specific scavenging rate
  ! p_N7fsol    [umolFe/m3]      Solubility concentration
  real(RLEN)  :: p_q10R6N7   ! Q10 temperature dependence
  real(RLEN)  :: p_sR6N7     ! Specific remineralization rate (d-1)
  real(RLEN)  :: p_sR1N7     ! Specific remineralization rate from chelated iron (d-1)
  real(RLEN)  :: p_scavN7f   ! Specific scavenging rate (d-1)
  real(RLEN)  :: p_N7fsol    ! Solubility concentration (umol Fe/m3)
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitPelChem
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPelChem()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /PelChem_parameters/ p_sN4N3, p_q10N4N3, p_q10R6N5, p_rOS, p_clO2o, &
    p_clN6r, p_sN3O4n, p_rPAo, p_sR6N5
#ifdef INCLUDE_PELFE
  namelist /PelChem_parameters_iron/ p_q10R6N7, p_sR6N7, p_sR1N7, p_scavN7f, &
    p_N7fsol
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading PelChem parameters.."
  open(NMLUNIT,file='Pelagic_Environment.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=PelChem_parameters,err=101)
#ifdef INCLUDE_PELFE
  read(NMLUNIT,nml=PelChem_parameters_iron,err=102)
#endif
  close(NMLUNIT)
  write(LOGUNIT,*) "#  Namelist is:"
  write(LOGUNIT,nml=PelChem_parameters)
#ifdef INCLUDE_PELFE
  write(LOGUNIT,nml=PelChem_parameters_iron)
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPelChem.f90","Pelagic_Environment.nml")
101 call error_msg_prn(NML_READ,"InitPelChem.f90","PelChem_parameters")
102 call error_msg_prn(NML_READ,"InitPelChem.f90","PelChem_parameters_iron")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPelChem
  end module mem_PelChem
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
