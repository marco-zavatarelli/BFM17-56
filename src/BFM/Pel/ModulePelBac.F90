!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelBac
!
! DESCRIPTION
!   Module containing the parameters for bacterioplankton and the 
!   initialization and consistency check
!
! !INTERFACE
  module mem_PelBac
!
! !USES:
  use global_mem
!  
!
! !AUTHORS
!   First ERSEM version by J.W. Baretta and H. Baretta-Bekker
!   Additional parametrizations by P. Ruardij and M. Vichi 
!   (Vichi et al., 2004; Vichi et al., 2007)
!   L. Polimene, I. Allen and M. Zavatarelli (Polimene et al., 2006)
!   Dynamical allocation by G. Mattia 
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! PelBac PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer     :: p_version  ! Switch for DOM uptake parameterization
  real(RLEN)  :: p_q10      ! Q10-value (temperature dependency)
  real(RLEN)  :: p_chdo     ! Michaelis const for O2 dependence (mmol/m3)
  real(RLEN)  :: p_sd       ! Independent specific mortality (1/d)
  real(RLEN)  :: p_sd2      ! Density dependent mortality (value: 0.009) (1/d)
  real(RLEN)  :: p_suhR1    ! Specific potential of rich DOM availability (1/d)
  real(RLEN)  :: p_sulR1    ! Specific potential sugar availability (1/d)
  real(RLEN)  :: p_suR2     ! Specific potential TEP availability (1/d)
  real(RLEN)  :: p_suR6     ! Availability of POM (1/d)
  real(RLEN)  :: p_suR7     ! Availability of semi-refractory DOC (1/d)
  real(RLEN)  :: p_sum      ! Specific potential uptake (1/d)
  real(RLEN)  :: p_pu_ra    ! Activity respiration (-)
  real(RLEN)  :: p_pu_ra_o  ! Decrease in Ass. efficiency at low O2 conc (-).
  real(RLEN)  :: p_srs      ! Specific rest respiration (1/day)
  real(RLEN)  :: p_qnc      ! Optimal N/C ratio (model units) 45:9:1
  real(RLEN)  :: p_qpc      ! Optimal P/C ratio (model units) C:N:P
  real(RLEN)  :: p_qlnc     ! Minimal N/C ratio (model units) 45:9:1 
  real(RLEN)  :: p_qlpc     ! Minimal P/C ratio (model units) C:N:P 
  real(RLEN)  :: p_qun      ! nutrient affinity ( mmol/mgC/day) 
  real(RLEN)  :: p_qup      ! nutrient affinity ( mmol/mgC/day) 
  real(RLEN)  :: p_lN4      ! ammonium conc. at which nitrate uptake are equal
  real(RLEN)  :: p_chn      ! half saturation ammonium conc. for uptake
  real(RLEN)  :: p_chp      ! half saturation phosphate conc. for uptake
  real(RLEN)  :: p_pu_ea_R7 ! excretion of semirefractory DOC (-) 

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitPelBac

  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Initialization of Pelagic Bacteria
  ! Definition of th enamelist with all the parameters
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPelBac()
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /PelBac_parameters/ p_version, p_q10, p_chdo, p_sd, p_sd2, p_suhR1, &
    p_sulR1, p_suR2, p_suR6, p_sum, p_pu_ra, p_pu_ra_o, p_pu_ea_R7, p_srs, &
    p_suR7, p_qpc, p_qlpc, p_qnc, p_qlnc, p_qun, p_qup, p_lN4, p_chn, p_chp

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading PelBac parameters.."
  open(NMLUNIT,file='Pelagic_Ecology.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=PelBac_parameters,err=101)
  close(NMLUNIT)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Check consistency of parameters according to the parametrization
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  write(LOGUNIT,*) "#  using p_version=",p_version
  select case ( p_version )
    case ( 1 ) ! Polimene et al. (2006)
      p_sulR1 = ZERO
      write(LOGUNIT,*) "#  forcing p_sulR1=0"
    case ( 2 ) ! Vichi et al. 2007
      p_sulR1 = ZERO
      p_suR2 = ZERO
      p_suR7 = ZERO
      write(LOGUNIT,*) "#  forcing p_sulR1,p_suR2,p_suR7=0"
    case ( 3 ) ! Vichi et al. 2004
      p_suR7 = ZERO
      write(LOGUNIT,*) "#  forcing p_suR7=0"
  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Write parameter list to the log
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  write(LOGUNIT,*) "#  Namelist is:"
  write(LOGUNIT,nml=PelBac_parameters)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPelBac.f90","Pelagic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitPelBac.f90","PelBac_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPelBac
  end module mem_PelBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
