!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: SeaiceBac
!
! DESCRIPTION
!   Module for seaice bacteria.
!
!
! !INTERFACE
  module mem_SeaiceBac
!
! !USES:

  use global_mem
  use mem,  ONLY: iiSeaiceBacteria, iiT1
!  
!
! !AUTHORS
!   Original version by M. Vichi and L. Tedesco
!
!
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2009 M. Vichi and L. Tedesco
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
  ! SeaiceBac PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: p_version  ! Switch for DOM uptake parameterization
  !  p_version=1 <LUCA> Polimenes version
  !  p_version=2 <BFM> version
  real(RLEN)  :: p_q10(iiSeaiceBacteria)  ! Q10-value (temperature dependency)
  real(RLEN)  :: p_chdo(iiSeaiceBacteria)  ! Michaelis const for O2 dependence (mmol/m3)
  real(RLEN)  :: p_sd(iiSeaiceBacteria)  ! Independent specific mortality (1/d)
  real(RLEN)  :: p_sd2(iiSeaiceBacteria)  ! Density dependent mortality (value: 0.009) (1/d)
  real(RLEN)  :: p_suhU1(iiSeaiceBacteria)  ! Specific potential of rich DOM availability (1/d)
  real(RLEN)  :: p_sulU1(iiSeaiceBacteria)  ! Specific potential sugar availability (1/d)
  real(RLEN)  :: p_suU6(iiSeaiceBacteria)  ! Availability of POM (1/d)
  real(RLEN)  :: p_sum(iiSeaiceBacteria)  ! Specific potential uptake (1/d)
  real(RLEN)  :: p_pu_ra(iiSeaiceBacteria)  ! Activity respiration (-)
  real(RLEN)  :: p_pu_ra_o(iiSeaiceBacteria)  ! Decrease in Ass. efficiency at low O2 conc (-).
  real(RLEN)  :: p_srs(iiSeaiceBacteria)  ! Specific rest respiration (1/day)
  real(RLEN)  :: p_qnc(iiSeaiceBacteria)  ! Optimal N/C ratio (model units) 45:9:1
  real(RLEN)  :: p_qpc(iiSeaiceBacteria)  ! Optimal P/C ratio (model units) C:N:P
  real(RLEN)  :: p_qlnc(iiSeaiceBacteria)  ! Minimal N/C ratio (model units) 45:9:1 <BFM>
  real(RLEN)  :: p_qlpc(iiSeaiceBacteria)  ! Minimal P/C ratio (model units) C:N:P (BFM>
  real(RLEN)  :: p_qun(iiSeaiceBacteria)  ! nutrient affinity ( mmol/mgC/day) <BFM>
  real(RLEN)  :: p_qup(iiSeaiceBacteria)  ! nutrient affinity ( mmol/mgC/day) <BFM>
  real(RLEN)  :: p_lI4(iiSeaiceBacteria)  ! ammonium conc. at which nutrate uptake are equal (BFM)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitSeaiceBac
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSeaiceBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /SeaiceBac_parameters/ p_version, p_q10, p_chdo, p_sd, p_sd2, p_suhU1, &
    p_sulU1, p_suU6, p_sum, p_pu_ra, p_pu_ra_o, p_srs, p_qpc, p_qlpc, &
    p_qnc, p_qlnc, p_qun, p_qup, p_lI4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading SeaiceBac parameters.."
  open(NMLUNIT,file='SeaiceBac.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=SeaiceBac_parameters,err=101)
  close(NMLUNIT)
  write(LOGUNIT,*) "#  Namelist is:"
  write(LOGUNIT,nml=SeaiceBac_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitSeaiceBac.f90","SeaiceBac.nml")
101 call error_msg_prn(NML_READ,"InitSeaiceBac.f90","SeaiceBac_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSeaiceBac

  end module mem_SeaiceBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
