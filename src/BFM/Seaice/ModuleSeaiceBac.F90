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
  !NAMELIST SeaiceBacteria_parameters
  !-------------------------------------------------------------------------!
  !  SEAICE BACTERIA
  !
  ! NAME         [UNIT]/KIND            DESCRIPTION
  ! p_version   integer         Switch for bacteria parameterization
  !                              1 : Baretta-Bekker et al. 1995;
  !                                  Vichi et al., 2007
  ! p_q10                        Q10-value (temperature dependency)
  ! p_chdo      [mmol/m3]        Half-saturation constant for O2 limitation
  ! p_sd        [1/d]            Specific mortality rate
  ! p_sd2       [1/d]            Density dependent specific mortality rate
  ! p_suhU1     [1/d]            Specific potential uptake for nutrient-rich DOM
  ! p_sulU1     [1/d]            Specific potential uptake for nutrient-poor DOM
  ! p_suU6      [1/d]            Specific potential uptake for POM (1/d)
  ! p_sum       [1/d]            Potential specific growth rate
  ! p_pu_ra     [-]              Activity respiration fraction
  ! p_pu_ra_o   [-]              Additional respiration fraction at low O2 conc
  ! p_srs       [1/d]            Specific rest respiration
  ! p_qncSBA    [mmolN/mgC]      Optimal N/C ratio 
  ! p_qpcSBA    [mmolP/mgC]      Optimal P/C ratio 
  ! p_chn       [mmolN/m3]       Half saturation ammonium conc. for uptake
  ! p_chp       [mmolP/m3]       Half saturation phosphate conc. for uptake
  ! p_ruen      [1/d]            Relaxation timescale for N uptake/remin.
  ! p_ruep      [1/d]            Relaxation timescale for P uptake/remin.
  integer  :: p_version(iiSeaiceBacteria)
  integer, parameter ::       BACT1=1,BACT2=2,BACT3=3
  real(RLEN)  :: p_q10(iiSeaiceBacteria)
  real(RLEN)  :: p_chdo(iiSeaiceBacteria)
  real(RLEN)  :: p_sd(iiSeaiceBacteria)
  real(RLEN)  :: p_sd2(iiSeaiceBacteria)
  real(RLEN)  :: p_suhU1(iiSeaiceBacteria)
  real(RLEN)  :: p_sulU1(iiSeaiceBacteria)
  real(RLEN)  :: p_suU6(iiSeaiceBacteria)
  real(RLEN)  :: p_sum(iiSeaiceBacteria)
  real(RLEN)  :: p_pu_ra(iiSeaiceBacteria)
  real(RLEN)  :: p_pu_ra_o(iiSeaiceBacteria)
  real(RLEN)  :: p_srs(iiSeaiceBacteria)
  real(RLEN)  :: p_qncSBA(iiSeaiceBacteria)
  real(RLEN)  :: p_qpcSBA(iiSeaiceBacteria)
  real(RLEN)  :: p_chn(iiSeaiceBacteria)
  real(RLEN)  :: p_chp(iiSeaiceBacteria)
  real(RLEN)  :: p_ruen(iiSeaiceBacteria)
  real(RLEN)  :: p_ruep(iiSeaiceBacteria)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitSeaiceBac
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSeaiceBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /SeaiceBac_parameters/ p_version, p_q10, p_chdo, p_sd, p_sd2, p_suhU1, &
    p_sulU1, p_suU6, p_sum, p_pu_ra, p_pu_ra_o, p_srs, p_qpcSBA, &
    p_qncSBA, p_chn, p_chp, p_ruen, p_ruep
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading SeaiceBac parameters.."
  open(NMLUNIT,file='Seaice_Ecology.nml',status='old',action='read',err=100)
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
100 call error_msg_prn(NML_OPEN,"ModuleSeaiceBac.f90","Seaice_Ecology.nml")
101 call error_msg_prn(NML_READ,"ModuleSeaiceBac.f90","SeaiceBac_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSeaiceBac

  end module mem_SeaiceBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
