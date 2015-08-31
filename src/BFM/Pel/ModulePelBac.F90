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
  use mem,  ONLY: iiPelBacteria
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
  integer,private     :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !NAMELIST PelBacteria_parameters
  !-------------------------------------------------------------------------!
  !  PELAGIC BACTERIA
  !
  ! NAME         [UNIT]/KIND            DESCRIPTION
  ! p_version   integer         Switch for bacteria parameterization
  !                              1 : Baretta-Bekker et al. 1995;
  !                                  Vichi et al., 2007
  !                              2 : Vichi et al., 2004
  !                              3 : Polimene et al., 2006
  ! p_q10                        Q10-value (temperature dependency)
  ! p_chdo      [mmol/m3]        Half-saturation constant for O2 limitation
  ! p_sd        [1/d]            Specific mortality rate
  ! p_sd2       [1/d]            Density dependent specific mortality rate
  ! p_suhR1     [1/d]            Specific potential uptake for nutrient-rich DOM
  ! p_sulR1     [1/d]            Specific potential uptake for nutrient-poor DOM
  ! p_suR2      [1/d]            Specific potential uptake for semi-labile DOC
  ! p_suR3      [1/d]            Specific potential uptake for semi-refractory DOC
  ! p_suR6      [1/d]            Specific potential uptake for POM (1/d)
  ! p_sum       [1/d]            Potential specific growth rate
  ! p_pu_ra     [-]              Activity respiration fraction
  ! p_pu_ra_o   [-]              Additional respiration fraction at low O2 conc
  ! p_srs       [1/d]            Specific rest respiration
  ! p_qncPBA    [mmolN/mgC]      Optimal N/C ratio 
  ! p_qpcPBA    [mmolP/mgC]      Optimal P/C ratio 
  ! p_qlnc      [mmolN/mgC]      Minimal N/C ratio 
  ! p_qlpc      [mmolP/mgC]      Minimal P/C ratio 
  ! p_qun       [mmolN/mgC/day]  Membrane affinity for N 
  ! p_qup       [mmolP/mgC/day]  Membrane affinity for N 
  ! p_chn       [mmolN/m3]       Half saturation ammonium conc. for uptake
  ! p_chp       [mmolP/m3]       Half saturation phosphate conc. for uptake
  ! p_ruen      [1/d]            Relaxation timescale for N uptake/remin.
  ! p_ruep      [1/d]            Relaxation timescale for P uptake/remin.
  ! p_rec       [1/d]            Relaxation timescale for semi-labile excretion
  ! p_pu_ea_R3  [-]              Excretion of semi-refractory DOC
  integer     :: p_version(iiPelBacteria)
  integer, parameter ::       BACT1=1,BACT2=2,BACT3=3
  real(RLEN)  :: p_q10(iiPelBacteria)
  real(RLEN)  :: p_chdo(iiPelBacteria)
  real(RLEN)  :: p_sd(iiPelBacteria)
  real(RLEN)  :: p_sd2(iiPelBacteria)
  real(RLEN)  :: p_suhR1(iiPelBacteria)
  real(RLEN)  :: p_sulR1(iiPelBacteria)
  real(RLEN)  :: p_suR2(iiPelBacteria)
  real(RLEN)  :: p_suR6(iiPelBacteria)
  real(RLEN)  :: p_suR3(iiPelBacteria)
  real(RLEN)  :: p_sum(iiPelBacteria)
  real(RLEN)  :: p_pu_ra(iiPelBacteria)
  real(RLEN)  :: p_pu_ra_o(iiPelBacteria)
  real(RLEN)  :: p_srs(iiPelBacteria)
  real(RLEN)  :: p_qncPBA(iiPelBacteria)
  real(RLEN)  :: p_qpcPBA(iiPelBacteria)
  real(RLEN)  :: p_qlnc(iiPelBacteria)
  real(RLEN)  :: p_qlpc(iiPelBacteria)
  real(RLEN)  :: p_qun(iiPelBacteria)
  real(RLEN)  :: p_qup(iiPelBacteria)
  real(RLEN)  :: p_chn(iiPelBacteria)
  real(RLEN)  :: p_chp(iiPelBacteria)
  real(RLEN)  :: p_ruen(iiPelBacteria)
  real(RLEN)  :: p_ruep(iiPelBacteria)
  real(RLEN)  :: p_rec(iiPelBacteria)
  real(RLEN)  :: p_pu_ea_R3(iiPelBacteria)

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
  namelist /PelBacteria_parameters/ p_version, p_q10, p_chdo, p_sd, p_sd2, p_suhR1, &
    p_sulR1, p_suR2, p_suR6, p_sum, p_pu_ra, p_pu_ra_o, p_pu_ea_R3, p_srs, &
    p_suR3, p_qpcPBA, p_qlpc, p_qncPBA, p_qlnc, p_qun, p_qup, p_chn, p_chp, &
    p_ruen, p_ruep, p_rec

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading PelBac parameters.."
  open(NMLUNIT,file='Pelagic_Ecology.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=PelBacteria_parameters,err=101)
  close(NMLUNIT)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Check consistency of parameters according to the parametrization
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  write(LOGUNIT,*) "#  using p_version=",p_version
  do i=1,iiPelBacteria
     write(LOGUNIT,*) "#  Checking PelBacteria parameters for group:",i
     select case ( p_version(i) )
       case ( BACT3 ) ! Polimene et al. (2006)
         p_sulR1(i) = ZERO
         write(LOGUNIT,*) "#  forcing p_sulR1=0"
         if (p_pu_ea_R3(i) + p_pu_ra(i) .GT. 0.3_RLEN) then
           write(LOGUNIT,*)"#  Warning: Bacterial growth efficiency is lower than 0.3!"
           write(LOGUNIT,*)"#  The release of capsular material is possibly larger than p_pu_ra/4"
         end if
       case ( BACT1 ) ! Vichi et al. 2007
         p_sulR1(i) = ZERO
         p_suR2(i) = ZERO
         p_suR3(i) = ZERO
         write(LOGUNIT,*) "#  forcing p_sulR1,p_suR2,p_suR3=0"
       case ( BACT2 ) ! Vichi et al. 2004
         p_suR3(i) = ZERO
         write(LOGUNIT,*) "#  forcing p_suR3=0"
     end select
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Write parameter list to the log
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  write(LOGUNIT,*) "#  Namelist is:"
  write(LOGUNIT,nml=PelBacteria_parameters)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPelBac.f90","Pelagic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitPelBac.f90","PelBacteria_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPelBac
  end module mem_PelBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
