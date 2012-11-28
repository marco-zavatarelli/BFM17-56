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

!  
!
! !AUTHORS
!   the ERSEM group, Marcello Vichi, JWB, HBB
!
!
!
! !REVISION_HISTORY
!   !
!
!
! COPYING
!   
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
  ! Phyto PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !  ---------------- Physiological parameters -----------------
  !
  real(RLEN)  :: p_q10(iiPhytoPlankton)  ! Doubling temperature
  real(RLEN)  :: p_temp(iiPhytoPlankton)=ZERO ! Cut-off for temperature factor
  real(RLEN)  :: p_sum(iiPhytoPlankton)  ! Maximal productivity at 10 degrees C
  real(RLEN)  :: p_srs(iiPhytoPlankton)  ! Respiration rate at 10 degrees C
  real(RLEN)  :: p_sdmo(iiPhytoPlankton)  ! Max.specific nutrient-stress lysis rate
  real(RLEN)  :: p_thdo(iiPhytoPlankton)  ! Half value for nutrient stress lysis
  real(RLEN)  :: p_seo(iiPhytoPlankton)  ! Extra lysis rate for P4
  real(RLEN)  :: p_pu_ea(iiPhytoPlankton)  ! Fraction of pp excreted as PLOC/PDET
  real(RLEN)  :: p_pu_ra(iiPhytoPlankton)  ! Activity respiration rate
  real(RLEN)  :: p_switchR1R2(iiPhytoPlankton)  ! Switch for R1-R2 excretion
  !
  !  ---------------- Nutrient parameters in phytoplankton -----------------
  !
  logical  :: p_netgrowth(iiPhytoPlankton)=.TRUE.  ! logical switch for nut. limitation growth
  integer  :: p_limnut(iiPhytoPlankton)  ! switch for nut. limitation (Liebig is default)
  real(RLEN)  :: p_qnlc(iiPhytoPlankton)
  real(RLEN)  :: p_qnRc(iiPhytoPlankton)
  real(RLEN)  :: p_xqn(iiPhytoPlankton)
  real(RLEN)  :: p_qplc(iiPhytoPlankton)
  real(RLEN)  :: p_qpRc(iiPhytoPlankton)
  real(RLEN)  :: p_xqp(iiPhytoPlankton)
  real(RLEN)  :: p_qslc(iiPhytoPlankton) ! Minimum quotum Si in PI
  real(RLEN)  :: p_qsRc(iiPhytoPlankton) ! Reference quotum Si in PI
  real(RLEN)  :: p_sheo(iiPhytoPlankton) 
  real(RLEN)  :: p_qus(iiPhytoPlankton)  ! affinity of PI for Si
  real(RLEN)  :: p_qun(iiPhytoPlankton)
  real(RLEN)  :: p_qup(iiPhytoPlankton)
  real(RLEN)  :: p_lN4(iiPhytoPlankton)
  real(RLEN)  :: p_chPs(iiPhytoPlankton) ! half-value of SIO4-lim (mmol Si m-3)
  real(RLEN)  :: p_Contois(iiPhytoplankton) ! parameter for Contois
  real(RLEN)  :: p_esNI(iiPhytoPlankton) ! Nutrient stress threshold for Sinking
  real(RLEN)  :: p_res(iiPhytoPlankton)  ! Sinking velocity (m/d)
  !
  !  ---------------- Light parameters in phytoplankton -----------------
  !
  real(RLEN)  :: p_alpha_chl(iiPhytoPlankton)  ! Initial slope P-I curve
  real(RLEN)  :: p_sdchl(iiPhytoPlankton)      ! Specific turnover rate for Chla [d-1]
  real(RLEN)  :: p_EpEk_or(iiPhytoplankton)    ! optimal E_PAR/E_K
  real(RLEN)  :: p_tochl_relt(iiPhytoplankton) ! relaxation rate for chl:C [d-1]
#ifdef INCLUDE_PELFE
  !
  !  ---------------- Iron parameters in phytoplankton -----------------
  !
  real(RLEN)  :: p_qflc(iiPhytoPlankton)
  real(RLEN)  :: p_qfRc(iiPhytoPlankton)
  real(RLEN)  :: p_xqf(iiPhytoPlankton)
  real(RLEN)  :: p_quf(iiPhytoPlankton)
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitPhyto

  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPhyto()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Phyto_parameters/ p_q10, p_sum, p_srs, p_sdmo, p_seo, p_pu_ea, &
                              p_temp, p_netgrowth,p_limnut, &
                              p_pu_ra, p_qnlc, p_qplc, p_qslc, &
                              p_qnRc, p_qpRc, p_qsRc, &
                              p_qun, p_qup, p_qus, &
                              p_xqn, p_xqp, p_sheo, &
                              p_esNI, p_thdo, p_res, p_lN4, p_chPs, &
                              p_Contois, p_EpEk_or, p_tochl_relt,   &
                              p_switchR1R2,                         &
                              p_alpha_chl, p_sdchl

#ifdef INCLUDE_PELFE
  namelist /Phyto_parameters_iron/ p_qflc, p_qfRc, p_xqf, p_quf
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  write(LOGUNIT,*) "#  Reading Phyto parameters.."
  open(NMLUNIT,file='Phyto.nml',status='old',action='read',err=100)
  read(NMLUNIT,nml=Phyto_parameters,err=101)
#ifdef INCLUDE_PELFE
  read(NMLUNIT,nml=Phyto_parameters_iron,err=101)
#endif
  close(NMLUNIT)
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
100 call error_msg_prn(NML_OPEN,"InitPhyto.f90","Phyto.nml")
101 call error_msg_prn(NML_READ,"InitPhyto.f90","Phyto_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPhyto
  end module mem_Phyto
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
