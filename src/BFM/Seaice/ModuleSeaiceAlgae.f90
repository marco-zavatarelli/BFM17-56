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
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
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
  ! Sea ice algae PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !  ---------------- Physiological parameters -----------------
  !
  real(RLEN)  :: p_q10(iiSeaiceAlgae)  ! Doubling temperature
  real(RLEN)  :: p_sum(iiSeaiceAlgae)  ! Maximal productivity at 10 degrees C
  real(RLEN)  :: p_srs(iiSeaiceAlgae)  ! Respiration rate at 10 degrees C
  real(RLEN)  :: p_sdmo(iiSeaiceAlgae)  ! Max.specific nutrient-stress lysis rate
  real(RLEN)  :: p_thdo(iiSeaiceAlgae)  ! Half value for nutrient stress lysis
  real(RLEN)  :: p_seo(iiSeaiceAlgae)  ! Extra lysis rate for P4
  real(RLEN)  :: p_pu_ea(iiSeaiceAlgae)  ! Fraction of pp excreted as PLOC/PDET
  real(RLEN)  :: p_pu_ra(iiSeaiceAlgae)  ! Activity respiration rate
  !
  !  ---------------- Nutrient parameters in sea ice algae -----------------
  !
  logical  :: p_netgrowth(iiSeaiceAlgae)=.TRUE.  ! logical switch for nut. limitation growth
  integer  :: p_limnut(iiSeaiceAlgae)  ! switch for nut. limitation (Liebig is default)
  real(RLEN)  :: p_qnlc(iiSeaiceAlgae)
  real(RLEN)  :: p_qnRc(iiSeaiceAlgae)
  real(RLEN)  :: p_xqn(iiSeaiceAlgae)
  real(RLEN)  :: p_qplc(iiSeaiceAlgae)
  real(RLEN)  :: p_qpRc(iiSeaiceAlgae)
  real(RLEN)  :: p_xqp(iiSeaiceAlgae)
  real(RLEN)  :: p_qslc(iiSeaiceAlgae)  ! Minimum quotum Si in PI
  real(RLEN)  :: p_qsRc(iiSeaiceAlgae)  ! Reference quotum Si in PI
  real(RLEN)  :: p_xqs(iiSeaiceAlgae)
  real(RLEN)  :: p_qus(iiSeaiceAlgae)  ! affinity of PI for Si
  real(RLEN)  :: p_qun(iiSeaiceAlgae)
  real(RLEN)  :: p_qup(iiSeaiceAlgae)
  real(RLEN)  :: p_lN4(iiSeaiceAlgae)
  real(RLEN)  :: p_chPs(iiSeaiceAlgae)
  real(RLEN)  :: p_esII(iiSeaiceAlgae)  ! Nutrient stress threshold for Sinking
  ! p_alpha_chl = 1.0e-5, 0.46e-5*2.0, 2.0e-5, 0.68e-5 # Initial slope P-I curve
  !  Thalassiosira sp. [0.48-0.63]
  !
  !  ------------- Chlorophyll parameters -----------
  real(RLEN)  :: p_alpha_chl(iiSeaiceAlgae)  ! Initial slope P-I curve
  real(RLEN)  :: p_sdchl(iiSeaiceAlgae)  ! Specific turnover rate for Chla [d-1]
  real(RLEN)  :: p_chlII(iiSeaiceAlgae)  ! Nutrient stress threshold for chl turnover
  real(RLEN)  :: p_qchlcSI(iiSeaiceAlgae)  ! Maximum quotum Chla:C [mg Chla (mg C)-1]
  


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitSeaiceAlgae
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSeaiceAlgae()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Seaicealgae_parameters/ p_q10, p_sum, p_srs, p_sdmo, p_seo, p_pu_ea, &
    p_pu_ra, p_qnlc, p_qplc, p_qslc, p_qnRc, p_qpRc, p_qsRc, p_qun, p_qup, &
    p_qus, p_xqn, p_xqp, p_xqs, p_chlII, p_thdo, p_lN4, p_chPs, &
    p_netgrowth, p_limnut, p_alpha_chl, p_sdchl, p_qchlcSI, p_esII
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading Sea ice algae parameters.."
    open(NMLUNIT,file='SeaiceAlgae.nml',status='old',action='read',err=100)
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
100 call error_msg_prn(NML_OPEN,"InitSeaicealgae.f90","Seaicealgae.nml")
101 call error_msg_prn(NML_READ,"InitSeaicealgae.f90","Seaicealgae_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSeaiceAlgae
  end module mem_SeaiceAlgae
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
