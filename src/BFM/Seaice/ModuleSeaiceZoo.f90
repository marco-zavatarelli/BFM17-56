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
  ! MicroZoo PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  real(RLEN)  :: p_q10(iiSeaiceZoo)  ! Q10 value
  real(RLEN)  :: p_srs(iiSeaiceZoo)  ! Respiration rate at 10 degrees Celsius
  real(RLEN)  :: p_sum(iiSeaiceZoo) ! Max. rel daily uptake as a fraction of biomass
  real(RLEN)  :: p_sdo(iiSeaiceZoo) ! Mortality due to oxygen limitation
  real(RLEN)  :: p_sd(iiSeaiceZoo)  ! Temperature independent mortality
  real(RLEN)  :: p_pu_ra(iiSeaiceZoo) ! Activity respiration
  real(RLEN)  :: p_pu_ea(iiSeaiceZoo)  ! Activity excretion
  real(RLEN)  :: p_chro(iiSeaiceZoo) ! Oxygen saturation where respiration is 0.5
  real(RLEN)  :: p_chuc(iiSeaiceZoo) ! Food concentration where total uptake rate is 0.5
  real(RLEN)  :: p_minfood(iiSeaiceZoo)  ! Concentration below which feeding on a particular
                                                !  foodsource is depressed
  real(RLEN)  :: p_suTI(iiSeaiceZoo,iiSeaiceBacteria) ! /day   #relative B1 uptake by zoo
  real(RLEN)  :: p_qnXI(iiSeaiceZoo)  ! Maximum quotum P
  real(RLEN)  :: p_qpXI(iiSeaiceZoo) ! Maximum quotum N
  real(RLEN)  :: p_suSI(iiSeaiceZoo,iiSeaiceAlgae)    ! /day   #relative P uptake by zoo
  real(RLEN)  :: p_suXI(iiSeaiceZoo,iiSeaiceZoo)! /day   #relative Z uptake by zoo
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitSeaiceZoo
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSeaiceZoo()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /SeaiceZoo_parameters/ p_q10, p_srs, p_sum, p_sdo, p_sd, p_pu_ra, &
    p_pu_ea, p_chro, p_chuc, p_minfood, p_suSI, p_suXI, p_suTI, p_qpXI, p_qnXI
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading SeaiceZoo parameters.."
open(NMLUNIT,file='SeaiceZoo.nml',status='old',action='read',err=100)
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
100 call error_msg_prn(NML_OPEN,"InitSeaiceZoo.f90","SeaiceZoo.nml")
101 call error_msg_prn(NML_READ,"InitSeaiceZoo.f90","SeaiceZoo_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSeaiceZoo
  end module mem_SeaiceZoo
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
