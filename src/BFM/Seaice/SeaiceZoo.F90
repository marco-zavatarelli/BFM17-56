#include "DEBUG.h"
#include "INCLUDE.h"

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
  subroutine SeaiceZooDynamics(zoo)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants,  ONLY:MW_C
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D2STATE_ICE, SeaiceAlgae, SeaiceZoo, SeaiceBacteria
  use mem, ONLY:  F2o, ppF2o, ppF3c, ppU1c, ppU6c, ppU6s, &
    ppU1n, ppU6n, ppU1p, ppU6p, ppI4n, ppI1p, ppSeaiceAlgae, ppSeaiceZoo, ppSeaiceBacteria, &
    ETB, qncSBA, qpcSBA, qncSAL, qpcSAL, qncSZO, qpcSZO, qlcSAL, qscSAL, &
    iiSeaiceBacteria, iiSeaiceAlgae, iiSeaiceZoo, iiS, iiC, iiN, iiP, iiL, &
    NO_BOXES_ICE, iiIce, iiIce, flux_vector,fixed_quota_flux_vector
#endif
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_small,check_fixed_quota
  use mem_SeaiceZoo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector,MM_power_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: zoo
!  
!
! !AUTHORS
!
!
! !REVISION_HISTORY
!
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer       :: ppzooc, ppzoon, ppzoop
  real(RLEN),dimension(NO_BOXES_ICE) :: zooc
  real(RLEN),dimension(NO_BOXES_ICE)  :: sut
  real(RLEN),dimension(NO_BOXES_ICE)  :: et
  real(RLEN),dimension(NO_BOXES_ICE)  :: eF2
  real(RLEN),dimension(NO_BOXES_ICE)  :: rumc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rumn
  real(RLEN),dimension(NO_BOXES_ICE)  :: rump
  real(RLEN),dimension(NO_BOXES_ICE)  :: rugc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rugn
  real(RLEN),dimension(NO_BOXES_ICE)  :: rugp
  real(RLEN),dimension(NO_BOXES_ICE)  :: runc
  real(RLEN),dimension(NO_BOXES_ICE)  :: runn
  real(RLEN),dimension(NO_BOXES_ICE)  :: runp
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrsc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrac
  real(RLEN),dimension(NO_BOXES_ICE)  :: reac
  real(RLEN),dimension(NO_BOXES_ICE)  :: rdc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrtc
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruSBAc
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruSALc
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruSZOc
  real(RLEN),dimension(NO_BOXES_ICE,iiSeaiceBacteria)  :: SBAc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rric
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr1c
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr6c
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr1p
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr1n
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrip
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr6p
  real(RLEN),dimension(NO_BOXES_ICE)  :: rep
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrin
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr6n
  real(RLEN),dimension(NO_BOXES_ICE)  :: ren, r
  real(RLEN),dimension(NO_BOXES_ICE,iiSeaiceAlgae)  :: SALc
  real(RLEN),dimension(NO_BOXES_ICE,iiSeaiceZoo)  :: SZOc
  real(RLEN),dimension(NO_BOXES_ICE)  :: tfluxc
  real(RLEN),dimension(NO_BOXES_ICE)  :: tfluxn
  real(RLEN),dimension(NO_BOXES_ICE)  :: tfluxp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppzooc = ppSeaiceZoo(zoo,iiC)
  ppzoon = ppSeaiceZoo(zoo,iiN)
  ppzoop = ppSeaiceZoo(zoo,iiP)
  zooc = D2STATE_ICE(ppzooc,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(ETB(:),  p_q10(zoo))
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eF2 = min(ONE, MM_vector(F2o(:), p_chro(zoo)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total potential food given the non-dim prey availability
  ! and capture efficiency with loops over all LFGs.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rumc = ZERO 
  rumn = ZERO
  rump = ZERO
 do i = 1 , iiSeaiceBacteria
    SBAc(:,i) = p_paSBA(zoo,i)*SeaiceBacteria(i,iiC)* &
                MM_vector(SeaiceBacteria(i,iiC), p_minfood(zoo))
    rumc = rumc + SBAc(:,i)
    rumn = rumn + SBAc(:,i)*qncSBA(i,:)
    rump = rump + SBAc(:,i)*qpcSBA(i,:)
 end do

  do i = 1 ,iiSeaiceAlgae
     SALc(:,i) = p_paSAL(zoo,i)*SeaiceAlgae(i,iiC)* &
                 MM_vector(SeaiceAlgae(i,iiC), p_minfood(zoo))
     rumc = rumc + SALc(:,i)
     rumn = rumn + SALc(:,i)*qncSAL(i,:)
     rump = rump + SALc(:,i)*qpcSAL(i,:)
  end do

  do i = 1, iiSeaiceZoo
     SZOc(:,i) = p_paSZO(zoo,i)* SeaiceZoo(i,iiC)* &
                MM_vector(SeaiceZoo(i,iiC), p_minfood(zoo))
     rumc = rumc + SZOc(:,i)
     rumn = rumn + SZOc(:,i)*qncSZO(i,:)
     rump = rump + SZOc(:,i)*qpcSZO(i,:)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake rate (eq 38 Vichi et al. 2007) and 
  ! specific uptake rate considering potentially available food (sut)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  = et*p_sum(zoo)*MM_vector(rumc, p_chuc(zoo))*zooc
  sut = rugc/rumc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptakes from every LFG
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugn = ZERO
  rugp = ZERO

  do i = 1 , iiSeaiceBacteria
    ruSBAc =  sut* SBAc(:, i)
    call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzooc, &
               ppSeaiceBacteria(i,iiC),ppzooc,ruSBAc,tfluxC)
    call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzoon, &
               ppSeaiceBacteria(i,iiN),ppzoon,ruSBAc*qncSBA(i,:),tfluxN)
    call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzoop, &
               ppSeaiceBacteria(i,iiP),ppzoop,ruSBAc*qpcSBA(i,:),tfluxP)
    rugn = rugn + ruSBAc*qncSBA(i,:)
    rugp = rugp + ruSBAc*qpcSBA(i,:)
  end do

  do i = 1, iiSeaiceAlgae
    ruSALc = sut*SALc(:,i)
    call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzooc, &
                ppSeaiceAlgae(i,iiC),ppzooc,ruSALc,tfluxC)
    call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzoon, &
                ppSeaiceAlgae(i,iiN),ppzoon,ruSALc*qncSAL(i,:),tfluxN)
    call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzoop, &
                ppSeaiceAlgae(i,iiP),ppzoop,ruSALc*qpcSAL(i,:),tfluxP)
    rugn = rugn + ruSALc*qncSAL(i,:)
    rugp = rugp + ruSALc*qpcSAL(i,:)
    ! Chl is transferred to the infinite sink
    call flux_vector(iiIce, ppSeaiceAlgae(i,iiL), &
               ppSeaiceAlgae(i,iiL),-ruSALc*qlcSAL(i,:))
    ! silicon constituent is transferred to biogenic silicate
    if (ppSeaiceAlgae(i,iiS) > 0) &
       call flux_vector(iiIce, ppSeaiceAlgae(i,iiS), ppU6s,-ruSALc*qscSAL(i,:))
  end do

  do i = 1, iiSeaiceZoo
    ruSZOc = sut* SZOc(:,i)
    ! Note that intra-group predation (cannibalism) is not added as a flux
    if ( i/= zoo) then
      call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzooc, &
                                   ppSeaiceZoo(i,iiC),ppzooc,ruSZOc,tfluxC)
      call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzoon, &
                                   ppSeaiceZoo(i,iiN),ppzoon,ruSZOc*qncSZO(i,:),tfluxN)
      call fixed_quota_flux_vector(check_fixed_quota,iiIce,ppzoop, &
                                   ppSeaiceZoo(i,iiP),ppzoop,ruSZOc*qpcSZO(i,:),tfluxP)
    end if
    rugn = rugn + ruSZOc*qncSZO(i,:)
    rugp = rugp + ruSZOc*qpcSZO(i,:)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !       Fluxes from seaice fauna
  ! The metabolic balance is the following:
  ! Ingestion = Growth + Excretion + Respiration
  ! Assimilation efficiency p_pu = G/I
  ! Excretion E = I*p_pu_ea
  ! therefore R = (1-p_pu-p_pu_ea)*I
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Rest, activity and total respiration fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrsc = p_srs(zoo)*et*zooc
  ! the activity respiration is derived from the other constant parameters
  rrac = rugc*(ONE - p_pu(zoo) - p_pu_ea(zoo))
  rrtc = rrsc + rrac
  call fixed_quota_flux_vector(check_fixed_quota, iiIce, ppzooc, &
                               ppzooc, ppF3c, rrtc, tfluxC)
  call flux_vector(iiIce, ppF2o, ppF2o, -rrtc/MW_C)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality (rdc) + Activity Excretion (reac)
  ! and partitioning between particulate and dissolved
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rdc  = ((ONE - eF2)*p_sdo(zoo) + p_sd(zoo))*zooc
  reac = rugc*(ONE - p_pu(zoo))*p_pu_ea(zoo)
  rric = reac + rdc
  rr1c = rric*p_pe_R1c
  rr6c = rric*(ONE - p_pe_R1c)
  call fixed_quota_flux_vector(check_fixed_quota, iiIce, ppzooc, ppzooc, &
                               ppU1c, rr1c, tfluxC)
  call fixed_quota_flux_vector(check_fixed_quota, iiIce, ppzooc, ppzooc, &
                               ppU6c, rr6c, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !     Nutrient dynamics in microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrin = rugn*p_pu_ea(zoo) + rdc*qncSZO(zoo,:)
  rr1n = rrin*p_pe_R1n
  rr6n = rrin - rr1n
  call fixed_quota_flux_vector(check_fixed_quota, iiIce,ppzoon, &
                               ppzoon, ppU1n, rr1n ,tfluxN)
  call fixed_quota_flux_vector(check_fixed_quota, iiIce,ppzoon, &
                               ppzoon, ppU6n, rr6n ,tfluxN)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Phosphorus dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrip = rugp*p_pu_ea(zoo) + rdc*qpcSZO(zoo,:)
  rr1p = rrip*p_pe_R1p
  rr6p = rrip - rr1p
  call fixed_quota_flux_vector(check_fixed_quota, iiIce,ppzoop, &
                               ppzoop, ppU1p, rr1p ,tfluxP)
  call fixed_quota_flux_vector(check_fixed_quota, iiIce,ppzoop, &
                               ppzoop, ppU6p, rr6p ,tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved nutrient dynamics
  ! Compare the quota of the net growth rates with the optimal quota
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  runc = max(ZERO, rugc*(ONE - p_pu_ea(zoo)) - rrac)
  runn = max(ZERO, rugn*(ONE - p_pu_ea(zoo)) + rrsc*qncSZO(zoo,:))
  runp = max(ZERO, rugp*(ONE - p_pu_ea(zoo)) + rrsc*qpcSZO(zoo,:))
  ren = max(ZERO,  runn/(p_small + runc) - p_qncSZO(zoo))* runc
  rep = max(ZERO,  runp/(p_small + runc) - p_qpcSZO(zoo))* runc
  call fixed_quota_flux_vector(check_fixed_quota, iiIce, ppzoon, &
                               ppzoon, ppI4n, ren ,tfluxN)
  call fixed_quota_flux_vector(check_fixed_quota, iiIce, ppzoop, &
                               ppzoop, ppI1p, rep ,tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! The following part is computed only if zooplankton has fixed 
  ! nutrient quota and check_fixed_quota is set to 1
  ! It controls all nutrient imbalances and gives a warning in case 
  ! there are nutrient leaks
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( check_fixed_quota == 1 ) then
     r=tfluxC*p_qncSZO(zoo)
     call fixed_quota_flux_vector( check_fixed_quota,-iiN,0,0,0,r,tfluxN)
     r=tfluxC*p_qpcSZO(zoo)
     call fixed_quota_flux_vector( check_fixed_quota,-iiP,0,0,0,r,tfluxP)
  end if

  end subroutine SeaiceZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
