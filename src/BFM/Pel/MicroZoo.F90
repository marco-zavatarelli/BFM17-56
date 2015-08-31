#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MicroZoo
!
! DESCRIPTION
!  Microzooplankton dynamics
!
! !INTERFACE
  subroutine MicroZooDynamics(zoo)
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
  use mem, ONLY: D3STATE, O2o, R1c, R6c, R1n, R6n, &
    R2c,R1p, R6p, N4n, N1p, PhytoPlankton, MicroZooPlankton, PelBacteria
  use mem, ONLY: ppPelBacteria, ppO2o, ppR1c, ppR6c, ppR6s, Depth,&
    ppR1n, ppR6n, ppR1p, ppR6p, ppN4n, ppN1p, ppPhytoPlankton, ppMicroZooPlankton, &
    ETW, eO2mO2, qncPBA, qpcPBA, qncPPY, qpcPPY, qncMIZ, qpcMIZ, iiPelBacteria, &
    qlcPPY, qscPPY, iiPhytoPlankton, iiMicroZooPlankton, iiC, iiN, iiP, iiL, iiS, &
    NO_BOXES, iiBen, iiPel, flux_vector,fixed_quota_flux_vector
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c
#endif
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF, qfcPPY, ppR6f
#endif
#endif
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_small,check_fixed_quota
  use mem_MicroZoo

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
!   First ERSEM version by H. Baretta-Bekker and J.W. Baretta
!   Additional parametrizations by P. Ruardij and M. Vichi 
!   Dynamical allocation by G. Mattia 
!
! !REVISION_HISTORY
!   !
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer       :: i
  integer       :: ppzooc, ppzoon, ppzoop
  integer, save :: first =0
  integer       :: AllocStatus, DeallocStatus
  real(RLEN),allocatable,save,dimension(:) :: sut,et,eO2,rumc,rumn,rump,  &
                                         rugc,rugn,rugp,runc,runn,runp, &
                                         rrsc,rrac,reac,rdc,rrtc,ruPBAc,ruPPYc,  &
                                         ruMIZc,rric,rr1c,rr6c,rr1p,rr1n, &
                                         rrip,rr6p,rep,rrin,zooc
  real(RLEN),allocatable,save,dimension(:)    :: rr6n,ren,pu_ra,r,tfluxc,tfluxn,tfluxp
  real(RLEN),allocatable,save,dimension(:,:)  :: PBAc,PPYc,MIZc
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if (first==0) then
     first=1
     allocate(PBAc(NO_BOXES,iiPelBacteria),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating PBAc"
     allocate(PPYc(NO_BOXES,iiPhytoPlankton),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating PPYc"
     allocate(MIZc(NO_BOXES,iiMicroZooPlankton),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating MIZc"
     allocate(sut(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sut,"
     allocate(et(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating et"
     allocate(eO2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eO2"
     allocate(rumc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumc"
     allocate(rumn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumn"
     allocate(rump(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rump"
     allocate(rugc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rugc"
     allocate(rugn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rugn"
     allocate(rugp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rugp"
     allocate(runc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runc"
     allocate(runn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runn"
     allocate(runp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runp"
     allocate(rrsc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrsc"
     allocate(rrac(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrac"
     allocate(reac(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating reac"
     allocate(rdc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rdc"
     allocate(rrtc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrtc"
     allocate(ruPBAc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruPBAc"
     allocate(ruPPYc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruPPYc"
     allocate(ruMIZc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruMIZc"
     allocate(rric(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rric"
     allocate(rr1c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr1c"
     allocate(rr6c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6c"
     allocate(rr1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr1p"
     allocate(rr1n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr1n"
     allocate(zooc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating zooc"
     allocate(rrip(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrip"
     allocate(rr6p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6p"
     allocate(rep(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rep"
     allocate(rrin(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrin"
     allocate(rr6n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6n"
     allocate(ren(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ren"
     allocate(pu_ra(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pu_ra"
     allocate(r(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating r"
     allocate(tfluxc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxc"
     allocate(tfluxn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxn"
     allocate(tfluxp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxp"
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppzooc = ppMicroZooPlankton(zoo,iiC)
  ppzoon = ppMicroZooPlankton(zoo,iiN)
  ppzoop = ppMicroZooPlankton(zoo,iiP)
  zooc = D3STATE(ppzooc,:)

  tfluxc=ZERO
  tfluxn=ZERO
  tfluxp=ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et = eTq_vector(ETW(:), p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2 = min(ONE, MM_vector(O2o(:), p_chro(zoo)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total potential food given the non-dim prey availability
  ! and capture efficiency with loops over all LFGs.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rumc   = ZERO
  rumn   = ZERO
  rump   = ZERO
  do i = 1 ,iiPelBacteria
     PBAc(:,i) = p_paPBA(zoo,i)*PelBacteria(i,iiC)* &
                 MM_vector(PelBacteria(i,iiC), p_minfood(zoo))
     rumc = rumc + PBAc(:,i)
     rumn = rumn + PBAc(:,i)*qncPBA(i,:)
     rump = rump + PBAc(:,i)*qpcPBA(i,:)
  end do

  do i = 1 ,iiPhytoPlankton
     PPYc(:,i) = p_paPPY(zoo,i)*PhytoPlankton(i,iiC)* &
                   MM_vector(PhytoPlankton(i,iiC), p_minfood(zoo))
    rumc = rumc + PPYc(:,i)
    rumn = rumn + PPYc(:,i)*qncPPY(i,:)
    rump = rump + PPYc(:,i)*qpcPPY(i,:)
  end do

  do i = 1, iiMicroZooPlankton
     MIZc(:,i) = p_paMIZ(zoo,i)*MicroZooPlankton(i,iiC)* &
                   MM_vector(MicroZooPlankton(i,iiC), p_minfood(zoo))
    rumc = rumc + MIZc(:,i)
    rumn = rumn + MIZc(:,i)*qncMIZ(i,:)
    rump = rump + MIZc(:,i)*qpcMIZ(i,:)
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
  ! Bacterioplankton
  rugn = ZERO
  rugp = ZERO

  do i = 1, iiPelBacteria
    ruPBAc = sut*PBAc(:,i)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
               ppPelBacteria(i,iiC), ppzooc, ruPBAc, tfluxc )
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, &
               ppPelBacteria(i,iiN), ppzoon, ruPBAc*qncPBA(i,:), tfluxn)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoop, &
               ppPelBacteria(i,iiP), ppzoop, ruPBAc*qpcPBA(i,:), tfluxp)
    rugn = rugn + ruPBAc*qncPBA(i,:)
    rugp = rugp + ruPBAc*qpcPBA(i,:)
  end do
  ! Phytoplankton
  do i = 1, iiPhytoPlankton
    ruPPYc = sut*PPYc(:,i)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
                ppPhytoPlankton(i,iiC), ppzooc, ruPPYc          , tfluxc)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, &
                ppPhytoPlankton(i,iiN), ppzoon, ruPPYc*qncPPY(i,:), tfluxn)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoop, &
                ppPhytoPlankton(i,iiP), ppzoop, ruPPYc*qpcPPY(i,:), tfluxp)
    rugn = rugn + ruPPYc*qncPPY(i,:)
    rugp = rugp + ruPPYc*qpcPPY(i,:)
    ! Chl is transferred to the infinite sink
    call flux_vector(iiPel, ppPhytoPlankton(i,iiL), &
               ppPhytoPlankton(i,iiL), -ruPPYc*qlcPPY(i,:))
    ! silicon constituent is transferred to biogenic silicate
    if ( ppPhytoPlankton(i,iiS) .gt. 0 ) & 
       call flux_vector(iiPel, ppPhytoPlankton(i,iiS), ppR6s, ruPPYc*qscPPY(i,:))                              
#ifdef INCLUDE_PELFE
    ! Fe constituent is transferred to particulate iron
    if ( ppPhytoPlankton(i,iiF) .gt. 0 ) & 
       call flux_vector(iiPel, ppPhytoPlankton(i,iiF), ppR6f, ruPPYc*qfcPPY(i,:))
#endif
  end do
  ! Microzooplankton
  do i = 1, iiMicroZooPlankton
    ruMIZc = sut*MIZc(:,i)
    ! Note that intra-group predation (cannibalism) is not added as a flux
    if ( i/= zoo) then
       call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
               ppMicroZooPlankton(i,iiC), ppzooc, ruMIZc           , tfluxc )
       call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, &
               ppMicroZooPlankton(i,iiN), ppzoon, ruMIZc*qncMIZ(i,:), tfluxn)
       call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoop, &
               ppMicroZooPlankton(i,iiP), ppzoop, ruMIZc*qpcMIZ(i,:), tfluxp)
    end if
    rugn = rugn + ruMIZc*qncMIZ(i,:)
    rugp = rugp + ruMIZc*qpcMIZ(i,:)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes from microzooplankton
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
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
                               ppzooc, ppO3c, rrtc, tfluxC)
  call flux_vector(iiPel, ppO2o, ppO2o, -rrtc/MW_C)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality (rdc) + Activity Excretion (reac)
  ! and partitioning between particulate and dissolved
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rdc  = ((ONE - eO2)*p_sdo(zoo) + p_sd(zoo))*zooc
  reac = rugc*(ONE - p_pu(zoo))*p_pu_ea(zoo)
  rric = reac + rdc
  rr1c = rric*p_pe_R1c
  rr6c = rric*(ONE - p_pe_R1c)
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, ppzooc, &
                               ppR1c, rr1c, tfluxC)
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, ppzooc, &
                               ppR6c, rr6c, tfluxC)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !     Nutrient dynamics in microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrin = rugn*p_pu_ea(zoo) + rdc*qncMIZ(zoo,:)
  rr1n = rrin*p_pe_R1n
  rr6n = rrin - rr1n
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, ppzoon, &
                               ppR1n, rr1n ,tfluxN)
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, ppzoon, &
                               ppR6n, rr6n ,tfluxN)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Phosphorus dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrip = rugp*p_pu_ea(zoo) + rdc*qpcMIZ(zoo,:)
  rr1p = rrip*p_pe_R1p
  rr6p = rrip - rr1p
  call fixed_quota_flux_vector(check_fixed_quota,iiPel, ppzoop, ppzoop, &
                               ppR1p, rr1p, tfluxP)
  call fixed_quota_flux_vector(check_fixed_quota,iiPel, ppzoop, ppzoop, &
                               ppR6p, rr6p, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved nutrient dynamics
  ! Compare the quota of the net growth rates with the optimal quota
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  runc = max(ZERO, rugc*(ONE - p_pu_ea(zoo)) - rrac)
  runn = max(ZERO, rugn*(ONE - p_pu_ea(zoo)) + rrsc*qncMIZ(zoo,:))
  runp = max(ZERO, rugp*(ONE - p_pu_ea(zoo)) + rrsc*qpcMIZ(zoo,:))
  ren  = max(ZERO, runn/(p_small + runc) - p_qncMIZ(zoo))* runc
  rep  = max(ZERO, runp/(p_small + runc) - p_qpcMIZ(zoo))* runc
  call fixed_quota_flux_vector(check_fixed_quota,iiPel, ppzoon, ppzoon, &
                               ppN4n, ren, tfluxN)
  call fixed_quota_flux_vector(check_fixed_quota,iiPel, ppzoop, ppzoop, &
                               ppN1p, rep, tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! The following part is computed only if zooplankton has fixed 
  ! nutrient quota and check_fixed_quota is set to 1
  ! It controls all nutrient imbalances and gives a warning in case 
  ! there are nutrient leaks
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( check_fixed_quota == 1 ) then
     r = tfluxC*p_qncMIZ(zoo)
     call fixed_quota_flux_vector(check_fixed_quota, -iiN, 0, 0, 0, r, tfluxN)
     r = tfluxC*p_qpcMIZ(zoo)
     call fixed_quota_flux_vector(check_fixed_quota, -iiP, 0, 0, 0, r, tfluxP)
  end if

  end subroutine MicroZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
