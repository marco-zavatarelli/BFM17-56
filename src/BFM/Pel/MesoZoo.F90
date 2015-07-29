#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MesoZoo
!
! DESCRIPTION
!   This submodel describes the carbon dynamics and associated
!   nutrient dynamics in mesozooplankton 
!
! !INTERFACE
  subroutine MesoZooDynamics(zoo)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN, ONE, ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D3STATE, O2o, N1p, N4n, R6c, R6p, R2c, &
    R6n, PhytoPlankton, MicroZooPlankton, MesoZooPlankton
  use mem, ONLY: Depth, ppO2o, ppN1p, ppN4n, ppR6c, ppR6n, ppR6p, ppR6s, &
    ppPhytoPlankton, ppMicroZooPlankton, ppMesoZooPlankton, ETW, &
    qncPPY, qpcPPY, qlcPPY, qscPPY, qncMIZ, qpcMIZ, qncMEZ, qpcMEZ, iiPhytoPlankton, &
    iiMicroZooPlankton, iiMesoZooPlankton, iiC, iiN, iiP, iiL, iiS, NO_BOXES, &
    iiBen, iiPel, flux_vector,fixed_quota_flux_vector
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c
#endif
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF, qfcPPY, ppR6f
#endif
#endif
#ifdef BFM_GOTM
  use mem, ONLY: jnetMeZc
#endif
  use mem_Param,  ONLY: p_small,check_fixed_quota
  use constants,ONLY: MIN_VAL_EXPFUN, MW_C
  use mem_MesoZoo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector, MM_power_vector

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
!   First ERSEM version by N. Broekhuizen and A.D. Bryant
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
  integer  :: i
  integer  :: ppzooc, ppzoon, ppzoop
  integer,dimension(NO_BOXES)  :: nut_lim
  integer, save :: first =0
  real(RLEN),allocatable,save,dimension(:) :: sut,temp_p,temp_n,rumc,rugc,eo,  &
                                       et,rrs_c,rrs_n,rrs_p,rut_c, &
                                       rut_n,rut_p,rd_c,rd_n,rd_p,sdo,rdo_c,  &
                                       rdo_n,rdo_p,ret_c,ret_n,ret_p,ru_c, &
                                       ru_n,ru_p,pu_e_n,pu_e_p,prI,pe_R6c

  real(RLEN),allocatable,save,dimension(:) :: pe_N1p,pe_N4n,ruPPYc,ruMIZc,ruMEZc,rq6c, &
                                       rq6n,rq6p,rrc,ren,rep,tfluxc,tfluxn,    &
                                       tfluxp,zooc,zoop,zoon
  real(RLEN),allocatable,save,dimension(:,:) :: PPYc,MIZc,MEZc
  real(RLEN),allocatable,save,dimension(:) :: net,r
  integer :: AllocStatus, DeallocStatus
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     first=1
     allocate(eo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eo"
     allocate(PPYc(NO_BOXES,4),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating PPYc"
       allocate(MIZc(NO_BOXES,2),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating MIZc"
       allocate(MEZc(NO_BOXES,2),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating MEZc"
       allocate(zooc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating zooc"
       allocate(zoop(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating zoop"
       allocate(zoon(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating zoon"
     allocate(tfluxp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxp"
     allocate(tfluxn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxn"
     allocate(tfluxc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxc"
     allocate(rep(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rep"
     allocate(ren(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ren"
     allocate(rrc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrc"
     allocate(rq6p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rq6p"
     allocate(rq6n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rq6n"
     allocate(rq6c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rq6c"
     allocate(ruMEZc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruMEZc"
     allocate(ruMIZc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruMIZc"
     allocate(ruPPYc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruPPYc"
     allocate(pe_N4n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pe_N4n"
     allocate(pe_N1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pe_N1p"
     allocate(pe_R6c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pe_R6c"
     allocate(prI(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating prI"
     allocate(pu_e_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pu_e_p"
     allocate(pu_e_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pu_e_n"
     allocate(ru_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ru_p"
     allocate(ru_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ru_n"
     allocate(ru_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ru_c"
     allocate(ret_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ret_p"
     allocate(ret_n(NO_BOXES),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating ret_n"
     allocate(ret_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ret_c"
     allocate(rdo_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rdo_p"
     allocate(rdo_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rdo_n"
     allocate(rdo_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rdo_c"
     allocate(sdo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sdo"
     allocate(rd_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rd_p"
     allocate(rd_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rd_n"
     allocate(rd_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rd_c"
     allocate(rut_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rut_p"
     allocate(rut_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rut_n"
     allocate(rut_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rut_c"
     allocate(et(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating et"
     allocate(sut(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sut"
     allocate(temp_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating temp_p"
     allocate(temp_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating temp_n"
     allocate(rumc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumc"
     allocate(rugc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rugc"
     allocate(net(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating net"
     allocate(r(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating r"
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Copy state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppzooc = ppMesoZooPlankton(zoo,iiC)
  ppzoon = ppMesoZooPlankton(zoo,iiN)
  ppzoop = ppMesoZooPlankton(zoo,iiP)
  zooc = D3STATE(ppzooc,:)
  zoon = D3STATE(ppzoon,:)
  zoop = D3STATE(ppzoop,:)

  ! temporary variables in case check_fixed_quota=1
  tfluxc = ZERO
  tfluxn = ZERO
  tfluxp = ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature and oxygen response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eo = MM_power_vector(max(p_small,O2o(:)), p_clO2o(zoo),3)
  et = eTq_vector(ETW(:), p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total potential food given the non-dim prey availability
  ! with loops over all LFGs.
  ! There is no parameter for capture efficiency in mesozooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rumc = ZERO
  do i = 1, iiPhytoPlankton
    PPYc(:,i) = p_paPPY(zoo,i)*PhytoPlankton(i,iiC)
    rumc = rumc + PPYc(:,i)
  end do
  do i = 1, iiMicroZooPlankton
    MIZc(:,i) = p_paMIZ(zoo,i)*MicroZooPlankton(i,iiC)
    rumc = rumc + MIZc(:,i)
  end do
  do i = 1, iiMesoZooPlankton
    MEZc(:,i) = p_paMEZ(zoo,i)*MesoZooPlankton(i,iiC)
    rumc = rumc + MEZc(:,i)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake rate (eq 38 Vichi et al. 2007) and 
  ! specific uptake rate considering potentially available food (sut)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  = et*p_sum(zoo)*MM_vector(p_vum(zoo)*rumc, p_sum(zoo))*zooc
  sut = rugc/(p_small + rumc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptakes from every LFG
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rut_c = ZERO
  rut_n = ZERO
  rut_p = ZERO
  ! Phytoplankton
  do i = 1, iiPhytoPlankton
    ruPPYc = sut*PPYc(:,i)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
                ppPhytoPlankton(i,iiC), ppzooc, ruPPYc          , tfluxc)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, &
                ppPhytoPlankton(i,iiN), ppzoon, ruPPYc*qncPPY(i,:), tfluxn)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoop, &
                ppPhytoPlankton(i,iiP), ppzoop, ruPPYc*qpcPPY(i,:), tfluxp)
    rut_c = rut_c + ruPPYc
    rut_n = rut_n + ruPPYc*qncPPY(i,:)
    rut_p = rut_p + ruPPYc*qpcPPY(i,:)
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
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
               ppMicroZooPlankton(i,iiC), ppzooc, ruMIZc           , tfluxc )
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, &
               ppMicroZooPlankton(i,iiN), ppzoon, ruMIZc*qncMIZ(i,:), tfluxn)
    call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoop, &
               ppMicroZooPlankton(i,iiP), ppzoop, ruMIZc*qpcMIZ(i,:), tfluxp)
    rut_c = rut_c + ruMIZc
    rut_n = rut_n + ruMIZc*qncMIZ(i,:)
    rut_p = rut_p + ruMIZc*qpcMIZ(i,:)
  end do
  ! Mesozooplankton
  do i = 1, iiMesoZooPlankton
    ruMEZc = sut*MEZc(:, i)
    ! Note that intra-group predation (cannibalism) is not added as a flux
    if ( i/= zoo ) then
      call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
                 ppMesoZooPlankton(i,iiC), ppzooc, ruMEZc          , tfluxc )
      call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                 ppMesoZooPlankton(i,iiN), ppzoon, ruMEZc*qncMEZ(i,:), tfluxn)
      call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
                 ppMesoZooPlankton(i,iiP), ppzoop, ruMEZc*qpcMEZ(i,:), tfluxp)
    end if
    rut_c = rut_c + ruMEZc
    rut_n = rut_n + ruMEZc*qncMEZ(i,:)
    rut_p = rut_p + ruMEZc*qpcMEZ(i,:)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity respiration and basal metabolism
  ! First compute the the energy cost of ingestion
  ! 1 - assimilation - egestion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  prI = ONE - p_puI(zoo) - p_peI(zoo)
  rrc = prI*rut_c + p_srs(zoo)*et*zooc
  call flux_vector(iiPel, ppO2o, ppO2o, -rrc/MW_C)
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzooc, &
                               ppzooc, ppO3c, rrc, tfluxc )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Specific rates of low oxygen mortality
  ! and Density dependent mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rdo_c = p_sdo(zoo)*(ONE-eo)*et*zooc
  rd_c  = p_sd(zoo)*zooc**p_sds(zoo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total egestion including pellet production 
  ! Eq. 40 and 44 Vichi et al. 2007
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rq6c = p_peI(zoo)*rut_c + rdo_c + rd_c
  rq6n = p_peI(zoo)*rut_n + qncMEZ(zoo,:)*(rdo_c + rd_c)
  rq6p = p_peI(zoo)*rut_p + qpcMEZ(zoo,:)*(rdo_c + rd_c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Check the assimilation rate for Carbon, Nitrogen and Phosphorus
  ! Note that activity respiration does not involve nutrient utilization
  ! so more nutrients than carbon are taken up.
  ! Then compute P:C and N:C ratios in the assimilation rate
  ! Eq 41 in Vichi et al. 2007 (there is an error in the denominator,
  ! the \Iota_c should be \Iota_i, with i=n,p)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ru_c = p_puI(zoo)*rut_c
  ru_n = (p_puI(zoo) + prI)* rut_n
  ru_p = (p_puI(zoo) + prI)* rut_p
  pu_e_n  =   ru_n/( p_small+ ru_c)
  pu_e_p  =   ru_p/( p_small+ ru_c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Eliminate the excess of the non-limiting constituent
  ! Determine whether C, P or N is the limiting element and assign the
  ! value to variable nut_lim
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  nut_lim = iiC
  temp_p  = pu_e_p/qpcMEZ(zoo,:)
  temp_n  = pu_e_n/qncMEZ(zoo,:)

  WHERE ( temp_p<temp_n .OR. abs(temp_p-temp_n)<p_small ) 
      WHERE ( pu_e_p< qpcMEZ(zoo,:) )
        nut_lim = iiP
      END WHERE
  ELSEWHERE
      WHERE ( pu_e_n<qncMEZ(zoo,:) )
        nut_lim = iiN
      END WHERE
  END WHERE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute the correction terms depending on the limiting constituent
  ! Eq. 42 Vichi et al 2007 for a combination of N and P limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  WHERE     ( nut_lim==iiC )
      pe_R6c = ZERO
      pe_N1p = max(ZERO, (ONE - p_peI(zoo))*rut_p - p_qpcMEZ(zoo)*ru_c)
      pe_N4n = max(ZERO, (ONE - p_peI(zoo))*rut_n - p_qncMEZ(zoo)*ru_c)
  ELSEWHERE ( nut_lim==iiP )
      pe_N1p = ZERO
      pe_R6c = max(ZERO, ru_c - (ONE - p_peI(zoo))*rut_p/p_qpcMEZ(zoo))
      pe_N4n = max( ZERO, (ONE - p_peI(zoo))*rut_n - p_qncMEZ(zoo)*(ru_c - pe_R6c))
  ELSEWHERE ( nut_lim==iiN )
      pe_N4n = ZERO
      pe_R6c = max(ZERO, ru_c - (ONE - p_peI(zoo))*rut_n/p_qncMEZ(zoo))
      pe_N1p = max(ZERO, (ONE - p_peI(zoo))*rut_p - p_qpcMEZ(zoo)*(ru_c - pe_R6c))
  END WHERE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient remineralization 
  ! basal metabolism + excess of non-limiting nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ren = p_srs(zoo)*et*eo*zoon + pe_N4n
  rep = p_srs(zoo)*et*eo*zoop + pe_N1p
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoop, &
                               ppzoop, ppN1p, rep, tfluxp)
  call fixed_quota_flux_vector(check_fixed_quota, iiPel, ppzoon, &
                               ppzoon, ppN4n, ren, tfluxn)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes to particulate organic matter
  ! Add the correction term for organic carbon release in case of
  ! nutrient limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rq6c = rq6c + pe_R6c
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                             ppzooc,ppR6c, rq6c ,tfluxc )
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
                             ppzoop,ppR6p, rq6p ,tfluxp)
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                             ppzoon,ppR6n, rq6n ,tfluxn)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! The following part is computed only if zooplankton has fixed 
  ! nutrient quota and check_fixed_quota is set to 1
  ! It controls all nutrient imbalances and gives a warning in case 
  ! there are nutrient leaks
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( check_fixed_quota == 1 ) then
     ren = tfluxC*p_qncMEZ(zoo)
     call fixed_quota_flux_vector( check_fixed_quota,-iiN,0,0,0,ren,tfluxN)
     rep = tfluxC*p_qpcMEZ(zoo)
     call fixed_quota_flux_vector( check_fixed_quota,-iiP,0,0,0,rep,tfluxP)
  end if

  end subroutine MesoZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
