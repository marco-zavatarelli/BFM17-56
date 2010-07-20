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
!   !
!
! !INTERFACE
  subroutine MicroZooDynamics(zoo,  ppzooc, ppzoon, ppzoop)
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
  use mem, ONLY: D3STATE, B1c, B1n, B1p, O2o, R1c, R6c, R1n, R6n, &
    R2c,R1p, R6p, N4n, N1p, PhytoPlankton, MicroZooPlankton
  use mem, ONLY: ppB1c, ppB1n, ppB1p, ppO2o, ppO3c, ppR1c, ppR6c, Depth,&
    ppR1n, ppR6n, ppR1p, ppR6p, ppN4n, ppN1p, ppPhytoPlankton, ppMicroZooPlankton, &
    flP1R6s, ETW, eO2mO2, qnB1c, qpB1c, qnPc, qpPc, qn_mz, qp_mz, jnetMiZc, &
    qlPc, qsPc, iiPhytoPlankton, iiMicroZooPlankton, iiP1, iiC, iiN, iiP, iiL, &
    NO_BOXES, iiBen, iiPel, flux_vector,fixed_quota_flux_vector
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
  integer,intent(IN) :: ppzooc
  integer,intent(IN) :: ppzoon
  integer,intent(IN) :: ppzoop
!
!
! !AUTHORS
!   ERSEM group, Hanneke Baretta-Bekker
!
! !REVISION_HISTORY
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer       :: i
  integer, save :: first =0
  integer       :: AllocStatus, DeallocStatus
  real(RLEN),allocatable,save,dimension(:) :: put_u,et,eO2,rumc,rumn,rump,  &
                                         rugc,rugn,rugp,runc,runn,runp,efood, &
                                         rrsc,rrac,reac,rdc,rrtc,ruB1c,ruPIc,  &
                                         ruZIc,rumB1c,rric,rr1c,rr6c,rr1p,rr1n, &
                                         rrip,rr6p,rep,rrin,zooc
  real(RLEN),allocatable,save,dimension(:)    :: rr6n,ren,pu_ra,r,tfluxc,tfluxn,tfluxp
  real(RLEN),allocatable,save,dimension(:,:)  :: rumPIc,rumZIc
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if (first==0) then
     first=1
     allocate(rumPIc(NO_BOXES,iiPhytoPlankton),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumPIc"
     allocate(rumZIc(NO_BOXES,iiMicroZooPlankton),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumZIc"
     allocate(put_u(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating put_u,"
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
     allocate(efood(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating efood"
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
     allocate(ruB1c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruB1c"
     allocate(ruPIc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruPIc"
     allocate(ruZIc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruZIc"
     allocate(rumB1c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumB1c"
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
  zooc = D3STATE(ppzooc,:)

  tfluxc=ZERO
  tfluxn=ZERO
  tfluxp=ZERO


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETW(:),  p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2  =   min(  ONE,  MM_vector(  O2o(:),   p_chro(zoo)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Available food, etc...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rumB1c  =   p_suB1(zoo)* B1c(:)* B1c(:)/( B1c(:)+ p_minfood(zoo))
  rumc  =   rumB1c
  rumn  =   rumB1c* qnB1c(:)
  rump  =   rumB1c* qpB1c(:)

  do i = 1 , ( iiPhytoPlankton)
    rumPIc(:, i) = p_suPI(zoo,i)* PhytoPlankton(i,iiC)* &
      PhytoPlankton(i,iiC)/( PhytoPlankton(i,iiC)+ p_minfood(zoo))
    rumc  =   rumc+ rumPIc(:, i)
    rumn  =   rumn+ rumPIc(:, i)* qnPc(i,:)
    rump  =   rump+ rumPIc(:, i)* qpPc(i,:)
  end do

  do i = 1 , ( iiMicroZooPlankton)

    rumZIc(:, i) = p_suZI(zoo,i)* &
      MicroZooPlankton(i,iiC)* MicroZooPlankton(i,iiC)/( MicroZooPlankton(i,iiC)+ &
      p_minfood(zoo))
    rumc  =   rumc+ rumZIc(:, i)
    rumn  =   rumn+ rumZIc(:, i)* qn_mz(i,:)
    rump  =   rump+ rumZIc(:, i)* qp_mz(i,:)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  efood  =   MM_vector(  rumc,  p_chuc(zoo))
  rugc  =   p_sum(zoo)* et* zooc* efood

  r  =   min(  rumn/ p_qn_mz(zoo),  rump/ p_qp_mz(zoo))

  pu_ra  =   max(  p_pu_ra(zoo),  ONE- r/ rumc)

  put_u  =   rugc/ rumc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes into microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruB1c  =   put_u* rumB1c
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc,ppB1c,ppzooc, ruB1c ,tfluxC)
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoon,ppB1n,ppzoon, &
                                                       ruB1c* qnB1c(:),tfluxN)
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoop,ppB1p,ppzoop, &
                                                       ruB1c* qpB1c(:),tfluxP)
  rugn  =   ruB1c* qnB1c(:)
  rugp  =   ruB1c* qpB1c(:)

  do i = 1 , ( iiPhytoPlankton)

    ruPIc  =   put_u* rumPIc(:, i)
    call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc,ppPhytoPlankton(i,iiC),&
                                                      ppzooc, ruPIc ,tfluxC)
    call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoon,ppPhytoPlankton(i,iiN),&
                                           ppzoon, ruPIc* qnPc(i,:) ,tfluxN)
    call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoop,ppPhytoPlankton(i,iiP),&
                                           ppzoop, ruPIc* qpPc(i,:) ,tfluxP)
    ! Chl is transferred to the sink
    call flux_vector( iiPel, ppPhytoPlankton(i,iiL),ppPhytoPlankton(i,iiL),-( &
      ruPIc* qlPc(i,:)) )
    if ( i== iiP1) then
      ! P1s is directly transferred to R6s
      ! PhytoPlankton[i].s -> R6.s = ruPIc * qsPc[i];
      flP1R6s(:)  =   flP1R6s(:)+ ruPIc* qsPc(i,:)
    end if

    rugn  =   rugn+ ruPIc* qnPc(i,:)
    rugp  =   rugp+ ruPIc* qpPc(i,:)
  end do

  do i = 1 , ( iiMicroZooPlankton)

    ruZIc  =   put_u* rumZIc(:, i)
    ! intra-group predation is not computed
    if ( i/= zoo) then
      call fixed_quota_flux_vector( check_fixed_quota,iiPel,ppzooc,ppMicroZooPlankton(i,iiC),&
                                                          ppzooc, ruZIc,tfluxC )
      call fixed_quota_flux_vector( check_fixed_quota,iiPel,ppzoon,ppMicroZooPlankton(i,iiN),&
                                              ppzoon, ruZIc* qn_mz(i,:) ,tfluxN)
      call fixed_quota_flux_vector( check_fixed_quota,iiPel,ppzoop,ppMicroZooPlankton(i,iiP),&
                                              ppzoop, ruZIc* qp_mz(i,:) ,tfluxP)
    end if

    rugn  =   rugn+ ruZIc* qn_mz(i,:)
    rugp  =   rugp+ ruZIc* qp_mz(i,:)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !       Fluxes from microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Rest, activity, total respiration fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrsc  =   p_srs(zoo)* et* zooc
  rrac  =   rugc* pu_ra
  rrtc  =   rrsc+ rrac

  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc,ppzooc,ppO3c, &
                                                  rrtc,tfluxC )
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrtc/ MW_C) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality (rdc) + Excretion (reac)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rdc  =  (( ONE- eO2)* p_sdo(zoo)+ p_sd(zoo))* zooc
  reac  =   rugc* p_pu_ea(zoo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes due to mortality and excetion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rric  =  ( reac+ rdc)
  rr1c  =   rric* p_pe_R1c
  rr6c  =   rric*( ONE- p_pe_R1c)

  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc,ppzooc,ppR1c, rr1c,tfluxC)
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc,ppzooc,ppR6c, rr6c,tfluxC)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !     Nutrient dynamics in microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrin  =   rugn* p_pu_ea(zoo)+ rdc* qn_mz(zoo,:)
  rr1n  =   rrin* p_pe_R1n
  rr6n  =   rrin- rr1n
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoon,ppzoon,ppR1n, rr1n ,tfluxN)
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoon,ppzoon,ppR6n, rr6n ,tfluxN)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Phosphorus dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrip  =   rugp* p_pu_ea(zoo)+ rdc* qp_mz(zoo,:)
  rr1p  =   rrip* p_pe_R1p
  rr6p  =   rrip- rr1p

  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoop,ppzoop,ppR1p, rr1p ,tfluxP)
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoop,ppzoop,ppR6p, rr6p ,tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved nutrient dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  runc  =   max(  ZERO,  rugc*( ONE- p_pu_ea(zoo))- rrac)
  runn  =   max(  ZERO,  rugn*( ONE- p_pu_ea(zoo))+ rrsc* qn_mz(zoo, :))
  runp  =   max(  ZERO,  rugp*( ONE- p_pu_ea(zoo))+ rrsc* qp_mz(zoo, :))

  ren  =   max(  ZERO,  runn/( p_small+ runc)- p_qn_mz(zoo))* runc
  rep  =   max(  ZERO,  runp/( p_small+ runc)- p_qp_mz(zoo))* runc
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoon,ppzoon,ppN4n, ren ,tfluxN)
  call fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoop,ppzoop,ppN1p, rep ,tfluxP)

#ifdef BFM_GOTM
  r=tfluxC*p_qn_mz(zoo)
  call fixed_quota_flux_vector( check_fixed_quota,-iiN,0,0,0,r,tfluxN)
  r=tfluxC*p_qp_mz(zoo)
  call fixed_quota_flux_vector( check_fixed_quota,-iiP,0,0,0,r,tfluxP)

  r=rugc-rrac-reac
  jnetMiZc(1)=jnetMiZc(1)+sum(Depth(:)*r)
#endif

  end subroutine MicroZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
