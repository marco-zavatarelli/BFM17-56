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
  subroutine SeaiceZooDynamics(zoo,  ppzooc, ppzoon, ppzoop)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants,  ONLY:MW_C
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE_ICE
#else
  use mem, ONLY: D2STATE_ICE, T1c, T1n, T1p, X1c, F2o, F3c, U1c, U6c, U1n, U6n, &
    U1p, U6p, I4n, I1p, SeaiceAlgae, SeaiceZoo, SeaiceBacteria
#endif
  use mem, ONLY: ppT1c, ppT1n, ppT1p, ppX1c, ppF2o, ppF3c, ppU1c, ppU6c, Depth,&
    ppU1n, ppU6n, ppU1p, ppU6p, ppI4n, ppI1p, ppSeaiceAlgae, ppSeaiceZoo, &
    ETB, eO2mO2, qncSBA, qpcSBA, qncSAL, qpcSAL, qncSZO, qpcSZO, qlcSAL, qscSAL, &
    iiSeaiceBacteria, iiSeaiceAlgae, iiSeaiceZoo, iiS1, iiC, iiN, iiP, iiL, &
    NO_BOXES_ICE, iiIce, iiPel, flux_vector,fixed_quota_flux_vector
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
  integer,intent(IN) :: ppzooc
  integer,intent(IN) :: ppzoon
  integer,intent(IN) :: ppzoop

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
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2007 the BFM team
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
  real(RLEN),dimension(NO_BOXES_ICE) :: zooc
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  real(RLEN),dimension(NO_BOXES_ICE)  :: CORROX
  real(RLEN),dimension(NO_BOXES_ICE)  :: put_u
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
  real(RLEN),dimension(NO_BOXES_ICE)  :: efood
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrsc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrac
  real(RLEN),dimension(NO_BOXES_ICE)  :: reac
  real(RLEN),dimension(NO_BOXES_ICE)  :: rdc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrtc
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruTIc
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruSIc
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruXIc
  real(RLEN),dimension(NO_BOXES_ICE,iiSeaiceBacteria)  :: rumTIc
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
  real(RLEN),dimension(NO_BOXES_ICE)  :: ren
  real(RLEN),dimension(NO_BOXES_ICE)  :: pu_ra
  real(RLEN),dimension(NO_BOXES_ICE)  :: r
  real(RLEN),dimension(NO_BOXES_ICE,iiSeaiceAlgae)  :: rumSIc
  real(RLEN),dimension(NO_BOXES_ICE,iiSeaiceZoo)  :: rumXIc
  real(RLEN),dimension(NO_BOXES_ICE)  :: flS1U6s
  real(RLEN),dimension(NO_BOXES_ICE)  :: tfluxc
  real(RLEN),dimension(NO_BOXES_ICE)  :: tfluxn
  real(RLEN),dimension(NO_BOXES_ICE)  :: tfluxp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  zooc = D2STATE_ICE(ppzooc,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETB(:),  p_q10(zoo))
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  CORROX  =   ONE+ p_chro(zoo)
  eF2  =   min(  ONE,  CORROX* MM_vector(  eO2mO2(:),   p_chro(zoo)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Available food, etc...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 do i = 1 , ( iiSeaiceBacteria)
  rumTIc(:, i)  =   p_suSBA(zoo, i)*SeaiceBacteria(i,iiC)* &
    SeaiceBacteria(i,iiC)/( SeaiceBacteria(i,iiC)+ p_minfood(zoo))
  rumc  =   rumTIc(:, i)
  rumn  =   rumTIc(:, i)* qncSBA(i, :)
  rump  =   rumTIc(:, i)* qpcSBA(i, :)
 end do

  do i = 1 , ( iiSeaiceAlgae)
    rumSIc(:, i) = p_suSAL(zoo,i)* SeaiceAlgae(i,iiC)* &
      SeaiceAlgae(i,iiC)/( SeaiceAlgae(i,iiC)+ p_minfood(zoo))
    rumc  =   rumc+ rumSIc(:, i)
    rumn  =   rumn+ rumSIc(:, i)* qncSAL(i,:)
    rump  =   rump+ rumSIc(:, i)* qpcSAL(i,:)
  end do

  do i = 1 , ( iiSeaiceZoo)
    rumXIc(:, i) = p_suSZO(zoo, i)* &
      SeaiceZoo(i,iiC)* SeaiceZoo(i,iiC)/( SeaiceZoo(i,iiC) + &
      p_minfood(zoo))
    rumc  =   rumc+ rumXIc(:, i)
    rumn  =   rumn+ rumXIc(:, i)* qncSZO(i, :)
    rump  =   rump+ rumXIc(:, i)* qpcSZO(i, :)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  efood  =   MM_vector(  rumc,  p_chuc(zoo))
  rugc  =   p_sum(zoo)* et* zooc* efood

  r  =   min(  rumn/ p_qncSZO(zoo),  rump/p_qpcSZO(zoo))
  pu_ra  =   max(  p_pu_ra(zoo),  ONE- r/ (rumc+ p_small))
  put_u  =   rugc/ (rumc+p_small)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes into microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  do i = 1 , ( iiSeaiceBacteria)

   ruTIc  =   put_u* rumTIc(:, i)
   call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzooc,ppT1c,ppzooc, ruTIc ,tfluxC)
   call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoon,ppT1n,ppzoon, &
                                                       ruTIc* qncSBA(i, :),tfluxN)
   call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoop,ppT1p,ppzoop, &
                                                       ruTIc* qpcSBA(i, :),tfluxP)
   rugn  =   ruTIc* qncSBA(i, :)
   rugp  =   ruTIc* qpcSBA(i, :)

  end do

  do i = 1 , ( iiSeaiceAlgae)

    ruSIc  =   put_u* rumSIc(:, i)
    call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzooc,ppSeaiceAlgae(i,iiC),&
                                                      ppzooc, ruSIc ,tfluxC)
    call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoon,ppSeaiceAlgae(i,iiN),&
                                           ppzoon, ruSIc* qncSAL(i,:) ,tfluxN)
    call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoop,ppSeaiceAlgae(i,iiP),&
                                           ppzoop, ruSIc* qpcSAL(i,:) ,tfluxP)
    ! Chl is transferred to the sink
    call flux_vector( iiIce, ppSeaiceAlgae(i,iiL),ppSeaiceAlgae(i,iiL),-( &
      ruSIc* qlcSAL(i,:)) )
    if ( i== iiS1) then
      ! S1s is directly transferred to U6s
      flS1U6s(:)  =   flS1U6s(:)+ ruSIc* qscSAL(i,:)
    end if

    rugn  =   rugn+ ruSIc* qncSAL(i,:)
    rugp  =   rugp+ ruSIc* qpcSAL(i,:)
  end do

  do i = 1 , ( iiSeaiceZoo)

    ruXIc  =   put_u* rumXIc(:, i)
    ! intra-group predation is not computed
    if ( i/= zoo) then
      call fixed_quota_flux_vector( check_fixed_quota,iiIce,ppzooc,ppSeaiceZoo(i,iiC),&
                                                          ppzooc, ruXIc,tfluxC )
      call fixed_quota_flux_vector( check_fixed_quota,iiIce,ppzoon,ppSeaiceZoo(i,iiN),&
                                              ppzoon, ruXIc* qncSZO(i,:) ,tfluxN)
      call fixed_quota_flux_vector( check_fixed_quota,iiIce,ppzoop,ppSeaiceZoo(i,iiP),&
                                              ppzoop, ruXIc* qpcSZO(i,:) ,tfluxP)
    end if

    rugn  =   rugn+ ruXIc* qncSZO(i,:)
    rugp  =   rugp+ ruXIc* qpcSZO(i,:)
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

  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzooc,ppzooc,ppF3c, &
                                                  rrtc,tfluxC )
  call flux_vector( iiIce, ppF2o,ppF2o,-( rrtc/ MW_C) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Mortality (rdc) + Excetion (reac)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rdc  =  (( ONE- eF2)* p_sdo(zoo)+ p_sd(zoo))* zooc
  reac  =   rugc* p_pu_ea(zoo)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes due to mortality and excetion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rric  =  ( reac+ rdc)
  rr1c  =   rric* p_pe_R1c
  rr6c  =   rric*( ONE- p_pe_R1c)

  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzooc,ppzooc,ppU1c, rr1c,tfluxC)
  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzooc,ppzooc,ppU6c, rr6c,tfluxC)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !     Nutrient dynamics in microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrin  =   rugn* p_pu_ea(zoo)+ rdc* qncSZO(zoo,:)
  rr1n  =   rrin* p_pe_R1n
  rr6n  =   rrin- rr1n

  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoon,ppzoon,ppU1n, rr1n ,tfluxN)
  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoon,ppzoon,ppU6n, rr6n ,tfluxN)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Phosphorus dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrip  =   rugp* p_pu_ea(zoo)+ rdc* qpcSZO(zoo,:)
  rr1p  =   rrip* p_pe_R1p
  rr6p  =   rrip- rr1p

  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoop,ppzoop,ppU1p, rr1p ,tfluxP)
  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoop,ppzoop,ppU6p, rr6p ,tfluxP)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved nutrient dynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  runc  =   max(  ZERO,  rugc*( ONE- p_pu_ea(zoo))- rrac)
  runn  =   max(  ZERO,  rugn*( ONE- p_pu_ea(zoo))+ rrsc* qncSZO(zoo, :))
  runp  =   max(  ZERO,  rugp*( ONE- p_pu_ea(zoo))+ rrsc* qpcSZO(zoo, :))

  ren  =   max(  ZERO,  runn/( p_small+ runc)- p_qncSZO(zoo))* runc
  rep  =   max(  ZERO,  runp/( p_small+ runc)- p_qpcSZO(zoo))* runc
  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoon,ppzoon,ppI4n, ren ,tfluxN)
  call fixed_quota_flux_vector( check_fixed_quota,iiIce, ppzoop,ppzoop,ppI1p, rep ,tfluxP)


  r=tfluxC*p_qncSZO(zoo)
  call fixed_quota_flux_vector( check_fixed_quota,-iiN,0,0,0,r,tfluxN)
  r=tfluxC*p_qpcSZO(zoo)
  call fixed_quota_flux_vector( check_fixed_quota,-iiP,0,0,0,r,tfluxP)


  end subroutine SeaiceZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
