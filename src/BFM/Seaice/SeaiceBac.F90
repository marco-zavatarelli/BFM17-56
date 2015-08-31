#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: SeaiceBac
!
! DESCRIPTION
!   This process describes the dynamics of bacteria in sea ice
!    
!
! !INTERFACE
  subroutine SeaiceBacDynamics(bac)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE_ICE
#else
  use mem, ONLY: D2STATE_ICE, iiSeaiceBacteria, ppSeaiceBacteria, &
    U1c, U6c, F2o, I4n, I1p, iiN, iiP, iiC
  use mem, ONLY: ppU6c, ppU6n, ppU6p, ppU1c, ppF3c, &
    ppU1n, ppU1p, ppF2o, ppI4n, ppI1p, ppI3n, qpcSOM, qncSOM,&
    ETB, qncSBA, qpcSBA, qpcSOM, qncSOM, NO_BOXES_ICE, iiIce, iiIce,  &
    flux_vector, iiU1, iiU6
#endif

  use constants,  ONLY: MW_C, ONE_PER_DAY
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro, p_small
  use mem_globalfun,   ONLY: eTq_vector, MM_power_vector, insw_vector, &
                             MM_vector
  use mem_SeaiceBac



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_power_vector, insw_vector

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
  integer,intent(IN)  :: bac

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer       :: ppbacc, ppbacn, ppbacp
  real(RLEN),dimension(NO_BOXES_ICE)  :: bacc
  real(RLEN),dimension(NO_BOXES_ICE)  :: et
  real(RLEN),dimension(NO_BOXES_ICE)  :: eO2
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rd
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruU1c
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruU1n
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruU1p
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruU6c
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruU6p
  real(RLEN),dimension(NO_BOXES_ICE)  :: ruU6n
  real(RLEN),dimension(NO_BOXES_ICE)  :: ren
  real(RLEN),dimension(NO_BOXES_ICE)  :: rep
  real(RLEN),dimension(NO_BOXES_ICE)  :: rut
  real(RLEN),dimension(NO_BOXES_ICE)  :: rum
  real(RLEN),dimension(NO_BOXES_ICE)  :: run
  real(RLEN),dimension(NO_BOXES_ICE)  :: rug
  real(RLEN),dimension(NO_BOXES_ICE)  :: cuU6
  real(RLEN),dimension(NO_BOXES_ICE)  :: cuU1
  real(RLEN),dimension(NO_BOXES_ICE)  :: iI1p
  real(RLEN),dimension(NO_BOXES_ICE)  :: iNIn
  real(RLEN),dimension(NO_BOXES_ICE)  :: eI4n,eI1p
  real(RLEN),dimension(NO_BOXES_ICE)  :: iN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppbacc = ppSeaiceBacteria(bac,iiC)
  ppbacn = ppSeaiceBacteria(bac,iiN)
  ppbacp = ppSeaiceBacteria(bac,iiP)
  bacc = D2STATE_ICE(ppbacc,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(ETB(:),  p_q10(bac))
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Oxygen environment: bacteria are both aerobic and anaerobic
  ! To provide a faster switching between the two metabolic pathways the
  ! oxygen regulating factor eO2 is cubic
  ! (eq. 19 in Vichi et al., 2004)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2 = MM_power_vector(max(p_small,F2o(:)),  p_chdo(bac),3)
 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! External nutrient limitation (used by some parametrizations)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eI4n = MM_vector(I4n(:), p_chn(bac))
  eI1p = MM_vector(I1p(:), p_chp(bac))
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Mortality:
  !   1. first order mortality: p_sd 
  !   2. density dependent mortality due to virus infection: p_sd2
  !
  !   It is assumed the mortality is distributed in the same way over
  !   LOC (U1) and detritus (U6) s for phytoplankton and microzooplankton.
  !   using the p_pe_R1x parameters defined in Param
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd  =  ( p_sd(bac)*et + (p_sd2(bac)*bacc))*bacc

  call flux_vector( iiIce, ppbacc,ppU6c, rd*( ONE- p_pe_R1c) )
  call flux_vector( iiIce, ppbacn,ppU6n, rd* qncSBA(bac,:)*( ONE- p_pe_R1n) )
  call flux_vector( iiIce, ppbacp,ppU6p, rd* qpcSBA(bac,:)*( ONE- p_pe_R1p) )

  call flux_vector( iiIce, ppbacc,ppU1c, rd* p_pe_R1c )
  call flux_vector( iiIce, ppbacn,ppU1n, rd* qncSBA(bac,:)* p_pe_R1n )
  call flux_vector( iiIce, ppbacp,ppU1p, rd* qpcSBA(bac,:)* p_pe_R1p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Substrate availability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version(bac) )

    case ( BACT1,BACT2 )  ! Vichi et al. (2004,2007), Lazzari et al. (2012) 

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
      ! Nutrient limitation (intracellular, eq. 51 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      iNIn  =   min(ONE, max(ZERO, qncSBA(bac,:)/ p_qncSBA(bac)))  !Nitrogen
      iI1p  =   min(ONE, max(ZERO, qpcSBA(bac,:)/ p_qpcSBA(bac)))  !Phosphorus
      iN  =   min(iI1p, iNIn)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria (eq. 50 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  = iN*et*p_sum(bac)*bacc

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! correction of substrate quality depending on nutrient content
      ! (eq. 52 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      cuU1 = min(ONE, qpcSOM(iiU1,:)/p_qpcSBA(bac), qncSOM(iiU1,:)/ p_qncSBA(bac))
      cuU6 = min(ONE, qpcSOM(iiU6,:)/p_qpcSBA(bac), qncSOM(iiU6,:)/ p_qncSBA(bac))

  end select


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate the realized substrate uptake rate depending on the
  ! type of detritus and quality (cuRx)
  ! See eq 27 in Vichi et al., 2004 for R2
  ! and eq 6 in Polimene et al., 2006 for R3 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ruU1c = (p_suhU1(bac)*cuU1 + p_sulU1(bac)*(ONE-cuU1))*U1c(:)
  ruU6c = p_suU6(bac)* cuU6* U6c(:)
  rut   = p_small + ruU6c+ ruU1c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria (eq. 50 Vichi et al. 2007)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rug = min( rum, rut )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruU6c = rug*ruU6c/rut
  ruU1c = rug*ruU1c/rut
  call flux_vector( iiIce, ppU6c, ppbacc, ruU6c )
  call flux_vector( iiIce, ppU1c, ppbacc, ruU1c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruU1n = qncSOM(iiU1,:)*ruU1c
  ruU6n = qncSOM(iiU6,:)*ruU6c
  call flux_vector( iiIce, ppU1n, ppbacn, ruU1n )
  call flux_vector( iiIce, ppU6n, ppbacn, ruU6n )

  ruU1p = qpcSOM(iiU1,:)*ruU1c
  ruU6p = qpcSOM(iiU6,:)*ruU6c
  call flux_vector( iiIce, ppU1p, ppbacp, ruU1p )
  call flux_vector( iiIce, ppU6p, ppbacp, ruU6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Aerobic Respiration 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrc = (p_pu_ra(bac)+ p_pu_ra_o(bac)*( ONE- eO2))* rug+ p_srs(bac)* bacc* et
  call flux_vector( iiIce, ppbacc, ppF3c, rrc )
  call flux_vector( iiIce, ppF2o, ppF2o,-rrc/ MW_C) 

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Net Production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  run = rug - rrc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !             Fluxes from bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version(bac))

    case ( BACT1 ) ! Vichi et al. 2007

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! There is no Carbon excretion in Vichi et al. 2007
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
      ! This rate is assumed to occur with a timescale p_ruen=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ren  =  (qncSBA(bac,:) - p_qncSBA(bac))*bacc*p_ruen(bac)
      call flux_vector(iiIce, ppbacn, ppI4n,       ren*insw_vector( ren))
      call flux_vector(iiIce, ppI4n, ppbacn, -eI4n*ren*insw_vector(-ren))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
      ! This rate is assumed to occur with a timescale of p_ruep=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rep  =  (qpcSBA(bac,:) - p_qpcSBA(bac))*bacc*p_ruep(bac)
      call flux_vector(iiIce, ppbacp, ppI1p,       rep*insw_vector( rep))
      call flux_vector(iiIce, ppI1p, ppbacp, -eI1p*rep*insw_vector(-rep))


  end select


  end subroutine SeaiceBacDynamics

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
