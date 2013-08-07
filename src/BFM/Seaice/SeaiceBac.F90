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
  subroutine SeaiceBacDynamics(bac,  ppbacc, ppbacn, ppbacp)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE_ICE
#else
  use mem, ONLY: D2STATE_ICE, T1c, U6c, T1n, U6n, T1p, U6p, U1c, U1n, U1p, F2o, F3c, &
    I4n, I1p, I3n
#endif
  use mem, ONLY: ppT1c, ppU6c, ppT1n, ppU6n, ppT1p, ppU6p, ppU1c, ppF3c, &
    ppU1n, ppU1p, ppF2o, ppN6r, ppI4n, ppI1p, ppI3n, Depth, qpcSDE, qncSDE,&
    ETB, qncSBA, qpcSBA, eO2mO2, qpcSDE, qncSDE, NO_BOXES_XY, iiIce, flux_vector, &
    iiU1, iiU6

  use constants,  ONLY: MW_C, ONE_PER_DAY
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro, p_small
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
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2009 M. Vichi and L. Tedesco
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
  integer,intent(IN)  :: bac
  integer,intent(IN) :: ppbacc
  integer,intent(IN) :: ppbacn
  integer,intent(IN) :: ppbacp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY)  :: bacc
  real(RLEN),dimension(NO_BOXES_XY)  :: runn
  real(RLEN),dimension(NO_BOXES_XY)  :: runp
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: eO2
  real(RLEN),dimension(NO_BOXES_XY)  :: r
  real(RLEN),dimension(NO_BOXES_XY)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY)  :: rd
  real(RLEN),dimension(NO_BOXES_XY)  :: ruU1c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruU1n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruU1p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruU6c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruU6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruU6n
  real(RLEN),dimension(NO_BOXES_XY)  :: cqun3
  real(RLEN),dimension(NO_BOXES_XY)  :: rump
  real(RLEN),dimension(NO_BOXES_XY)  :: rumn
  real(RLEN),dimension(NO_BOXES_XY)  :: rumn3
  real(RLEN),dimension(NO_BOXES_XY)  :: rumn4
  real(RLEN),dimension(NO_BOXES_XY)  :: misp
  real(RLEN),dimension(NO_BOXES_XY)  :: misn
  real(RLEN),dimension(NO_BOXES_XY)  :: rupp
  real(RLEN),dimension(NO_BOXES_XY)  :: rupn
  real(RLEN),dimension(NO_BOXES_XY)  :: ren
  real(RLEN),dimension(NO_BOXES_XY)  :: rep
  real(RLEN),dimension(NO_BOXES_XY)  :: rut
  real(RLEN),dimension(NO_BOXES_XY)  :: rum
  real(RLEN),dimension(NO_BOXES_XY)  :: run
  real(RLEN),dimension(NO_BOXES_XY)  :: sun
  real(RLEN),dimension(NO_BOXES_XY)  :: rug
  real(RLEN),dimension(NO_BOXES_XY)  :: suU1
  real(RLEN),dimension(NO_BOXES_XY)  :: suU1n
  real(RLEN),dimension(NO_BOXES_XY)  :: suU1p
  real(RLEN),dimension(NO_BOXES_XY)  :: cuU6
  real(RLEN),dimension(NO_BOXES_XY)  :: cuU1
  real(RLEN),dimension(NO_BOXES_XY)  :: iI1p
  real(RLEN),dimension(NO_BOXES_XY)  :: iNIn
  real(RLEN),dimension(NO_BOXES_XY)  :: iN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  bacc = D2STATE_ICE(ppbacc,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETB(:),  p_q10(bac))
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Mortality:
  !   1. first order mortality: old definition
  !   2. density dependent mortality due to virus infection
  !
  !   It is assumed the mortality is distributed in the same way over
  !   LOC (U1) and detritus (U6) s for phytoplankton and microzooplankton.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd  =  ( p_sd(bac)* et+( p_sd2(bac)* T1c(:)))* T1c(:)
  call flux_vector( iiIce, ppT1c,ppU6c, rd*( ONE- p_pe_R1c) )
  call flux_vector( iiIce, ppT1n,ppU6n, rd* qncSBA(bac, :)*( ONE- p_pe_R1n) )
  call flux_vector( iiIce, ppT1p,ppU6p, rd* qpcSBA(bac, :)*( ONE- p_pe_R1p) )

  call flux_vector( iiIce, ppT1c,ppU1c, rd* p_pe_R1c )
  call flux_vector( iiIce, ppT1n,ppU1n, rd* qncSBA(bac, :)* p_pe_R1n )
  call flux_vector( iiIce, ppT1p,ppU1p, rd* qpcSBA(bac, :)* p_pe_R1p )

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate quota in U1c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  qncSDE (iiU1, :) =   U1n(:)/ (p_small + U1c(:))
  qpcSDE (iiU1, :) =   U1p(:)/ (p_small + U1c(:))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Substrate availability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version)

    case ( 1 )  !LUCA

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rum  =   p_sum(bac)* et* T1c(:)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! No correction of food avilabilities:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      cuU1  =   min(  ONE, qpcSDE(iiU1,:)/ p_qpc(bac),  qncSDE(iiU1,:)/ p_qnc(bac))
      cuU6  =   ONE

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! oxygen environment:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      eO2  =   min(  ONE,  eO2mO2(:))


    case ( 2 )  ! BFM option

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
      ! Nutrient limitation (intracellular)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      iNIn  =   min(  ONE,  max(  ZERO,   qncSBA(bac, :)/ p_qnc(bac)))  !Nitrogen
      iI1p  =   min(  ONE,  max(  ZERO,   qpcSBA(bac, :)/ p_qpc(bac)))  !Phosphorus

      iN  =   min(  iI1p,  iNIn)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  =   p_sum(bac)* iN* et* T1c(:)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! correction of food avilabilities dependent on internal quota
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      cuU1  =   min(  ONE, qpcSDE(iiU1,:)/ p_qpc(bac),  qncSDE(iiU1,:)/ p_qnc(bac))
      cuU6  =   min(  ONE, qpcSDE(iiU6,:)/ p_qpc(bac),  qncSDE(iiU6,:)/ p_qnc(bac))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! oxygen environment:
      ! To provide a faster switching between the two metabolic pathways the
      ! oxygen dependence eO2 has been changed from the standard
      !     eO2 = MM(O2.o, p_chdo) to the cubic one written below.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      eO2  =   MM_power_vector(max(p_small,F2o(:)),  p_chdo(bac),3)

    case ( 3 )  ! Piets,BFM option

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
      ! Nutrient limitation (intracellular)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      iN  =   ONE

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  =   p_sum(bac)* iN* et* T1c(:)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! correction of food avilabilities dependent on internal quota
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      cuU1  =   ONE
      cuU6  =   ONE

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! oxygen environment:
      ! To provide a faster switching between the two metabolic pathways the
      ! oxygen dependence eO2 has been changed from the standard
      !     eO2 = MM(O2.o, p_chdo) to the cubic one written below.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      eO2  =   MM_power_vector(max(p_small,F2o(:)),  p_chdo(bac),3)

  end select


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate amount for U1, U6, and U2 and total amount of substrate avilable
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  ruU1c  =   (p_suhU1(bac)* cuU1(:) + p_sulU1(bac)*(ONE-cuU1(:))) * U1c(:)
  ruU6c  =   p_suU6(bac)* cuU6* U6c(:)
  !ruU2c  =   p_suU2* U2c(:)
  !rut  =   p_small + ruU6c+ ruU2c+ ruU1c
   rut  =   p_small + ruU6c+ ruU1c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rug  =   min(  rum,  rut)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruU6c  =   rug* ruU6c/ rut
  !ruU2c  =   rug* ruU2c/ rut
  ruU1c  =   rug* ruU1c/ rut

  call flux_vector( iiIce, ppU6c,ppT1c, ruU6c )
  !call flux_vector( iiIce, ppU2c,ppT1c, ruU2c )
  call flux_vector( iiIce, ppU1c,ppT1c, ruU1c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruU6n  =   qncSDE(iiU6,:)* ruU6c
  ruU1n  =   qncSDE(iiU1,:)* ruU1c

  call flux_vector( iiIce, ppU6n,ppT1n, ruU6n )
  call flux_vector( iiIce, ppU1n,ppT1n, ruU1n )

  ruU6p  =   qpcSDE(iiU6,:)* ruU6c
  ruU1p  =   qpcSDE(iiU1,:)* ruU1c

  call flux_vector( iiIce, ppU6p,ppT1p, ruU6p )
  call flux_vector( iiIce, ppU1p,ppT1p, ruU1p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Aerobic Respiration calculation + flux
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrc  =  ( p_pu_ra(bac)+ p_pu_ra_o(bac)*( ONE- eO2))* rug+ p_srs(bac)* T1c(:)* et
  call flux_vector( iiIce, ppT1c,ppF3c, rrc )
  call flux_vector( iiIce, ppF2o,ppF2o,-rrc/ MW_C) 

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  run  =   rug- rrc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !             Fluxes from bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version)
    case ( 1 )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Only uptake of ammonium possible
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ren  =  ( qncSBA(bac, :)- p_qnc(bac))* T1c(:)* ONE_PER_DAY
      call flux_vector( iiIce, ppT1n,ppI4n, ren* insw_vector( ren) )
      call flux_vector(iiIce, ppI4n,ppT1n,- ren* insw_vector( - ren)* I4n(:)/( &
        ONE+ I4n(:)))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rep  =  ( qpcSBA(bac, :)- p_qpc(bac))* T1c(:)* ONE_PER_DAY
      call flux_vector( iiIce, ppT1p,ppI1p,  rep* insw_vector( rep) )
      call flux_vector( iiIce, ppI1p,ppT1p,- rep* insw_vector( - rep)* I1p(:)/( &
        0.5+ I1p(:)))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Activity exrecetion (defined as reU7c) + stress excetion (defined as &
      ! reUc)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !reU7c  =   p_pu_ea_R7* run

      r  =   max(  ONE- qpcSBA(bac, :)/ p_qpc(bac),  ONE- qncSBA(bac, :)/ p_qnc(bac))
      !reU2c  =   ONE_PER_DAY* r* insw_vector(  r)* T1c(:)

      !run  =   run- reU7c- reU2c



    case ( 2,3 )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Nutrient uptake
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cqun3  =  p_lI4(bac)/( p_lI4(bac)+ I4n(:))
      rumn3  =  max(ZERO,p_qun(bac)* I3n(:)* T1c(:)* cqun3)  ! max pot. uptake of N3
      rumn4  =  max(ZERO,p_qun(bac)* I4n(:)* T1c(:))  ! max pot. uptake of N4
      rumn  =   rumn3+ rumn4
     ! misn  =   run/T1c(:)*( p_qnc(bac)* T1c(:)- T1n(:))  ! intracellular missing amount of N
      misn  =   run/(T1c(:)+p_small)*( p_qnc(bac)* T1c(:)- T1n(:))  ! intracellular missing amount of N
      rupn  =   run* p_qnc(bac)  ! N uptake based on C uptake
      runn=     min(rumn,rupn+misn)

      rump  =   max(ZERO,p_qup(bac)* I1p(:)* T1c(:))  ! max pot. uptake
      !misp  =   run/T1c(:)*( p_qpc(bac)* T1c(:)- T1p(:))  ! intracellular missing amount of P
      misp  =   run/(T1c(:)+p_small)*( p_qpc(bac)* T1c(:)- T1p(:))  ! intracellular missing amount of P
      rupp  =   run* p_qpc(bac)  ! P uptake based on C uptake
      !!LETI limitation for phosphorous
      !rupp  =   run* p_qpc(bac)* I1p(:)/(I1p(:)+0.025_RLEN)
      runp=     min(rump,rupp+misp)



      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Only stress excetion (defined as reU7c) , no other excretion (reU2c=0)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r      =  min(  run, ( ruU6n+ ruU1n+ runn)/ p_qlnc(bac))
      !reU7c  =  run- min(  r, ( ruU6p+ ruU1p+ runp)/ p_qlpc)
      !reU7c  =  max(  ZERO,  reU7c)

      !reU2c  =  ZERO
      !run    =  run- reU7c

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! insw: No excretion if net. growth <0
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      ren  =   max(ruU6n+ruU1n-run*p_qnc(bac),-runn) *insw_vector(run)
      ! excess of nutrients : ren > 0
      call flux_vector( iiIce, ppT1n,ppI4n,  ren*insw_vector(ren) )

      ! shortage of nutrients : ren < 0 --> Nutrient uptake
      r=-ren*insw_vector(-ren)
      call flux_vector(iiIce, ppI4n,ppT1n, r* rumn4/( p_small+ rumn))
      call flux_vector(iiIce, ppI3n,ppT1n, r* rumn3/( p_small+ rumn))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! insw: No excretion if net. growth <0
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      rep  =   max(ruU6p+ruU1p-run*p_qpc(bac),-runp) *insw_vector(run)
      
      ! excess of nutrients : rep > 0
      call flux_vector( iiIce, ppT1p,ppI1p, rep* insw_vector(rep) )

      ! shortage of nutrients : rep < 0 --> Nutrient uptake
      call flux_vector( iiIce, ppI1p,ppT1p, -rep*insw_vector(-rep) )

  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion fluxes + correction net prod.:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !call flux_vector( iiIce, ppT1c,ppU2c, reU2c )
  !call flux_vector( iiIce, ppT1c,ppU7c, reU7c )

  end subroutine SeaiceBacDynamics

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
