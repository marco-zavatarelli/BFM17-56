#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Seaicealgae
!
! DESCRIPTION
!   This process describes the dynamics of all sea ice algae
!    groups. The differences in behaviour
!    are expressed by differences in parameter-values only.
!    
! !INTERFACE
  subroutine SeaiceAlgaeDynamics(alg)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants, ONLY: MW_C
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE_ICE
#else
  use mem, ONLY: iiC,iiN,iiP,iiS,iiL
  use mem, ONLY: D2STATE_ICE, I3n, I4n, I1p, I5s
  use mem, ONLY: ppU1c, ppU6c, ppF2o, ppF3c, ppI3n, ppI4n, ppI1p, ppU1n, &
    ppU6n, ppU1p, ppU6p, ppU6s, ppI5s, ETB, EIB, &
    eiSAL, iiS1, qncSAL, qpcSAL, qlcSAL, NO_BOXES_ICE, &
    iiIce, flux_vector,ppSeaiceAlgae
#endif
  use constants,  ONLY: SEC_PER_DAY, E2W, HOURS_PER_DAY
  use mem_Param,  ONLY: p_small
  use mem_SeaiceAlgae


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector, insw_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: alg

!  
!
! !AUTHORS
! L. Tedesco and M. Vichi
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer :: ppalgc, ppalgn, ppalgp, ppalgs, ppalgl
  real(RLEN),dimension(NO_BOXES_ICE) :: algc
  real(RLEN),dimension(NO_BOXES_ICE) :: algn
  real(RLEN),dimension(NO_BOXES_ICE) :: algp
  real(RLEN),dimension(NO_BOXES_ICE) :: algs
  real(RLEN),dimension(NO_BOXES_ICE) :: algl
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_ICE)  :: r
  real(RLEN),dimension(NO_BOXES_ICE)  :: et
  real(RLEN),dimension(NO_BOXES_ICE)  :: sum
  real(RLEN),dimension(NO_BOXES_ICE)  :: sadap
  real(RLEN),dimension(NO_BOXES_ICE)  :: sea
  real(RLEN),dimension(NO_BOXES_ICE)  :: sdo
  real(RLEN),dimension(NO_BOXES_ICE)  :: rugc
  real(RLEN),dimension(NO_BOXES_ICE)  :: sra
  real(RLEN),dimension(NO_BOXES_ICE)  :: srs
  real(RLEN),dimension(NO_BOXES_ICE)  :: srt
  real(RLEN),dimension(NO_BOXES_ICE)  :: slc
  real(RLEN),dimension(NO_BOXES_ICE)  :: run
  real(RLEN),dimension(NO_BOXES_ICE)  :: pe_U6
  real(RLEN),dimension(NO_BOXES_ICE)  :: rupp
  real(RLEN),dimension(NO_BOXES_ICE)  :: rump
  real(RLEN),dimension(NO_BOXES_ICE)  :: misp
  real(RLEN),dimension(NO_BOXES_ICE)  :: rupn
  real(RLEN),dimension(NO_BOXES_ICE)  :: rumn3
  real(RLEN),dimension(NO_BOXES_ICE)  :: rumn4
  real(RLEN),dimension(NO_BOXES_ICE)  :: rumn
  real(RLEN),dimension(NO_BOXES_ICE)  :: misn
  real(RLEN),dimension(NO_BOXES_ICE)  :: cqun3
  real(RLEN),dimension(NO_BOXES_ICE)  :: iI
  real(RLEN),dimension(NO_BOXES_ICE)  :: iI1p
  real(RLEN),dimension(NO_BOXES_ICE)  :: iINn
  real(RLEN),dimension(NO_BOXES_ICE)  :: eI5s
  real(RLEN),dimension(NO_BOXES_ICE)  :: rrc
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr1c
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr1n
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr1p
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr6c
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr6n
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr6p
  real(RLEN),dimension(NO_BOXES_ICE)  :: rr6s
  real(RLEN),dimension(NO_BOXES_ICE)  :: runn
  real(RLEN),dimension(NO_BOXES_ICE)  :: runn3
  real(RLEN),dimension(NO_BOXES_ICE)  :: runn4
  real(RLEN),dimension(NO_BOXES_ICE)  :: runp
  real(RLEN),dimension(NO_BOXES_ICE)  :: runs
  real(RLEN),dimension(NO_BOXES_ICE)  :: Irr
  real(RLEN),dimension(NO_BOXES_ICE)  :: rho_Chl
  real(RLEN),dimension(NO_BOXES_ICE)  :: rate_Chl
  real(RLEN),dimension(NO_BOXES_ICE)  :: flSIU2c,flS1U6s
  real(RLEN),dimension(NO_BOXES_ICE)  :: seo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Reset flux to biogenic silica
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  flS1U6s = ZERO
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppalgc = ppSeaiceAlgae(alg,iiC)
  ppalgn = ppSeaiceAlgae(alg,iiN)
  ppalgp = ppSeaiceAlgae(alg,iiP)
  ppalgs = ppSeaiceAlgae(alg,iiS)
  ppalgl = ppSeaiceAlgae(alg,iiL)
  algc = D2STATE_ICE(ppalgc,:)
  algn = D2STATE_ICE(ppalgn,:)
  algp = D2STATE_ICE(ppalgp,:)
  algl = D2STATE_ICE(ppalgl,:)
  if ( ppalgs > 0 )  algs = D2STATE_ICE(ppalgs,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular) N, P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iI1p = min( ONE, max( p_small, ( qpcSAL(alg, &
    :)- p_qplc(alg))/( p_qpcSAL(alg)- p_qplc(alg))))
  iINn = min( ONE, max( p_small, ( qncSAL(alg, &
    :)- p_qnlc(alg))/( p_qncSAL(alg)- p_qnlc(alg))))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Sea ice algae growth is limited by nitrogen and phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_limnut(alg))
    case ( 0 )
      iI  =   (iI1p* iINn)**(0.5_RLEN)  ! geometric mean
    case ( 1 )
      iI  =   min(  iI1p,  iINn)  ! Liebig rule
    case ( 2 )
      iI  =   2.0_RLEN/( ONE/ iI1p+ ONE/ iINn)  ! combined
  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation due to extracellular silicate.
  ! eI5s controls externally nutrient limitation.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( ppalgs > 0 ) then
     eI5s = min( ONE,I5s/(I5s + p_chsSAL(alg))); ! Michaelis-Menten quota
  else
    eI5s  =   ONE
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Sea ice algae
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETB(:),  p_q10(alg))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Photosynthesis (Irradiance EIB is in uE m-2 s-1)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Light is already at the middle of the cell in the BAL
  Irr =  max(p_small, EIB(:))*SEC_PER_DAY;
  eiSAL(alg,:) = ( ONE- exp( - qlcSAL(alg, :)* p_alpha_chl(alg)/ &
                 p_sum(alg)* Irr))
  sum  =   p_sum(alg)* et* eiSAL(alg,:) * eI5s

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sdo  =  ( p_thdo(alg)/( iI+ p_thdo(alg)))* p_sdmo(alg)  ! nutr. -stress lysis
  sea  =   sum* p_pu_ea(alg)  ! activity excretion
  ! nutrient stress excretion
  seo = sum*(ONE-p_pu_ea(alg))*(ONE- iI) 

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apportioning over R1 and R6:
  ! Cell lysis generates both DOM and POM.
  ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
  ! Assuming that this structural part is not easily degradable,
  ! at least a fraction equal to the minimum quota is released as POM.
  ! Therefore, nutrients (and C) in the structural part go to R6.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pe_U6 = min( p_qplc(alg)/( qpcSAL(alg, :)+ p_small), p_qnlc(alg)/( &
             qncSAL(alg, :)+ p_small))
  pe_U6  =   min(  ONE,  pe_U6)
  rr6c  =   pe_U6* sdo* algc
  rr1c  =  ( ONE- pe_U6)* sdo* algc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration rate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sra  =   p_pu_ra(alg)*( sum- sea)  ! activity
  srs  =   et* p_srs(alg)  ! rest
  srt  =   sra+ srs  ! total
  rrc  =   srt* algc  ! total actual respiration

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Production, productivity and C flows
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  =   sum* algc  ! gross production
  slc  =   sea + seo + srt+ sdo  ! specific loss terms
  ! Activity excretion is assigned to U1
  rr1c = rr1c + sea*algc
  ! Nutrient-stress excretion is assigned to a temporary U2
  flSIU2c  =   seo*algc

  call flux_vector( iiIce,ppF3c,ppalgc,rugc )
  call flux_vector( iiIce, ppalgc,ppU6c, rr6c )
  call flux_vector( iiIce, ppalgc,ppU1c, rr1c )

  call flux_vector( iiIce, ppalgc,ppF3c,rrc )
  call flux_vector( iiIce, ppF2o,ppF2o,-(rrc/MW_C) )
  call flux_vector( iiIce, ppF2o,ppF2o, rugc/MW_C )


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential-Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sadap  =   max(  srs,  p_sum(alg))
  run  =   max(  ZERO, ( sum- slc)* algc)  ! net production

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Uptake: calculate maximal uptake of N, P
  ! Check if C-fixation is # larger to make of all C new biomass
  ! Assumed is that Si-depletion directly the growth rate in contradiction
  ! to N and P.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cqun3  =   p_lN4(alg)/( p_lN4(alg)+ I4n(:))
  rumn3  =   p_qun(alg)* I3n(:)* algc* cqun3  ! max pot. uptake of I3
  rumn4  =   p_qun(alg)* I4n(:)* algc  ! max pot. uptake of I4
  rumn  =   rumn3+ rumn4  ! max pot. uptake of SAL

  rump  =   p_qup(alg)* I1p(:)* algc  ! max pot. uptake

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! fluxes to dissolved organic carbon
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Polysaccharides U2 are not defined in sea ice yet
  ! hence both rr1c and flSIU2c are assigned to U1
  call flux_vector( iiIce, ppalgc,ppU1c, rr1c+flSIU2c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misn  =   sadap*( p_xqn(alg)* p_qncSAL(alg)* algc- algn)  ! Intracellular missing amount of N
  rupn  =   p_xqn(alg)* p_qncSAL(alg)* run-( srs+ sdo)* algn  ! N uptake based on net assimilat. C
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of SAL

  r  =   insw_vector(  runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of In
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of In
  call flux_vector( iiIce, ppI3n,ppalgn, runn3 )  ! source/sink.n
  call flux_vector( iiIce, ppI4n,ppalgn, runn4 )  ! source/sink.n
  call flux_vector(iiIce, ppalgn,ppI4n,- runn*( ONE- r))  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misp  =   sadap*( p_xqp(alg)* p_qpcSAL(alg)* algc- algp)  ! intracellular missing amount of P
  rupp  =   p_xqp(alg)* run* p_qpcSAL(alg)-( sdo+ srs)* algp  ! P uptake based on C uptake
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

  r  =   insw_vector(  runp)
  call flux_vector( iiIce, ppI1p,ppalgp, runp* r )  ! source/sink.p
  call flux_vector(iiIce, ppalgp,ppI1p,- runp*( ONE- r))  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr6n  =   pe_U6* sdo* algn
  rr1n  =   sdo* algn- rr6n

  rr6p  =   pe_U6* sdo* algp
  rr1p  =   sdo* algp- rr6p

  call flux_vector( iiIce, ppalgn,ppU1n, rr1n )  ! source/sink.n
  call flux_vector( iiIce, ppalgn,ppU6n, rr6n )  ! source/sink.n

  call flux_vector( iiIce, ppalgp,ppU1p, rr1p )  ! source/sink.p
  call flux_vector( iiIce, ppalgp,ppU6p, rr6p )  ! source/sink.p


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( ppalgs > 0 )  then
     runs = max(ZERO, p_qscSAL(alg) * run );          ! net uptake
     call flux_vector( iiIce, ppI5s,ppalgs, runs)  ! source/sink.c
     ! The fixed loss rate for basal respiration is maintained to have 
     ! constant Si:C quota in the absence of production
     flS1U6s(:)  =   flS1U6s(:)+ srs*algs
              
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Losses of Si
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr6s  =   sdo* algs  ! Lysis, particulate

    ! Collect first all fluxes of P-->silica
    flS1U6s(:)  =   flS1U6s(:)+ rr6s

!MAV: check this
    ! fluxes to biogenic particulate silica are assigned here.
    ! In the pelagic this is done in PelChem since the predation of zooplankton
    ! also release biogenic particulate silica
    ! Since there is no zooplankton here, we need to assign the flux
    call flux_vector( iiIce, ppalgs,ppU6s, flS1U6s(:) )
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Chl-a synthesis and photoacclimation
  ! Note that differently from phytoplankton, chl in sea ice algae 
  ! must be a variable
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rho_Chl = p_qlcSAL( alg)* min(ONE, p_sum(alg)* eiSAL(alg,:) *  &
            algc/(p_alpha_chl(alg)*( algl+ p_small)* Irr))
  rate_Chl = rho_Chl*(sum - sea + seo) * algc - (srt+ sdo)*algl
  call flux_vector( iiIce, ppalgl,ppalgl, rate_Chl )

  ! End of computation section for process SeaiceAlgaeDynamics
  end subroutine SeaiceAlgaeDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
