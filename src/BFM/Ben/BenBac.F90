#include "DEBUG.h"
#include "INCLUDE.h"
#ifdef INCLUDE_BEN
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenBac
!
! DESCRIPTION
!   !    This submodel describes the carbon dynamics and associated
!    nutrient dynamics in benthic bacteria (represented
!    by state variables H1-H2) 
! 
! !INTERFACE
  subroutine BenBacDynamics(hx,  pphxc, pphxn, pphxp)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN, ZERO, ONE
  use constants, ONLY: MW_C
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D2STATE_BEN, Q6c, Q6n, Q6p, G2o, K16r, D6m, &
    D7m, D8m, D1m, D2m, BenDetritus, BenthicAmmonium, BenthicPhosphate
  use mem, ONLY: ppQ6c, ppQ6n, ppQ6p, ppG3c,ppG13c, ppG2o, ppK16r, ppD6m, &
    ppD7m, ppD8m, ppD1m, ppD2m, ppBenDetritus, ppBenthicAmmonium, &
    ppBenthicPhosphate, rrBTo, reBTn, reBTp, rrATo, reATn, reATp, ETW_Ben, ruHI, &
    iiH1, iiH2, iiC, iiN, iiP, NO_BOXES_XY_BEN, iiBen, iiPel, flux_vector, reHI
#endif
  use mem_Param,  ONLY: p_d_tot, p_small, p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro, &
                        p_clD1D2m
  use mem_BenBac


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, &
  ! PartQ_vector, eramp_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun, ONLY: eTq_vector, MM_vector, PartQ_vector, eramp_vector, &
    insw_vector


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: hx
  integer,intent(IN) :: pphxc
  integer,intent(IN) :: pphxn
  integer,intent(IN) :: pphxp

!  
!
! !AUTHORS
!   W. Ebenhoh and C. Kohlmeier,  
! 
!
!
! !REVISION_HISTORY
!   
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY_BEN) :: hxc
  real(RLEN),dimension(NO_BOXES_XY_BEN) :: hxn
  real(RLEN),dimension(NO_BOXES_XY_BEN) :: hxp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: clm
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: cm
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: chm
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: cmm
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: et
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: eN
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: eo
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: availQ6_c
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: availQ6_p
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: availQ6_n
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: suQ1
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: ruQ1c
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: ruQ1n
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: ruQ1p
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: ruQ6c
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: ruQ6n
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: ruQ6p
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rumn
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rump
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: sm
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: misn
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: misp
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: ren
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rep
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: sHc
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: sHn
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: sHp
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rq6c
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rqt6c
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rqt6n
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rqt6p
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: qnQ6c
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: qpQ6c
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rut
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rum
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: run
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rug
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: netgrowth
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: r
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: runp
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: runn
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rupp
  real(RLEN),dimension(NO_BOXES_XY_BEN)  :: rupn
  integer,dimension(NO_BOXES_XY_BEN)  :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  hxc = D2STATE_BEN(pphxc,:)
  hxn = D2STATE_BEN(pphxn,:)
  hxp = D2STATE_BEN(pphxp,:)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assign functional group-dependent parameters:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( hx)

    case ( iiH1 )
      clm  =   ZERO
      cm  =   D1m(:)
      chm  =   D1m(:)

    case ( iiH2 )
      clm  =   D1m(:)
      cm  =   D2m(:)
      chm  =   p_d_tot

  end select

  ! Determine where the Q6 center of mass is in the range from clm to D1m
  cmm  =   clm-log(0.5_RLEN*(RLEN+exp(-(chm-clm)/D6m(:))))*D6m(:)*0.5_RLEN


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW_Ben(:),  p_q10(hx))
  eo  =   MM_vector(  max(cm-clm,p_clD1D2m), p_cdm(hx))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus (if eaten): calculate available amount
  ! and add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  availQ6_c  =   Q6c(:)* PartQ_vector(  D6m(:),  clm,  chm,  p_d_tot)
  availQ6_n  =   Q6n(:)* PartQ_vector(  D7m(:),  clm,  chm,  p_d_tot)
  availQ6_p  =   Q6p(:)* PartQ_vector(  D8m(:),  clm,  chm,  p_d_tot)

  qnQ6c  =   availQ6_n/( availQ6_c+ 1.0D-30)
  qpQ6c  =   availQ6_p/( availQ6_c+ 1.0D-30)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Growth is controlled by quality of detritus (N and P content):
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  eN  =   eramp_vector(qnQ6c, p_qnc(hx))* eramp_vector(qpQ6c, p_qpc(hx))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total substrate availability:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6c  =  p_suhQ6(hx)* availQ6_c* eN+ p_sulQ6(hx)* availQ6_c
  ruQ1c  =  p_suQ1(hx)* BenDetritus(p_iQ1(hx),iiC)

  rut  =   ruQ6c+ ruQ1c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rum  =   p_sum(hx)* et* eo* hxc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rug  = min(rum, rut)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6c = rug* ruQ6c/ rut
  ruQ1c = rug* ruQ1c/ rut

  call flux_vector( iiBen, ppQ6c,pphxc, ruQ6c )
  call flux_vector( iiBen, ppBenDetritus(p_iQ1(hx),iiC),pphxc, ruQ1c )
  ruHI(hx,:)  =   ruQ1c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient fluxes into bacteria from carbon fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6n  =   ruQ6c* qnQ6c
  ruQ6p  =   ruQ6c* qpQ6c
  ruQ1n  =   p_suQ1(hx)* BenDetritus(p_iQ1(hx),iiN)* rug/ rut
  ruQ1p  =   p_suQ1(hx)* BenDetritus(p_iQ1(hx),iiP)* rug/ rut

  call flux_vector( iiBen, ppQ6n,pphxn, ruQ6n )
  call flux_vector( iiBen, ppQ6p,pphxp, ruQ6p )

  call flux_vector( iiBen, ppBenDetritus(p_iQ1(hx),iiN),pphxn, ruQ1n )
  call flux_vector( iiBen, ppBenDetritus(p_iQ1(hx),iiP),pphxp, ruQ1p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrc  =   p_srr(hx)* hxc* et+( ruQ1c+ ruQ6c)* p_pur(hx)

  run  =   max(ZERO, rug- rrc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of potential nutrient uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rumn =   BenthicAmmonium(p_iK4(hx),iiN) * p_sumKIn(hx) *hxc
  rump =   BenthicPhosphate(p_iK1(hx),iiP)* p_sumKIp(hx) *hxc
   
  ! a too high nutrient content of food is related to the food upake ( r < 0) 
  ! a too low nutrient content of food is related to net growth of the organisms   ( r > 0)
  r=p_qnc(hx)* hxc(:)-hxn(:)
  misn = (run*r*insw_vector(r)+max(rrc,rug)*r*insw_vector(-r))/hxc(:)
  r=p_qpc(hx)* hxc(:)-hxp(:)
  misp = (run*r*insw_vector(r)+max(rrc,rug)*r*insw_vector(-r))/hxc(:)
   
  rupn=run*p_qnc(hx)
  rupp=run*p_qpc(hx)

  runn=min(rumn-min(ZERO,misn),rupn+max(ZERO,misn))
  runp=min(rump-min(ZERO,misp),rupp+max(ZERO,misp))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon correction: all C which cannot be used for growth due to
  ! lack of nutrients is excreted!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  netgrowth  =   min(  run,       ( ruQ6n+ ruQ1n+ runn)/ p_qlnc(hx))
  netgrowth  =   min(  netgrowth, ( ruQ6p+ ruQ1p+ runp)/ p_qlpc(hx))
  netgrowth  =   max(  netgrowth,  ZERO)

  call flux_vector( iiBen, pphxc,pphxc,-( run- netgrowth) )
  run  =   netgrowth

  ren  =   max(ruQ6n+ruQ1n-run*p_qnc(hx),-runn) *insw_vector(run) -min(ZERO,misn)
  rep  =   max(ruQ6p+ruQ1p-run*p_qpc(hx),-runp) *insw_vector(run) -min(ZERO,misp)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  sm  =   p_sd(hx)*( ONE- eo)

  rqt6c  =   hxc(:)* sm*( ONE- p_pe_R1c)
  rqt6n  =   hxn(:)* sm*( ONE- p_pe_R1n)
  rqt6p  =   hxp(:)* sm*( ONE- p_pe_R1p)

  call flux_vector( iiBen, pphxc,ppBenDetritus(p_iQ1(hx),iiC), hxc(:)* sm* p_pe_R1c )
  call flux_vector( iiBen, pphxn,ppBenDetritus(p_iQ1(hx),iiN), hxn(:)* sm* p_pe_R1n )
  call flux_vector( iiBen, pphxp,ppBenDetritus(p_iQ1(hx),iiP), hxp(:)* sm* p_pe_R1p )

  reHI(hx,:)  =  hxc(:)* sm* p_pe_R1c;

  call flux_vector( iiBen, pphxc,ppQ6c, rqt6c )
  call flux_vector( iiBen, pphxn,ppQ6n, rqt6n )
  call flux_vector( iiBen, pphxp,ppQ6p, rqt6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector(iiBen, pphxn,ppBenthicAmmonium(p_iK4(hx),iiN),  ren*insw_vector(ren))
  call flux_vector(iiBen, ppBenthicAmmonium(p_iK4(hx),iiN),pphxn, -ren*insw_vector(-ren))

  call flux_vector(iiBen, pphxp,ppBenthicPhosphate(p_iK1(hx),iiP),  rep* insw_vector(rep) )
  call flux_vector(iiBen, ppBenthicPhosphate(p_iK1(hx),iiP),pphxp, -rep* insw_vector(-rep))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes depends on type of bacteria and layers in which they occur:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( hx)

    case ( iiH1 )
      call flux_vector( iiBen, pphxc,ppG3c,rrc )
      call flux_vector(iiBen, ppG2o,ppG2o,-( rrc/ MW_C))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Add respiration and excretion to benthic totals:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rrBTo(:)  =   rrBTo(:)+ rrc/ MW_C
      reBTn(:)  =   reBTn(:)+ ren
      reBTp(:)  =   reBTp(:)+ rep

    case ( iiH2 )
      call flux_vector( iiBen, pphxc,ppG13c,rrc )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Respiration in anoxic circumstances produces reduced material
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux_vector( iiBen, ppK16r,ppK16r, rrc/ MW_C* p_qro )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Add respiration and excretion to benthic totals:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rrATo(:)  =   rrATo(:)+ rrc/ MW_C
      reATn(:)  =   reATn(:)+ ren
      reATp(:)  =   reATp(:)+ rep

  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to uptake of detritus in distribution
  ! of state variables (Dx.m is a undetermined source)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector(iiBen, ppD6m,ppD6m,( cmm- D6m(:))*( ruQ6c- rqt6c)/ Q6c(:))
  call flux_vector(iiBen, ppD7m,ppD7m,( cmm- D7m(:))*( ruQ6n- rqt6n)/ Q6n(:))
  call flux_vector(iiBen, ppD8m,ppD8m,( cmm- D8m(:))*( ruQ6p- rqt6p)/ Q6p(:))

  end subroutine BenBacDynamics
#endif
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
