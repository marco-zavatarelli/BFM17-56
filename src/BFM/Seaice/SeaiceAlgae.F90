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
  subroutine SeaiceAlgaeDynamics(phyto, ppphytoc, ppphyton, ppphytop, ppphytos, &
    ppphytol)
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
  use mem, ONLY: D2STATE_ICE, U1c, U6c, F2o, F3c, &
                 I3n, I4n, I1p, U1n, U6n, U1p, U6p, I5s
#endif
  use mem, ONLY: ppU1c, ppU6c, ppF2o, ppF3c, ppI3n, ppI4n, ppI1p, ppU1n, &
    ppU6n, ppU1p, ppU6p, ppU6s, ppI5s, SUNQ, ThereIsLight, ETB, EIB, &
    EHB, eiSI, iiS1, qnSc, qpSc, qsSc, qlSc, sediPPY, sunPPY, NO_BOXES_XY_ICE, &
    iiIce, flux_vector
  use constants,  ONLY: SEC_PER_DAY, E2W, HOURS_PER_DAY
  use mem_Param,  ONLY: p_small, ChlDynamicsFlag, LightPeriodFlag 
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
  integer,intent(IN)  :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol

!  
!
! !AUTHORS
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2007 the BFM team
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
  real(RLEN),dimension(NO_BOXES_XY_ICE) :: phytoc
  real(RLEN),dimension(NO_BOXES_XY_ICE) :: phyton
  real(RLEN),dimension(NO_BOXES_XY_ICE) :: phytop
  real(RLEN),dimension(NO_BOXES_XY_ICE) :: phytos
  real(RLEN),dimension(NO_BOXES_XY_ICE) :: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                         :: silica_control
  integer,dimension(NO_BOXES_XY_ICE)     :: i
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: r
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: et
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: sum
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: sadap
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: sea
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: sdo
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rugc
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: sra
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: srs
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: srt
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: slc
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: run
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: pe_U6
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rupp
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rump
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: misp
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rupn
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rumn3
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rumn4
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rumn
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: netgrowth
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: misn
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: cqun3
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rums
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rups
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: miss
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: tI
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: iI
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: iI1p
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: iIIn
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: eI5s
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rr1c
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rr1n
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rr1p
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rr6c
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rr6n
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rr6p
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rr6s
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: runn
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: runn3
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: runn4
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: runp
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: runs
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: Irr
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rho_Chl
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: rate_Chl
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: Photo_max
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: flSIU2c,flS1U6s
  real(RLEN),dimension(NO_BOXES_XY_ICE)  :: seo
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  silica_control =0 : no silica component present in cell
  !  silica_control =1 : external regulation of silica limitation & limitation of
  !                      carbon fixation under silica depletion
  !  silica_control =2 : internal regulation of silica limitation & excretion
  !                      of fixed carbon under nutrient stress
  !                      Process description based on:
  !                      Growth physiology and fate of diatoms in the ocean: a review
  !                      G.Sarthou, K.R. Timmermans, S. Blain, & P. Treguer
  !                      JSR 53 (2005) 25-42
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   silica_control=0
   if ( p_qus(phyto) > 0.0 )  then
      silica_control=2
   elseif ( p_chPs(phyto) > 0.0 ) then
      silica_control=1
   endif
   flS1U6s = ZERO
  
   ! force external regulation with nutrient-stress excretion
   if ( (.not.p_netgrowth(phyto)).and.(ppphytos > 0))  silica_control=1
 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  phytoc = D2STATE_ICE(ppphytoc,:)
  phyton = D2STATE_ICE(ppphyton,:)
  phytop = D2STATE_ICE(ppphytop,:)
  phytol = D2STATE_ICE(ppphytol,:)
  if ( ppphytos > 0 )  phytos = D2STATE_ICE(ppphytos,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular) N, P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iI1p = min( ONE, max( p_small, ( qpSc(phyto, &
    :)- p_qplc(phyto))/( p_qpcPPY(phyto)- p_qplc(phyto))))
  iIIn = min( ONE, max( p_small, ( qnSc(phyto, &
    :)- p_qnlc(phyto))/( p_qncPPY(phyto)- p_qnlc(phyto))))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Sea ice algae growth is limited by nitrogen and phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_limnut(phyto))
    case ( 0 )
      iI  =   (iI1p* iIIn)**(0.5D+00)  ! geometric mean

    case ( 1 )
      iI  =   min(  iI1p,  iIIn)  ! Liebig rule

    case ( 2 )
      iI  =   2.0D+00/( ONE/ iI1p+ ONE/ iIIn)  ! combined
  end select


  ! tI controls sedimentation of sea ice algae
  tI= iI;

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation due to intra- extracellular silicate.
  ! eI5s limit externally nutrient limitation.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !calculate internal quota
  if ( silica_control > 0 ) then
    select case (silica_control) 
      case(1)
        eI5s = min( ONE,I5s/(I5s + p_chPs(phyto))); ! Michaelis-Menten quota
        tI=min(iI,eI5s);
      case(2)
        eI5s=ONE
        r = max( p_small, ( qsSc(phyto,:)- p_qslc(phyto))/ &
                                    (( p_qscPPY(phyto)- p_qslc(phyto))))
        tI=min(r,iI);
      end select
  else
    eI5s  =   ONE
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Sea ice algae
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETB(:),  p_q10(phyto))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Photosynthesis (Irradiance EIB is in uE m-2 s-1)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Light is already at the middle of the cell in the BAL
  Irr =  max(p_small, EIB(:))*SEC_PER_DAY;
  eiSI(phyto,:) = ( ONE- exp( - qlSc(phyto, :)* p_alpha_chl(phyto)/ &
      p_sum(phyto)* Irr))

  select case ( LightPeriodFlag)
    case ( 1 )
      sum  =   p_sum(phyto)* et* eiSI(phyto,:)   *  eI5s

    case ( 2 )
      sum  =   p_sum(phyto)* et* eiSI(phyto,:)*( SUNQ/ HOURS_PER_DAY) * eI5s

    case ( 3 )
      sum  =   p_sum(phyto)* et* eiSI(phyto,:)* ThereIsLight * eI5s

  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sdo  =  ( p_thdo(phyto)/( iI+ p_thdo(phyto)))* p_sdmo(phyto)  ! nutr. -stress lysis
  sea  =   sum* p_pu_ea(phyto)  ! activity excretion

  if (p_netgrowth(phyto)) then
     seo = ZERO
  else 
     ! nutrient stress excretion
     seo = sum*(ONE-p_pu_ea(phyto))*(ONE- iI) 
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apportioning over R1 and R6:
  ! Cell lysis generates both DOM and POM.
  ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
  ! Assuming that this structural part is not easily degradable,
  ! at least a fraction equal to the minimum quota is released as POM.
  ! Therefore, nutrients (and C) in the structural part go to R6.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pe_U6 = min( p_qplc(phyto)/( qpSc(phyto, :)+ p_small), p_qnlc(phyto)/( &
    qnSc(phyto, :)+ p_small))
  pe_U6  =   min(  ONE,  pe_U6)
  rr6c  =   pe_U6* sdo* phytoc
  rr1c  =  ( ONE- pe_U6)* sdo* phytoc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration rate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sra  =   p_pu_ra(phyto)*( sum- sea)  ! activity
  srs  =   et* p_srs(phyto)  ! rest
  srt  =   sra+ srs  ! total
  rrc  =   srt* phytoc  ! total actual respiration

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Production, productivity and C flows
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  =   sum* phytoc  ! gross production
  slc  =   sea + seo + srt+ sdo  ! specific loss terms
  if (p_netgrowth(phyto)) then
      ! Activity excretion is assigned to R2
      flSIU2c  =   sea* phytoc
   else
      ! Activity excretion is assigned to R1
      rr1c = rr1c + sea*phytoc
      ! Nutrient-stress excretion is assigned to R2
      flSIU2c  =   seo*phytoc
   end if

  !call flux_vector( iiIce, ppphytoc,ppphytoc, rugc )
  call flux_vector( iiIce,ppF3c,ppphytoc,rugc )
  call flux_vector( iiIce, ppphytoc,ppU6c, rr6c )
  call flux_vector( iiIce, ppphytoc,ppU1c, rr1c )

  !call flux_vector( iiIce, ppphytoc,ppphytoc, rrc )
  call flux_vector( iiIce, ppphytoc,ppF3c,rrc )
  call flux_vector( iiIce, ppF2o,ppF2o,-( rrc/ MW_C) )
  call flux_vector( iiIce, ppF2o,ppF2o, rugc/ MW_C )


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential-Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (p_netgrowth(phyto)) then
     sadap  =   max(  srs,  sum)
  else
     sadap  =   max(  srs,  p_sum(phyto))
  end if
  run  =   max(  ZERO, ( sum- slc)* phytoc)  ! net production

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Uptake: calculate maximal uptake of N, P
  ! Check if C-fixation is # larger to make of all C new biomass
  ! Assumed is that Si-depletion directly the growth rate in contradiction
  ! to N and P.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cqun3  =   p_lN4(phyto)/( p_lN4(phyto)+ I4n(:))
  rumn3  =   p_qun(phyto)* I3n(:)* phytoc* cqun3  ! max pot. uptake of I3
  rumn4  =   p_qun(phyto)* I4n(:)* phytoc  ! max pot. uptake of I4
  rumn  =   rumn3+ rumn4  ! max pot. uptake of II

  rump  =   p_qup(phyto)* I1p(:)* phytoc  ! max pot. uptake

  if (p_netgrowth(phyto)) then
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Check which fraction of fixed C can be used for new biomass
   ! by comparing the potential nutrient availability
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth = min( run, ( rumn+ max( ZERO, 0.05D+00* &
      rugc*( qnSc(phyto, :)- p_qnlc(phyto))))/ p_qnlc(phyto))
      netgrowth = min( netgrowth, ( rump+ max( ZERO, &
       0.05D+00* rugc*( qpSc(phyto, :)- p_qplc(phyto))))/ p_qplc(phyto))
      if ( silica_control  ==  2) then
          rums  =   p_qus(phyto)* I5s(:)* phytoc  ! max pot uptake
          netgrowth = min( netgrowth, ( rums+ max( ZERO, &
            1.00D+00* rugc*( qsSc(phyto, :)- p_qslc(phyto))))/ p_qslc(phyto))
      endif
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Excrete C which can not be used for growth as carbo-hydrates:
   ! Correct net C-uptake
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth  =   max(  netgrowth,  ZERO)
      flSIU2c  =   flSIU2c+ run- netgrowth
      run  =   netgrowth
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! fluxes to dissolved organic carbon
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Polysaccharides U2 are not defined in sea ice yet
  ! hence both rr1c and flSIU2c are assigned to U1
  !call flux_vector( iiIce, ppphytoc, ppU2c, flSIU2c )
  call flux_vector( iiIce, ppphytoc,ppU1c, rr1c+flSIU2c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misn  =   sadap*( p_xqn(phyto)* p_qncPPY(phyto)* phytoc- phyton)  ! Intracellular missing amount of N
  rupn  =   p_xqn(phyto)* p_qncPPY(phyto)* run-( srs+ sdo)* phyton  ! N uptake based on net assimilat. C
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of II

  r  =   insw_vector(  runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of In
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of In
  call flux_vector( iiIce, ppI3n,ppphyton, runn3 )  ! source/sink.n
  call flux_vector( iiIce, ppI4n,ppphyton, runn4 )  ! source/sink.n
  call flux_vector(iiIce, ppphyton,ppI4n,- runn*( ONE- r))  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misp  =   sadap*( p_xqp(phyto)* p_qpcPPY(phyto)* phytoc- phytop)  ! intracellular missing amount of P
  rupp  =   p_xqp(phyto)* run* p_qpcPPY(phyto)-( sdo+ srs)* phytop  ! P uptake based on C uptake
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

  r  =   insw_vector(  runp)
  call flux_vector( iiIce, ppI1p,ppphytop, runp* r )  ! source/sink.p
  call flux_vector(iiIce, ppphytop,ppI1p,- runp*( ONE- r))  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr6n  =   pe_U6* sdo* phyton
  rr1n  =   sdo* phyton- rr6n

  rr6p  =   pe_U6* sdo* phytop
  rr1p  =   sdo* phytop- rr6p

  call flux_vector( iiIce, ppphyton,ppU1n, rr1n )  ! source/sink.n
  call flux_vector( iiIce, ppphyton,ppU6n, rr6n )  ! source/sink.n

  call flux_vector( iiIce, ppphytop,ppU1p, rr1p )  ! source/sink.p
  call flux_vector( iiIce, ppphytop,ppU6p, rr6p )  ! source/sink.p


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( silica_control > 0 )  then
    select case (silica_control)

    case (1)

     runs = max(ZERO, p_qscPPY(phyto) * run );          ! net uptake
     call flux_vector( iiIce, ppI5s,ppphytos, runs)  ! source/sink.c
     ! The fixed loss rate for basal respiration is maintained to have 
     ! constant Si:C quota in the absence of production
     flS1U6s(:)  =   flS1U6s(:)+ srs*phytos
    case (2)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Nutrient uptake
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      miss  =   sadap*( p_xqs(phyto)* p_qscPPY(phyto)* phytoc- phytos)  ! intracellular missing Si
      rups  =   (run* p_qscPPY(phyto)-( sdo+ srs)* phytos)  ! Si uptake based on C uptake
      runs  =   min(  rums,  rups+ miss)  ! actual uptake

      call flux_vector( iiIce, ppI5s,ppphytos, runs* insw_vector(runs) )  ! source/sink.c
      call flux_vector(iiIce, ppphytos,ppI5s,- runs*insw_vector(-runs))  ! source/sink.c
    end select
              
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Losses of Si
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr6s  =   sdo* phytos  ! Lysis, particulate

    ! Collect first all fluxes of P-->silica
    flS1U6s(:)  =   flS1U6s(:)+ rr6s

  ! fluxes to biogenic particulate silica are assigned here.
  ! In the pelagic this is done in PelChem since the predation of zooplankton
  ! also release biogenic particulate silica
  ! Since there is no zooplankton here, we need to assign the flux
  !call flux_vector( iiPel, ppphytos,ppU6s, flS1U6s(:) )
  call flux_vector( iiIce, ppphytos,ppU6s, flS1U6s(:) )
 endif

  if ( ChlDynamicsFlag== 2) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Chl-a synthesis and photoacclimation
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rho_Chl = p_qlcPPYSI( phyto)* min(ONE, p_sum(phyto)* eiSI(phyto,:)* phytoc/( &
          p_alpha_chl(phyto)*( phytol+ p_small)* Irr))
      rate_Chl = rho_Chl*(sum - sea + seo) * phytoc - (srt+ sdo)*phytol
    call flux_vector( iiIce, ppphytol,ppphytol, rate_Chl )
  end if

  ! End of computation section for process PhytoDynamics


  end
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
