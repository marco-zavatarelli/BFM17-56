#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   This process describes the dynamics of all phytoplankton
!    groups in the ERSEM model. The differences in behaviour
!    are expressed by differences in parameter-values only.
!    
! !INTERFACE
  subroutine PhytoDynamics(phyto, ppphytoc, ppphyton, ppphytop, ppphytos, &
    ppphytol)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants, ONLY: MW_C
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D3STATE, R1c, R6c, O2o, R2c, &
                 N3n, N4n, N1p, R1n, R6n, R1p, R6p, N5s
  use mem, ONLY: ppR1c, ppR6c, ppO2o, ppO3c, ppR2c, ppN3n, ppN4n, ppN1p, ppR1n, &
    ppR6n, ppR1p, ppR6p, ppN5s, SUNQ, ThereIsLight, flP1R6s, ETW, EIR, xEPS, &
    Depth, eiPI, sediPI, sunPI, qpPc, qnPc, qsPc, qlPc, iiP1, iiP4, NO_BOXES, &
    iiBen, iiPel, flux_vector, sourcesink_flux_vector
#endif
  use constants,  ONLY: SEC_PER_DAY, E2W, HOURS_PER_DAY
  use mem_Param,  ONLY: p_small, ChlLightFlag, LightForcingFlag, p_qchlc, &
                        LightLocationFlag
  use mem_Phyto


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
!   ERSEM group + J.G. Baretta-Bekker + W.Ebenhoeh
!     P. Ruardij (NIOZ)
!
!
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
  real(RLEN),dimension(NO_BOXES) :: phytoc
  real(RLEN),dimension(NO_BOXES) :: phyton
  real(RLEN),dimension(NO_BOXES) :: phytop
  real(RLEN),dimension(NO_BOXES) :: phytos
  real(RLEN),dimension(NO_BOXES) :: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                         :: silica_control
  real(RLEN),dimension(NO_BOXES)  :: r,tmp
  real(RLEN),dimension(NO_BOXES)  :: et
  real(RLEN),dimension(NO_BOXES)  :: sum
  real(RLEN),dimension(NO_BOXES)  :: sadap
  real(RLEN),dimension(NO_BOXES)  :: sea
  real(RLEN),dimension(NO_BOXES)  :: sdo
  real(RLEN),dimension(NO_BOXES)  :: rugc
  real(RLEN),dimension(NO_BOXES)  :: sra
  real(RLEN),dimension(NO_BOXES)  :: srs
  real(RLEN),dimension(NO_BOXES)  :: srt
  real(RLEN),dimension(NO_BOXES)  :: slc
  real(RLEN),dimension(NO_BOXES)  :: run
  real(RLEN),dimension(NO_BOXES)  :: pe_R6
  real(RLEN),dimension(NO_BOXES)  :: rupp
  real(RLEN),dimension(NO_BOXES)  :: rump
  real(RLEN),dimension(NO_BOXES)  :: misp
  real(RLEN),dimension(NO_BOXES)  :: rupn
  real(RLEN),dimension(NO_BOXES)  :: rumn3
  real(RLEN),dimension(NO_BOXES)  :: rumn4
  real(RLEN),dimension(NO_BOXES)  :: rumn
  real(RLEN),dimension(NO_BOXES)  :: netgrowth
  real(RLEN),dimension(NO_BOXES)  :: misn
  real(RLEN),dimension(NO_BOXES)  :: cqun3
  real(RLEN),dimension(NO_BOXES)  :: rums
  real(RLEN),dimension(NO_BOXES)  :: rups
  real(RLEN),dimension(NO_BOXES)  :: miss
  real(RLEN),dimension(NO_BOXES)  :: tN
  real(RLEN),dimension(NO_BOXES)  :: iN
  real(RLEN),dimension(NO_BOXES)  :: iN1p
  real(RLEN),dimension(NO_BOXES)  :: iNIn
  real(RLEN),dimension(NO_BOXES)  :: eN5s
  real(RLEN),dimension(NO_BOXES)  :: rrc
  real(RLEN),dimension(NO_BOXES)  :: rr1c
  real(RLEN),dimension(NO_BOXES)  :: rr1n
  real(RLEN),dimension(NO_BOXES)  :: rr1p
  real(RLEN),dimension(NO_BOXES)  :: rr6c
  real(RLEN),dimension(NO_BOXES)  :: rr6n
  real(RLEN),dimension(NO_BOXES)  :: rr6p
  real(RLEN),dimension(NO_BOXES)  :: rr6s
  real(RLEN),dimension(NO_BOXES)  :: runn
  real(RLEN),dimension(NO_BOXES)  :: runn3
  real(RLEN),dimension(NO_BOXES)  :: runn4
  real(RLEN),dimension(NO_BOXES)  :: runp
  real(RLEN),dimension(NO_BOXES)  :: runs
  real(RLEN),dimension(NO_BOXES)  :: Irr
  real(RLEN),dimension(NO_BOXES)  :: rho_Chl
  real(RLEN),dimension(NO_BOXES)  :: rate_Chl
  real(RLEN),dimension(NO_BOXES)  :: flPIR2c

  real(RLEN),dimension(NO_BOXES)  :: seo
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
  
   ! force external regulation with nutrient-stress excretion
   if ( (.not.p_netgrowth(phyto)).and.(ppphytos > 0))  silica_control=1
 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  phytoc = D3STATE(ppphytoc,:)
  phyton = D3STATE(ppphyton,:)
  phytop = D3STATE(ppphytop,:)
  phytol = D3STATE(ppphytol,:)
  if ( ppphytos > 0 )  phytos = D3STATE(ppphytos,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular) N, P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iN1p = min( ONE, max( p_small, ( qpPc(phyto,:) &
         - p_qplc(phyto))/( p_qpRc(phyto)- p_qplc(phyto))))
  iNIn = min( ONE, max( p_small, ( qnPc(phyto,:) &
         - p_qnlc(phyto))/( p_qnRc(phyto)- p_qnlc(phyto))))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Phytoplankton growth is limited by nitrogen and phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_limnut(phyto))
    case ( 0 )
      iN  =   (iN1p* iNIn)**(0.5D+00)  ! geometric mean

    case ( 1 )
      iN  =   min(  iN1p,  iNIn)  ! Liebig rule

    case ( 2 )
      iN  =   2.0D+00/( ONE/ iN1p+ ONE/ iNIn)  ! combined

  end select

  ! tN controls sedimentation of phytoplankton
  tN= iN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation due to intra- extracellular silicate.
  ! eN5s limit externally nutrient limitation.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !calculate internal quota
  if ( silica_control > 0 ) then
    select case (silica_control) 
      case(1)
        eN5s = min( ONE,N5s(:)/(N5s(:) + p_chPs(phyto)));
        tN=min(iN,eN5s);
      case(2)
        eN5s=ONE
        r = max( p_small, ( qsPc(phyto,:)- p_qslc(phyto))/ &
                                    (( p_qsRc(phyto)- p_qslc(phyto))))
        tN=min(r,iN);
      end select
  else
    eN5s  =   ONE
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Phytoplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector(  ETW(:),  p_q10(phyto))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Photosynthesis 
  ! Irradiance EIR is in uE m-2 s-1, 
  ! Irr is top, middle or average irradiance in uE m-2 day-1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( ChlLightFlag== 2) then

    select case ( LightLocationFlag)
       case ( 1 )
          ! Light at the top of the cell
          Irr =  max(p_small,EIR(:))*SEC_PER_DAY;
       case ( 2 )
          ! Light in the middle of the cell
          Irr = max(p_small,EIR(:))*exp(-xEPS(:)*0.5_RLEN*Depth(:)) &
                *SEC_PER_DAY
       case ( 3 ) ! default
          ! Average Light in the cell
          r = xEPS(:)* Depth(:)
          r = EIR(:)/xEPS(:)/Depth(:)*(ONE-exp(-r))
          Irr = max(p_small,r*SEC_PER_DAY)
    end select

    r(:) = qlPc(phyto, :)* p_alpha_chl(phyto)/p_sum(phyto)* Irr
    eiPI(phyto,:) = ( ONE- exp( - r))

  end if

  select case ( LightForcingFlag)
    case ( 1 )
      sum  =   p_sum(phyto)* et* eiPI(phyto,:)   *  eN5s

    case ( 2 )
      sum  =   p_sum(phyto)* et* eiPI(phyto,:)*( SUNQ/ HOURS_PER_DAY) * eN5s

    case ( 3 )
      sum  =   p_sum(phyto)* et* eiPI(phyto,:)* ThereIsLight * eN5s

  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sdo  =  ( p_thdo(phyto)/( iN+ p_thdo(phyto)))* p_sdmo(phyto)  ! nutr. -stress lysis
  !MAV: TO BE REMOVED ASAP !
  ! extra lysis for high-density
  sdo  =   sdo+ p_seo(phyto)* MM_vector(  phytoc,  100.0_RLEN)

  sea  =   sum* p_pu_ea(phyto)  ! activity excretion

  if (p_netgrowth(phyto)) then
     seo = ZERO
  else 
     ! nutrient stress excretion
     seo = sum*(ONE-p_pu_ea(phyto))*(ONE- iN) 
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apportioning over R1 and R6:
  ! Cell lysis generates both DOM and POM.
  ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
  ! Assuming that this structural part is not easily degradable,
  ! at least a fraction equal to the minimum quota is released as POM.
  ! Therefore, nutrients (and C) in the structural part go to R6.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pe_R6 = min( p_qplc(phyto)/( qpPc(phyto, :)+ p_small), p_qnlc(phyto)/ &
          ( qnPc(phyto, :)+ p_small))
  pe_R6  =   min(  ONE,  pe_R6)
  rr6c  =   pe_R6* sdo* phytoc
  rr1c  =  ( ONE- pe_R6)* sdo* phytoc

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
     flPIR2c  =   sea* phytoc
  else
     ! Activity excretion is assigned to R1
     rr1c = rr1c + sea*phytoc
     ! Nutrient-stress excretion is assigned to R2
     flPIR2c  =   seo*phytoc
  end if

  call sourcesink_flux_vector( iiPel, ppO3c,ppphytoc, rugc )  
  call flux_vector( iiPel, ppphytoc,ppR1c, rr1c )
  call flux_vector( iiPel, ppphytoc,ppR6c, rr6c )


  call sourcesink_flux_vector( iiPel, ppphytoc,ppO3c, rrc )
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrc/ MW_C) )
  call flux_vector( iiPel, ppO2o,ppO2o, rugc/ MW_C ) 


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

  cqun3  =   p_lN4(phyto)/( p_lN4(phyto)+ N4n(:))
  rumn3  =   p_qun(phyto)* N3n(:)* phytoc* cqun3  ! max pot. uptake of N3
  rumn4  =   p_qun(phyto)* N4n(:)* phytoc  ! max pot. uptake of N4
  rumn  =   rumn3+ rumn4  ! max pot. uptake of NI

  rump  =   p_qup(phyto)* N1p(:)* phytoc  ! max pot. uptake

  if (p_netgrowth(phyto)) then
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Check which fraction of fixed C can be used for new biomass
   ! by comparing the potential nutrient availability
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth = min( run, ( rumn+ max( ZERO, 0.05D+00* &
      rugc*( qnPc(phyto, :)- p_qnlc(phyto))))/ p_qnlc(phyto))
      netgrowth = min( netgrowth, ( rump+ max( ZERO, &
       0.05D+00* rugc*( qpPc(phyto, :)- p_qplc(phyto))))/ p_qplc(phyto))
      if ( silica_control  ==  2) then
          rums  =   p_qus(phyto)* N5s(:)* phytoc  ! max pot uptake
          netgrowth = min( netgrowth, ( rums+ max( ZERO, &
            1.00D+00* rugc*( qsPc(phyto, :)- p_qslc(phyto))))/ p_qslc(phyto))
      endif
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Excrete C which can not be used for growth as carbo-hydrates:
   ! Correct net C-uptake
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth  =   max(  netgrowth,  ZERO)
      flPIR2c  =   flPIR2c+ run- netgrowth
      run  =   netgrowth
  end if

  call flux_vector( iiPel, ppphytoc,ppR2c, flPIR2c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apparent Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  sunPI(phyto,:)  =   run/( p_small+ phytoc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misn  =   sadap*( p_xqn(phyto)* p_qnRc(phyto)* phytoc- phyton)  ! Intracellular missing amount of N
  rupn  =   p_xqn(phyto)* p_qnRc(phyto)* run-( srs+ sdo)* phyton  ! N uptake based on net assimilat. C
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of NI

  r  =   insw_vector(  runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of Nn
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of Nn
  call flux_vector( iiPel, ppN3n,ppphyton, runn3 )  ! source/sink.n
  call flux_vector( iiPel, ppN4n,ppphyton, runn4 )  ! source/sink.n
  tmp = - runn*( ONE- r)
  call flux_vector(iiPel, ppphyton,ppN4n,tmp)  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misp  =   sadap*( p_xqp(phyto)* p_qpRc(phyto)* phytoc- phytop)  ! intracellular missing amount of P
  rupp  =   p_xqp(phyto)* run* p_qpRc(phyto)-( sdo+ srs)* phytop  ! P uptake based on C uptake
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

  r  =   insw_vector(  runp)
  tmp = runp*r
  call flux_vector( iiPel, ppN1p,ppphytop, tmp )  ! source/sink.p
  tmp = - runp*( ONE- r)
  call flux_vector(iiPel, ppphytop,ppN1p, tmp)  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr6n  =   pe_R6* sdo* phyton
  rr1n  =   sdo* phyton- rr6n

  rr6p  =   pe_R6* sdo* phytop
  rr1p  =   sdo* phytop- rr6p

  call flux_vector( iiPel, ppphyton,ppR1n, rr1n )  ! source/sink.n
  call flux_vector( iiPel, ppphyton,ppR6n, rr6n )  ! source/sink.n

  call flux_vector( iiPel, ppphytop,ppR1p, rr1p )  ! source/sink.p
  call flux_vector( iiPel, ppphytop,ppR6p, rr6p )  ! source/sink.p


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( silica_control > 0 )  then
    select case (silica_control)

    case (1)

     ! Net uptake of silicate
     runs = max(ZERO, p_qsRc(phyto) * run )
     call flux_vector( iiPel, ppN5s,ppphytos, runs)
     ! The fixed loss rate for basal respiration is maintained to have 
     ! constant Si:C quota in the absence of production
     flP1R6s(:)  =   flP1R6s(:)+ srs*phytos

    case (2)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Nutrient uptake
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      miss  =   sadap*( p_xqs(phyto)* p_qsRc(phyto)* phytoc- phytos)  ! intracellular missing Si
      rups  =   (run* p_qsRc(phyto)-( sdo+ srs)* phytos)  ! Si uptake based on C uptake
      runs  =   min(  rums,  rups+ miss)  ! actual uptake

      call flux_vector( iiPel, ppN5s,ppphytos, runs* insw_vector(runs) )  ! source/sink.c
      call flux_vector(iiPel, ppphytos,ppN5s,- runs*insw_vector(-runs))  ! source/sink.c
    end select
              
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Losses of Si (lysis)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr6s  =   (sdo+srs) * phytos  ! Lysis, particulate

    ! Collect first all fluxes of P-->silica
    flP1R6s(:)  =   flP1R6s(:)+ rr6s
  endif


  if ( ChlLightFlag== 2) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Chl-a synthesis and photoacclimation
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rho_Chl = p_qchlc( phyto)* min(ONE, p_sum(phyto)* eiPI(phyto,:)* phytoc/( &
          p_alpha_chl(phyto)*( phytol+ p_small)* Irr))
    ! total synthesis, only when there is net production (run > 0)
    ! The fixed loss rate due to basal respiration is introduced to have 
    ! mortality in the absence of light (< 1 uE/m2/s)
    rate_Chl = rho_Chl*run - p_sdchl(phyto)*phytol*max( ZERO, ( p_esNI(phyto)-tN)) &
                  -srs * phytol * ONE/(Irr+ONE)

    call flux_vector( iiPel, ppphytol,ppphytol, rate_Chl )
  end if



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Sedimentation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( p_res(phyto)> ZERO) then
    sediPI(phyto,:) = sediPI(phyto,:) &
                   + p_res(phyto)* max( ZERO, ( p_esNI(phyto)-tN))
  end if


  ! End of computation section for process PhytoDynamics




  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
