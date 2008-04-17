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
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: D2STATE, U1c, U6c, F2o, F3c, &
                 I3n, I4n, I1p, U1n, U6n, U1p, U6p, I5s
#endif
  use mem, ONLY: ppU1c, ppU6c, ppF2o, ppF3c, ppI3n, ppI4n, ppI1p, ppU1n, &
    ppU6n, ppU1p, ppU6p, ppI5s, SUNQ, ThereIsLight, flP1R6s, ETB, EIB, &
    EHB, eiSI, iiS1, qnSc, qpSc, qsSc, qlSc, sediPI, sunPI, NO_BOXES_XY, &
    iiBen, flux_vector, sourcesink_flux_vector
  use constants,  ONLY: SEC_PER_DAY, E2W, HOURS_PER_DAY
  use mem_Param,  ONLY: p_small, ChlLightFlag, LightForcingFlag 
  use mem_Seaicealgae


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
  real(RLEN),dimension(NO_BOXES_XY) :: phytoc
  real(RLEN),dimension(NO_BOXES_XY) :: phyton
  real(RLEN),dimension(NO_BOXES_XY) :: phytop
  real(RLEN),dimension(NO_BOXES_XY) :: phytos
  real(RLEN),dimension(NO_BOXES_XY) :: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                         :: silica_control
  integer,dimension(NO_BOXES_XY)     :: i
  real(RLEN),dimension(NO_BOXES_XY)  :: r
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: sum
  real(RLEN),dimension(NO_BOXES_XY)  :: sadap
  real(RLEN),dimension(NO_BOXES_XY)  :: sea
  real(RLEN),dimension(NO_BOXES_XY)  :: sdo
  real(RLEN),dimension(NO_BOXES_XY)  :: rugc
  real(RLEN),dimension(NO_BOXES_XY)  :: sra
  real(RLEN),dimension(NO_BOXES_XY)  :: srs
  real(RLEN),dimension(NO_BOXES_XY)  :: srt
  real(RLEN),dimension(NO_BOXES_XY)  :: slc
  real(RLEN),dimension(NO_BOXES_XY)  :: run
  real(RLEN),dimension(NO_BOXES_XY)  :: pe_U6
  real(RLEN),dimension(NO_BOXES_XY)  :: rupp
  real(RLEN),dimension(NO_BOXES_XY)  :: rump
  real(RLEN),dimension(NO_BOXES_XY)  :: misp
  real(RLEN),dimension(NO_BOXES_XY)  :: rupn
  real(RLEN),dimension(NO_BOXES_XY)  :: rumn3
  real(RLEN),dimension(NO_BOXES_XY)  :: rumn4
  real(RLEN),dimension(NO_BOXES_XY)  :: rumn
  real(RLEN),dimension(NO_BOXES_XY)  :: netgrowth
  real(RLEN),dimension(NO_BOXES_XY)  :: misn
  real(RLEN),dimension(NO_BOXES_XY)  :: cqun3
  real(RLEN),dimension(NO_BOXES_XY)  :: rums
  real(RLEN),dimension(NO_BOXES_XY)  :: rups
  real(RLEN),dimension(NO_BOXES_XY)  :: miss
  real(RLEN),dimension(NO_BOXES_XY)  :: tI
  real(RLEN),dimension(NO_BOXES_XY)  :: iI
  real(RLEN),dimension(NO_BOXES_XY)  :: iI1p
  real(RLEN),dimension(NO_BOXES_XY)  :: iIIn
  real(RLEN),dimension(NO_BOXES_XY)  :: eI5s
  real(RLEN),dimension(NO_BOXES_XY)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY)  :: rr1c
  real(RLEN),dimension(NO_BOXES_XY)  :: rr1n
  real(RLEN),dimension(NO_BOXES_XY)  :: rr1p
  real(RLEN),dimension(NO_BOXES_XY)  :: rr6c
  real(RLEN),dimension(NO_BOXES_XY)  :: rr6n
  real(RLEN),dimension(NO_BOXES_XY)  :: rr6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rr6s
  real(RLEN),dimension(NO_BOXES_XY)  :: runn
  real(RLEN),dimension(NO_BOXES_XY)  :: runn3
  real(RLEN),dimension(NO_BOXES_XY)  :: runn4
  real(RLEN),dimension(NO_BOXES_XY)  :: runp
  real(RLEN),dimension(NO_BOXES_XY)  :: runs
  real(RLEN),dimension(NO_BOXES_XY)  :: Irr
  real(RLEN),dimension(NO_BOXES_XY)  :: rho_Chl
  real(RLEN),dimension(NO_BOXES_XY)  :: rate_Chl
  real(RLEN),dimension(NO_BOXES_XY)  :: Photo_max
  real(RLEN),dimension(NO_BOXES_XY)  :: flSIU2c

  real(RLEN),dimension(NO_BOXES_XY)  :: seo
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
  phytoc = D2STATE(ppphytoc,:)
  phyton = D2STATE(ppphyton,:)
  phytop = D2STATE(ppphytop,:)
  phytol = D2STATE(ppphytol,:)
  if ( ppphytos > 0 )  phytos = D2STATE(ppphytos,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular) N, P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iI1p = min( ONE, max( p_small, ( qpSc(phyto, &
    :)- p_qplc(phyto))/( p_qpRc(phyto)- p_qplc(phyto))))
  iIIn = min( ONE, max( p_small, ( qnSc(phyto, &
    :)- p_qnlc(phyto))/( p_qnRc(phyto)- p_qnlc(phyto))))

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


  ! tI controls sedimentation of phytoplankton
  tI= iI;

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation due to intra- extracellular silicate.
  ! eI5s limit externally nutrient limitation.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !calculate internal quota
  if ( silica_control > 0 ) then
    select case (silica_control) 
      case(1)
        eI5s = min( ONE,I5s/(I5s + p_chPs(phyto)));
        tI=min(iI,eI5s);
      case(2)
        eI5s=ONE
        r = max( p_small, ( qsSc(phyto,:)- p_qslc(phyto))/ &
                                    (( p_qsRc(phyto)- p_qslc(phyto))))
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
    ! Light is already at the middle of the cell
    Irr =  max(p_small, EIB(:))*SEC_PER_DAY;
    eiSI(phyto,:) = ( ONE- exp( - qlSc(phyto, :)* p_alpha_chl(phyto)/ &
      p_sum(phyto)* Irr))

!  end if

  select case ( LightForcingFlag)
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
!   if (p_netgrowth(phyto)) then
!      ! Activity excretion is assigned to R2
!      flPIR2c  =   sea* phytoc
!   else
!      ! Activity excretion is assigned to R1
!      rr1c = rr1c + sea*phytoc
!      ! Nutrient-stress excretion is assigned to R2
!      flPIR2c  =   seo*phytoc
!   end if

  !call flux_vector( iiBen, ppphytoc,ppphytoc, rugc )
  call sourcesink_flux_vector( iiBen,ppF3c,ppphytoc,rugc )
  call flux_vector( iiBen, ppphytoc,ppU1c, rr1c )
  call flux_vector( iiBen, ppphytoc,ppU6c, rr6c )


  !call flux_vector( iiBen, ppphytoc,ppphytoc, rrc )
  call sourcesink_flux_vector( iiBen, ppphytoc,ppF3c,rrc )
  call flux_vector( iiBen, ppF2o,ppF2o,-( rrc/ MW_C) )
  call flux_vector( iiBen, ppF2o,ppF2o, rugc/ MW_C )


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
!       flPIR2c  =   flPIR2c+ run- netgrowth
      run  =   netgrowth
  end if

!   call flux_vector( iiBen, ppphytoc, ppU2c, flPIR2c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misn  =   sadap*( p_xqn(phyto)* p_qnRc(phyto)* phytoc- phyton)  ! Intracellular missing amount of N
  rupn  =   p_xqn(phyto)* p_qnRc(phyto)* run-( srs+ sdo)* phyton  ! N uptake based on net assimilat. C
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of II

  r  =   insw_vector(  runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of In
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of In
  call flux_vector( iiBen, ppI3n,ppphyton, runn3 )  ! source/sink.n
  call flux_vector( iiBen, ppI4n,ppphyton, runn4 )  ! source/sink.n
  call flux_vector(iiBen, ppphyton,ppI4n,- runn*( ONE- r))  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misp  =   sadap*( p_xqp(phyto)* p_qpRc(phyto)* phytoc- phytop)  ! intracellular missing amount of P
  rupp  =   p_xqp(phyto)* run* p_qpRc(phyto)-( sdo+ srs)* phytop  ! P uptake based on C uptake
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

  r  =   insw_vector(  runp)
  call flux_vector( iiBen, ppI1p,ppphytop, runp* r )  ! source/sink.p
  call flux_vector(iiBen, ppphytop,ppI1p,- runp*( ONE- r))  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr6n  =   pe_U6* sdo* phyton
  rr1n  =   sdo* phyton- rr6n

  rr6p  =   pe_U6* sdo* phytop
  rr1p  =   sdo* phytop- rr6p

  call flux_vector( iiBen, ppphyton,ppU1n, rr1n )  ! source/sink.n
  call flux_vector( iiBen, ppphyton,ppU6n, rr6n )  ! source/sink.n

  call flux_vector( iiBen, ppphytop,ppU1p, rr1p )  ! source/sink.p
  call flux_vector( iiBen, ppphytop,ppU6p, rr6p )  ! source/sink.p


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( silica_control > 0 )  then
    select case (silica_control)

    case (1)

     runs = max(ZERO, p_qsRc(phyto) * run );          ! net uptake
     call flux_vector( iiBen, ppI5s,ppphytos, runs)  ! source/sink.c

    case (2)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Nutrient uptake
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      miss  =   sadap*( p_xqs(phyto)* p_qsRc(phyto)* phytoc- phytos)  ! intracellular missing Si
      rups  =   (run* p_qsRc(phyto)-( sdo+ srs)* phytos)  ! Si uptake based on C uptake
      runs  =   min(  rums,  rups+ miss)  ! actual uptake

      call flux_vector( iiBen, ppI5s,ppphytos, runs* insw_vector(runs) )  ! source/sink.c
      call flux_vector(iiBen, ppphytos,ppI5s,- runs*insw_vector(-runs))  ! source/sink.c
    end select
              
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Losses of Si
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr6s  =   sdo* phytos  ! Lysis, particulate

    ! Collect first all fluxes of P-->silica
    flP1R6s(:)  =   flP1R6s(:)+ rr6s
  endif


  if ( ChlLightFlag== 2) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Chl-a synthesis and photoacclimation
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rho_Chl = p_qchlcSI( phyto)* p_sum(phyto)* eiSI(phyto,:)* phytoc/( &
          p_alpha_chl(phyto)*( phytol+ p_small)* Irr)
! total synthesis, only when there is net production (run > 0)
!       rate_Chl = rho_Chl*( sum- sea- sra)* phytoc- sdo* phytol+ min( &
!         ZERO, sum- slc+ sdo)* max( ZERO, phytol- p_qchlcSI( phyto)* phytoc)
        rate_Chl = rho_Chl*( sum- sea- sra)* phytoc- sdo* phytol+ min( &
          -p_sdchl(phyto), sum- slc+ sdo)* max( ZERO, phytol- p_qchlcSI( phyto)* phytoc)
!     rate_Chl = rho_Chl*( max(srs,sum-slc) )* phytoc- sdo* phytol+ min( &
!         -srs -p_sdchl, sum- slc+ sdo)* max( ZERO, phytol- p_qchlcSI( phyto)* phytoc)

        rate_Chl = rho_Chl*run - p_sdchl(phyto)*phytol*max( ZERO, ( p_chlNI(phyto)-tI))

    call flux_vector( iiBen, ppphytol,ppphytol, rate_Chl )
  end if

  ! End of computation section for process PhytoDynamics


  end
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
