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
!    groups. The differences in behaviour
!    are expressed by differences in parameter-values only.
!    
! !INTERFACE
  subroutine PhytoDynamics(phyto)
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
  use mem, ONLY: iiC,iiN,iiP,iiS,iiL
  use mem, ONLY: D3STATE, R1c, R6c, O2o, R2c, &
                 N3n, N4n, N1p, R1n, R6n, R1p, R6p, N5s
  use mem, ONLY: ppR1c, ppR6c, ppO2o, ppO3c, ppR2c, ppN3n, ppN4n, ppN1p, ppR1n, &
    ppR6n, ppR1p, ppR6p, ppN5s, ppR6s, SUNQ, ThereIsLight, ETW, EIR, &
    xEPS, Depth, eiPI, sediPI, sunPI, qpPc, qnPc, qsPc, qlPc, NO_BOXES, &
    iiBen, iiPel, flux_vector, sourcesink_flux_vector
  use mem, ONLY: ppPhytoPlankton
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF,N7f,qfPc,ppN7f,ppR6f,ppR1f
#endif
#endif
  use constants,  ONLY: SEC_PER_DAY, E2W, HOURS_PER_DAY
  use mem_Param,  ONLY: p_small, ChlLightFlag, ProductionLightFlag, p_qchlc, &
                        LightLocationFlag, ChlSynthesisFlag
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
  integer                         :: silica_control
  integer, save :: first=0
  integer :: ppphytoc, ppphyton, ppphytop, ppphytos, ppphytol 
  real(RLEN),allocatable,save,dimension(:) :: phytoc,phyton,phytop,phytos,phytol
                                                                                                                                                             
  real(RLEN),allocatable,save,dimension(:) :: r,tmp,et,sum,sadap,sea,sdo,rugc,sra,srs, &
                                       srt,slc,run,pe_R6,rupp,rump,misp,rupn, &
                                       rumn3,rumn4,rumn,netgrowth,misn,cqun3
  real(RLEN),allocatable,save,dimension(:) :: rums,rups,miss,tN,fpplim,iN,iN1p,iNIn,eN5s,rrc,rr1c, &
                                       rr1n,rr1p,rr6c,rr6n,rr6p,rr6s,runn,runn3, &
                                       runn4,runp,runs,Irr,rho_Chl,rate_Chl,seo,flPIR2c
  real(RLEN),allocatable,save,dimension(:) :: iN5s,chl_opt
#ifdef INCLUDE_PELFE
  integer :: ppphytof
  real(RLEN),allocatable,save,dimension(:) :: phytof
  real(RLEN),allocatable,save,dimension(:) :: iN7f,misf,rr1f,rr6f,rupf,rumf,runf
#endif
  integer :: AllocStatus, DeallocStatus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if (first==0) then
     first=1
     allocate(phytoc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating phytoc"
     allocate(phyton(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating phyton"
     allocate(phytop(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating phytop"
     allocate(phytos(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating phytos"
     allocate(phytol(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating phytol"
     allocate(r(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating r"
     allocate(tmp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tmp"
     allocate(et(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating et"
     allocate(sum(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sum"
     allocate(sadap(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sadap"
     allocate(sea(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sea"
     allocate(sdo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sdo"
     allocate(rugc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rugc"
     allocate(sra(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sra"
     allocate(srs(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating srs"
     allocate(srt(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating srt"
     allocate(slc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating slc"
     allocate(run(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating run"
     allocate(pe_R6(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pe_R6"
     allocate(rupp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rupp"
     allocate(rump(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rump"
     allocate(misp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating misp"
     allocate(rupn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rupn"
     allocate(rumn3(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumn3"
     allocate(rumn4(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumn4"
     allocate(rumn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumn"
     allocate(netgrowth(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating netgrowth"
     allocate(misn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating misn"
     allocate(cqun3(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating cqun3"
     allocate(rums(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rums"
     allocate(rups(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rups"
     allocate(miss(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating miss"
     allocate(tN(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tN"
     allocate(fpplim(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fpplim"
     allocate(iN(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iN"
     allocate(iN1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iN1p"
     allocate(iNIn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iNIn"
     allocate(eN5s(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eN5s"
     allocate(rrc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrc"
     allocate(rr1c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr1c"
     allocate(rr1n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr1n"
     allocate(rr1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr1p"
     allocate(rr6c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6c"
     allocate(rr6n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6n"
     allocate(rr6p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6p"
     allocate(rr6s(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6s"
     allocate(runn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runn"
     allocate(runn3(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runn3"
     allocate(runn4(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runn4"
     allocate(runp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runp"
     allocate(runs(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runs"
     allocate(Irr(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating Irr"
     allocate(rho_Chl(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rho_Chl"
     allocate(rate_Chl(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rate_Chl"
     allocate(flPIR2c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating flPIR2c"
     allocate(seo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating seo"
     allocate(iN5s(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iN5s"
     allocate(chl_opt(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating chl_opt"
#ifdef INCLUDE_PELFE
     allocate(phytof(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating phytof"
     allocate(misf(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating misf"
     allocate(iN7f(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iN7f"
     allocate(rr1f(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr1f"
     allocate(rr6f(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rr6f"
     allocate(rupf(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rupf"
     allocate(rumf(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumf"
     allocate(runf(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runf"
#endif
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppphytoc = ppPhytoPlankton(phyto,iiC)
  ppphyton = ppPhytoPlankton(phyto,iiN)
  ppphytop = ppPhytoPlankton(phyto,iiP)
  ppphytos = ppPhytoPlankton(phyto,iiS)
  ppphytol = ppPhytoPlankton(phyto,iiL)
  phytoc(:) = D3STATE(ppphytoc,:)
  phyton(:) = D3STATE(ppphyton,:)
  phytop(:) = D3STATE(ppphytop,:)
  phytol(:) = D3STATE(ppphytol,:)
  if ( ppphytos > 0 )  phytos(:) = D3STATE(ppphytos,:)
#ifdef INCLUDE_PELFE
  ppphytof = ppPhytoPlankton(phyto,iiF)
  phytof(:) = D3STATE(ppphytof,:)
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  silica_control =1 : external regulation of silica limitation 
  !  silica_control =2 : internal regulation of silica limitation 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   silica_control=1
   if ( p_qus(phyto) > ZERO )  silica_control=2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitations (intracellular and extracellular)
  ! fpplim is the combined non-dimensional factor limiting photosynthesis
  ! Note for silicate limitation:
  ! The standard Michaelis-Menten formulation contains the Contois parameter
  ! p_Contois=0: standard Michaelis Menten Formulation
  ! 0<p_Contois<=1: The Contois formulation is active. 
  !                 The limiting role of the population size (intraspecific 
  !                 competition) can be tuned by increasing p_Contois 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  iN1p = min( ONE, max( p_small, ( qpPc(phyto,:) &
         - p_qplc(phyto))/( p_qpRc(phyto)- p_qplc(phyto))))
  iNIn = min( ONE, max( p_small, ( qnPc(phyto,:) &
         - p_qnlc(phyto))/( p_qnRc(phyto)- p_qnlc(phyto))))
  iN5s = min(ONE, max( p_small, ( qsPc(phyto,:) &
         - p_qslc(phyto))/( p_qsRc(phyto)- p_qslc(phyto))))
  eN5s = min( ONE,N5s(:)/(N5s(:) + p_chPs(phyto)));
  eN5s = min( ONE, N5s(:)/(N5s(:) + p_chPs(phyto)+(p_Contois(phyto)*phytos(:))))
  select case (silica_control) 
    case (1)  ! external control
      fpplim = eN5s
    case (2) ! internal control
      fpplim = iN5s
  end select
#ifdef INCLUDE_PELFE
  iN7f = min( ONE, max( p_small, ( qfPc(phyto,:) &
         - p_qflc(phyto))/( p_qfRc(phyto)- p_qflc(phyto))))
  fpplim = fpplim*iN7f
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Multiple nutrient limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_limnut(phyto))
    case ( 0 )
      iN  =   (iN1p* iNIn)**(0.5_RLEN)  ! geometric mean

    case ( 1 )
      iN  =   min(  iN1p,  iNIn)  ! Liebig rule

    case ( 2 )
      iN  =   2.0_RLEN/( ONE/ iN1p+ ONE/ iNIn)  ! combined

  end select

  ! tN only controls sedimentation of phytoplankton (Liebig)
  tN= min(iN,fpplim)

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
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute exponent E_PAR/E_K = alpha0/PBmax
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  r(:) = qlPc(phyto, :)*p_alpha_chl(phyto)/p_sum(phyto)* Irr

  select case ( ProductionLightFlag)
    case ( 1 ) ! instantaneous light
      ! no other factors needed
    case ( 2 ) ! daylight average is used
      ! recompute r and photsynthesis limitation using daylight scaling
      fpplim  =   fpplim*SUNQ/HOURS_PER_DAY
      r(:) = r(:)*HOURS_PER_DAY/SUNQ
    case ( 3 ) ! on-off
      fpplim  =   fpplim*ThereIsLight
  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Light limitation factor and total photosynthesis
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eiPI(phyto,:) = ( ONE- exp( - r))
  sum  =   p_sum(phyto)*et*eiPI(phyto,:)*fpplim

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sdo  =  ( p_thdo(phyto)/( iN+ p_thdo(phyto)))* p_sdmo(phyto)  ! nutr. -stress lysis
  ! extra lysis for high-density
  sdo  =   sdo+ p_seo(phyto)* MM_vector(phytoc, p_sheo(phyto))

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
  sra  =   p_pu_ra(phyto)*( sum - sea - seo)  ! activity
  srs  =   et* p_srs(phyto)                   ! basal
  srt  =   sra+ srs                           ! total
  rrc  =   srt* phytoc                        ! total actual respiration

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
     rr1c = rr1c + p_switchR1R2*sea*phytoc
     ! Nutrient-stress excretion is assigned to R2
     flPIR2c  =  seo*phytoc + (ONE-p_switchR1R2)*sea*phytoc
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
     sadap  =   max(  0.05_RLEN,  sum- slc)
  else
     sadap  =   et*p_sum(phyto)
  end if
  run  =   max(  ZERO, ( sum- slc)* phytoc)  ! net production

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Uptake: calculate maximum uptake of N, P
  ! based on affinity
  ! Ammonium preference is considered if p_lN4 /= 0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (p_netgrowth(phyto)) then
     cqun3  =   p_lN4(phyto)/( p_lN4(phyto)+ N4n(:))
  else
     cqun3 = p_lN4(phyto)
  end if
  rumn3  =   p_qun(phyto)* N3n(:)* phytoc* cqun3  ! max pot. uptake of N3
  rumn4  =   p_qun(phyto)* N4n(:)* phytoc  ! max pot. uptake of N4
  rumn  =   rumn3+ rumn4  ! max pot. uptake of DIN

  rump  =   p_qup(phyto)* N1p(:)* phytoc  ! max pot. uptake of PO4

  if (p_netgrowth(phyto)) then
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Check which fraction of fixed C can be used for new biomass
   ! given the internal storage.
   ! N and P uptake are compared sequentially
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth = min( run, ( rumn+ max( ZERO, 0.05_RLEN* &
      rugc*( qnPc(phyto, :)- p_qnlc(phyto))))/ p_qnlc(phyto))
      netgrowth = min( netgrowth, ( rump+ max( ZERO, &
       0.05_RLEN* rugc*( qpPc(phyto, :)- p_qplc(phyto))))/ p_qplc(phyto))
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Excrete C that cannot be used for growth as carbo-hydrates:
   ! Correct the net C uptake
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      netgrowth  =   max(  netgrowth,  ZERO)
      flPIR2c  =   flPIR2c+ run- netgrowth
      run  =   netgrowth
  end if

  call flux_vector( iiPel, ppphytoc,ppR2c, flPIR2c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Specific net growth rate (d-1)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sunPI(phyto,:)  =   run/( p_small+ phytoc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misn  =   sadap*( p_xqn(phyto)* p_qnRc(phyto)* phytoc- phyton)  ! Intracellular missing amount of N
  rupn  =   p_xqn(phyto)* p_qnRc(phyto)* run  ! N uptake based on net assimilat. C
#ifdef EXTRACOST
  rupn  =   p_xqn(phyto)* p_qnRc(phyto)* run-( srs+ sdo)* phyton  ! N uptake based on net assimilat. C
#endif
  runn  =   min(  rumn,  rupn+ misn)  ! actual uptake of NI

  r  =   insw_vector(  runn)
  runn3  =   r* runn* rumn3/( p_small+ rumn)  ! actual uptake of Nn
  runn4  =   r* runn* rumn4/( p_small+ rumn)  ! actual uptake of Nn
  call flux_vector( iiPel, ppN3n,ppphyton, runn3 )  ! source/sink.n
  call flux_vector( iiPel, ppN4n,ppphyton, runn4 )  ! source/sink.n
  tmp = - runn*( ONE- r)
  call flux_vector(iiPel, ppphyton,ppR1n,tmp)  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  misp  =   sadap*( p_xqp(phyto)* p_qpRc(phyto)* phytoc- phytop)  ! intracellular missing amount of P
  rupp  =   p_xqp(phyto)* run* p_qpRc(phyto)  ! P uptake based on C uptake
#ifdef EXTRACOST
  rupp  =   p_xqp(phyto)* run* p_qpRc(phyto)-( sdo+ srs)* phytop  ! P uptake based on C uptake
#endif
  runp  =   min(  rump,  rupp+ misp)  ! actual uptake

  r  =   insw_vector(  runp)
  tmp = runp*r
  call flux_vector( iiPel, ppN1p,ppphytop, tmp )  ! source/sink.p
  tmp = - runp*( ONE- r)
  call flux_vector(iiPel, ppphytop,ppR1p, tmp)  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     rr6n  =  pe_R6* sdo* phyton
     rr1n  =  sdo* phyton- rr6n

     rr6p  =  pe_R6* sdo* phytop
     rr1p  =  sdo* phytop- rr6p

  call flux_vector( iiPel, ppphyton,ppR6n, rr6n )  ! source/sink.n
  call flux_vector( iiPel, ppphyton,ppR1n, rr1n )  ! source/sink.n

  call flux_vector( iiPel, ppphytop,ppR6p, rr6p )  ! source/sink.p
  call flux_vector( iiPel, ppphytop,ppR1p, rr1p )  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( ppphytos > 0 )  then
    select case (silica_control)
    case (1)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Gross uptake of silicate excluding respiratory costs
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      runs = max(ZERO, p_qsRc(phyto) * (sum-srs) * phytoc)
    case (2)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Silicate uptake based on intracellular needs (note, no luxury)
      !  There can be efflux of dissolved silicate (M-J et al., 2000)
      !  however this generates fake remineralization and it is not implemented
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rums  =   p_qus(phyto)* N5s(:)* phytoc  ! max pot uptake based on affinity
      miss  =   max(ZERO, p_qsRc(phyto)*phytoc - phytos) ! intracellular missing Si
      rups  =   run* p_qsRc(phyto)* phytos  ! Si uptake based on net C uptake
      runs  =   min(  rums,  rups+ miss)  ! actual uptake
    end select
              
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Uptake and Losses of Si (only lysis)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    call flux_vector( iiPel, ppN5s,ppphytos, runs)
    call flux_vector( iiPel, ppphytos, ppR6s, sdo*phytos )
  endif

#ifdef INCLUDE_PELFE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: IRON
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Nutrient uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rumf  =   p_quf(phyto)* N7f(:)* phytoc  ! max potential uptake for iron
  misf  =   sadap*max(ZERO,p_xqf(phyto)*p_qfRc(phyto)*phytoc - phytof)  ! intracellular missing amount of F
  rupf  =   p_xqp(phyto)* run* p_qfRc(phyto)  ! Fe uptake based on C uptake
  runf  =   min(  rumf,  rupf+ misf)  ! actual uptake

  r  =   insw_vector(  runf)
  call flux_vector( iiPel, ppN7f,ppphytof, runf* r )  ! source/sink.p
  call flux_vector(iiPel, ppphytof,ppR1f,- runf*( ONE- r))  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Losses of Fe
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rr6f  =   rr6c* p_qflc(phyto)
  rr1f  =   sdo* phytof- rr6f

  call flux_vector( iiPel, ppphytof,ppR1f, rr1f )  ! source/sink.fe
  call flux_vector( iiPel, ppphytof,ppR6f, rr6f )  ! source/sink.fe
#endif

  if ( ChlLightFlag== 2) then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Chl-a synthesis and photoacclimation
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case (ChlSynthesisFlag)
      case (1) ! PELAGOS
           rho_Chl = p_qchlc( phyto)* min(ONE, p_sum(phyto)* eiPI(phyto,:)* phytoc/( &
                     p_alpha_chl(phyto)*( phytol+ p_small)* Irr))
           rate_Chl = rho_Chl*(sum - seo - sea - sra) * phytoc - sdo*phytol
      case (2) ! OPATM-BFM
           rho_Chl  =   p_qchlc(phyto)* sum/( p_alpha_chl(phyto)* qlPc(phyto,:)* Irr)
           rate_Chl = iN* rho_Chl* run- max( p_sdchl(phyto)*( ONE - iN), sdo)* &
               phytol+ min( ZERO, sum- slc+ sdo)* max( ZERO, phytol- p_qchlc(phyto)* phytoc)
      case (3) ! UNIBO
           rho_Chl = p_qchlc(phyto)*min(ONE,          &
                     (sum-seo-sea-sra) *phytoc /          &
                     (p_alpha_chl(phyto)*(phytol+p_small) *Irr))
           ! The "optimal" chl concentration corresponds to the chl that
           ! (given the actual C biomass) would give (Epar/Ek)=p_EpEk
           chl_opt = p_EpEk_or(phyto)*p_sum(phyto)*phytoc/  &
                     (p_alpha_chl(phyto)*Irr+p_small)
           !  Actual chlorophyll concentration exceeding the "optimal" value is 
           !  discarded with a p_tochl_relt relaxation.
           rate_Chl = rho_Chl*(sum-seo-sea-sra)*phytoc-(sdo*srs)*phytol - &
                      max(ZERO,(phytol-chl_opt))*p_tochl_relt(phyto)
      case (4) ! NIOZ
          ! total synthesis, only when there is net production (run > 0)
          ! The fixed loss rate due to basal respiration is introduced to have 
          ! mortality in the absence of light (< 1 uE/m2/s)
           rho_Chl = p_qchlc( phyto)* min(ONE, p_sum(phyto)* eiPI(phyto,:)* phytoc/( &
                     p_alpha_chl(phyto)*( phytol+ p_small)* Irr))
           rate_Chl = rho_Chl*run - p_sdchl(phyto)*phytol*max( ZERO, ( p_esNI(phyto)-tN)) &
                     -srs * phytol * ONE/(Irr+ONE)
    end select
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

  end subroutine PhytoDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
