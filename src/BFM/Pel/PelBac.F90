#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelBac
!
! DESCRIPTION
!   This process describes the dynamics of bacterioplankton
!    
!
! !INTERFACE
  subroutine PelBacDynamics(bac)
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY:  R6c, R6n, R6p, R1c, R1n, R1p, R2c, O2o, N6r, &
    N4n, N1p, N3n, R3c, iiR1, iiR6, D3STATE
  use mem, ONLY: iiPelBacteria, ppPelBacteria, iiC, iiN, iiP, ppR6c, &
    ppR6n, ppR6p, ppR1c, ppR1n, ppR1p, &
    ppR2c, ppO2o, ppN6r, ppN4n, ppN1p, ppN3n, ppR3c, flPTN6r, Depth, ETW, &
    qncPBA, qpcPBA, eO2mO2, qpcOMT, qncOMT, NO_BOXES, iiBen, iiPel, flux_vector
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c
#endif
#endif
  use constants,  ONLY: MW_C, ONE_PER_DAY
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, p_qro, p_small
  use mem_PelBac
  use mem_globalfun,   ONLY: eTq_vector, MM_power_vector, insw_vector, &
                             MM_vector
!  
!
! !AUTHORS
!   First ERSEM version by J.W. Baretta and H. Baretta-Bekker
!   Additional parametrizations by P. Ruardij and M. Vichi 
!   (Vichi et al., 2004; Vichi et al., 2007)
!   L. Polimene, I. Allen and M. Zavatarelli (Polimene et al., 2006)
!   Dynamical allocation by G. Mattia 
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
! !INPUT:
  integer,intent(IN)  :: bac

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer       :: i
  integer       :: ppbacc, ppbacn, ppbacp
  integer, save :: first =0
  real(RLEN),allocatable,save,dimension(:) :: runn,runp,et,eO2,r,flN6rPBA,rrc,  &
                                          rd,ruR1c,ruR1n,ruR1p,ruR2c,ruR3c,  &
                                          ruR6c,ruR6p,ruR6n,cqun3,rump,  &
                                          rumn,rumn3,rumn4,ren,rep,reR2c, &
                                          reR3c,rut,rum,run,sun,rug,suR1, &
                                          suR1n,suR1p,suR2,cuR6,cuR1,iN1p, &
                                          iNIn,iN,eN1p,eN4n, &
                                          huln, hulp, bacc
  real(RLEN),allocatable,save,dimension(:) ::  misn,misp,rupp,rupn
  integer :: AllocStatus
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Local memory allocation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (first==0) then
     first=1
     allocate(misn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating misn"
     allocate(misp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating misp"
     allocate(rupp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rupp"
     allocate(rupn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rupn"
     allocate(runn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runn"
     allocate(runp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating runp"
     allocate(et(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating et"
     allocate(eO2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eO2"
     allocate(r(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating r"
     allocate(flN6rPBA(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating flN6rPBA"
     allocate(rrc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrc"
     allocate(rd(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rd"
     allocate(ruR1c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR1c"
     allocate(ruR1n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR1n"
     allocate(ruR1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR1p"
     allocate(ruR2c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR2c"
     allocate(ruR6c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR6c"
     allocate(ruR6p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR6p"
     allocate(ruR6n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR6n"
     allocate(ruR3c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruR3c"
     allocate(cqun3(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating cqun3"
     allocate(rump(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rump"
     allocate(rumn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumn"
     allocate(rumn3(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumn3"
     allocate(rumn4(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumn4"
     allocate(ren(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ren"
     allocate(rep(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rep"
     allocate(reR2c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating reR2c"
     allocate(reR3c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating reR3c"
     allocate(rut(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rut"
     allocate(rum(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rum"
     allocate(run(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating run"
     allocate(sun(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sun"
     allocate(rug(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rug"
     allocate(suR1(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating suR1"
     allocate(suR1n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating suR1n"
     allocate(suR1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating suR1p"
     allocate(suR2(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating suR2"
     allocate(cuR6(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating cuR6"
     allocate(cuR1(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating cuR1"
     allocate(iN1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iN1p"
     allocate(iNIn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iNIn"
     allocate(iN(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating iN"
     allocate(eN1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eN1p"
     allocate(eN4n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eN4n"
     allocate(huln(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating huln"
     allocate(hulp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating hulp"
     allocate(bacc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating bacc"
  end if

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Copy state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ppbacc = ppPelBacteria(bac,iiC)
  ppbacn = ppPelBacteria(bac,iiN)
  ppbacp = ppPelBacteria(bac,iiP)
  bacc = D3STATE(ppbacc,:)
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et = eTq_vector(ETW(:), p_q10(bac))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Oxygen environment: bacteria are both aerobic and anaerobic
  ! To provide a faster switching between the two metabolic pathways the
  ! oxygen regulating factor eO2 is cubic
  ! (eq. 19 in Vichi et al., 2004)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eO2 = MM_power_vector(max(p_small,O2o(:)),  p_chdo(bac),3)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! External nutrient limitation (used by some parametrizations)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eN4n = MM_vector(N4n(:), p_chn(bac))
  eN1p = MM_vector(N1p(:), p_chp(bac))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Mortality:
  !   1. first order mortality: p_sd 
  !   2. density dependent mortality due to virus infection: p_sd2
  !
  !   It is assumed that mortality is distributed in the same way over
  !   LOC (R1) and detritus (R6) as for phytoplankton and microzooplankton
  !   using the p_pe_R1x parameters defined in Param
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd  =  ( p_sd(bac)*et + p_sd2(bac)*bacc ) * bacc

  call flux_vector( iiPel, ppbacc,ppR6c, rd*(ONE-p_pe_R1c) )
  call flux_vector( iiPel, ppbacn,ppR6n, rd*qncPBA(bac,:)*(ONE-p_pe_R1n) )
  call flux_vector( iiPel, ppbacp,ppR6p, rd*qpcPBA(bac,:)*(ONE-p_pe_R1p) )

  call flux_vector( iiPel, ppbacc,ppR1c, rd*p_pe_R1c )
  call flux_vector( iiPel, ppbacn,ppR1n, rd*qncPBA(bac,:)*p_pe_R1n )
  call flux_vector( iiPel, ppbacp,ppR1p, rd*qpcPBA(bac,:)*p_pe_R1p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Substrate availability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version(bac) )

    case ( BACT3 )  ! Polimene et al. (2006)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria
      ! Note: oxygen control in eq 5 (Polimene et al. 2006) is not included
      !        as bacteria are both aerobic and anaerobic
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  =  p_sum(bac)*et*bacc

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! No correction of organic material quality
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      cuR1 = ONE
      cuR6 = ONE

    case ( BACT1,BACT2 )  ! Vichi et al. (2004,2007), Lazzari et al. (2012) 

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-==--=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-
      ! Nutrient limitation (intracellular, eq. 51 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      iNIn = min(ONE, max(ZERO, qncPBA(bac,:)/p_qncPBA(bac)))  !Nitrogen
      iN1p = min(ONE, max(ZERO, qpcPBA(bac,:)/p_qpcPBA(bac)))  !Phosphorus
      iN   = min(iN1p, iNIn)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Potential uptake by bacteria (eq. 50 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rum  = iN*et*p_sum(bac)*bacc

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! correction of substrate quality depending on nutrient content
      ! (eq. 52 Vichi et al. 2007)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      cuR1 = min(ONE, qpcOMT(iiR1,:)/p_qpcPBA(bac), qncOMT(iiR1,:)/ p_qncPBA(bac))
      cuR6 = min(ONE, qpcOMT(iiR6,:)/p_qpcPBA(bac), qncOMT(iiR6,:)/ p_qncPBA(bac))

  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate the realized substrate uptake rate depending on the
  ! type of detritus and quality (cuRx)
  ! See eq 27 in Vichi et al., 2004 for R2
  ! and eq 6 in Polimene et al., 2006 for R3 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ruR1c = (p_suhR1(bac)*cuR1 + p_sulR1(bac)*(ONE-cuR1))*R1c(:)
  ruR2c = p_suR2(bac)*R2c(:)
  ruR3c = p_suR3(bac)*R3c(:)
  ruR6c = p_suR6(bac)*cuR6*R6c(:)
  rut   = p_small + ruR1c + ruR2c + ruR6c + ruR3c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria (eq. 50 Vichi et al. 2007)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rug = min( rum, rut )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR1c = rug*ruR1c/rut
  ruR2c = rug*ruR2c/rut
  ruR3c = rug*ruR3c/rut
  ruR6c = rug*ruR6c/rut
  call flux_vector( iiPel, ppR1c, ppbacc, ruR1c )
  call flux_vector( iiPel, ppR2c, ppbacc, ruR2c )
  call flux_vector( iiPel, ppR3c, ppbacc, ruR3c )
  call flux_vector( iiPel, ppR6c, ppbacc, ruR6c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruR1n = qncOMT(iiR1,:)*ruR1c
  ruR6n = qncOMT(iiR6,:)*ruR6c
  call flux_vector( iiPel, ppR1n, ppbacn, ruR1n )
  call flux_vector( iiPel, ppR6n, ppbacn, ruR6n )

  ruR1p = qpcOMT(iiR1,:)*ruR1c
  ruR6p = qpcOMT(iiR6,:)*ruR6c
  call flux_vector( iiPel, ppR1p,ppbacp, ruR1p )
  call flux_vector( iiPel, ppR6p,ppbacp, ruR6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Aerobic and anaerobic respiration 
  ! Pelagic bacteria are a wide functional group comprising both aerobic and
  ! anaerobic bacteria. At (very) low Oxygen concentrations bacteria use
  ! N6r as electron acceptor in the respiration process. 
  ! However, the carbon cost is higher and an additional term is used.
  ! If nitrate is present, the rate of consumption of N6r is converted to N3n
  ! consumption (eq 19 Vichi et al., 2004 and PelChem.F90)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrc = (p_pu_ra(bac)+ p_pu_ra_o(bac)*(ONE-eO2) )*rug + p_srs(bac)* bacc* et
  call flux_vector( iiPel, ppbacc, ppO3c, rrc )
  call flux_vector( iiPel, ppO2o, ppO2o, -eO2*rrc/MW_C )
  flN6rPBA = (ONE- eO2)*rrc/ MW_C* p_qro
  call flux_vector( iiPel, ppN6r, ppN6r, flN6rPBA )
  ! Update the total rate of formation of reduction equivalent
  flPTN6r(:) = flPTN6r(:) + flN6rPBA

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Net Production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  run = rug - rrc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !             Fluxes from bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( p_version(bac))
    case ( BACT3 ) ! Polimene et al. (2006)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Carbon excretion as Semi-Labile (R2) and Semi-Refractory (R3) DOC 
      ! The R2 rate is assumed to occur with a timescale of 1 day
      ! (eq 8 Polimene et al., 2006)
      ! The renewal of capsular material is a constant rate, equivalent
      ! to about 1/4 of the respiration rate, ~5% of uptake 
      ! (Stoderegger and Herndl, 1998)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      reR2c = max((ONE-(qpcPBA(bac,:)/p_qpcPBA(bac))), &
              (ONE-(qncPBA(bac,:)/p_qncPBA(bac))))*p_rec(bac)
      reR2c = max(ZERO,reR2c)*bacc
      reR3c = rug*(ONE-p_pu_ra(bac))*(p_pu_ra(bac)*p_pu_ea_R3(bac))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
      ! This rate is assumed to occur with a timescale p_ruen=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ren = (qncPBA(bac,:) - p_qncPBA(bac))*bacc*p_ruen(bac)
      call flux_vector(iiPel, ppbacn, ppN4n,       ren*insw_vector( ren))
      call flux_vector(iiPel, ppN4n, ppbacn, -eN4n*ren*insw_vector(-ren))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
      ! This rate is assumed to occur with a timescale of p_ruep=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rep  =  (qpcPBA(bac,:) - p_qpcPBA(bac))*bacc*p_ruep(bac)
      call flux_vector(iiPel, ppbacp, ppN1p,       rep*insw_vector( rep))
      call flux_vector(iiPel, ppN1p, ppbacp, -eN1p*rep*insw_vector(-rep))

    case ( BACT1 ) ! Vichi et al. 2007

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! There is no Carbon excretion in Vichi et al. 2007
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      reR2c = ZERO
      reR3c = ZERO

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Nitrogen dynamics
      ! Direct uptake of ammonium if N excretion rate is negative (ren < 0)
      ! This rate is assumed to occur with a timescale p_ruen=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ren  =  (qncPBA(bac,:) - p_qncPBA(bac))*bacc*p_ruen(bac)
      call flux_vector(iiPel, ppbacn, ppN4n,       ren*insw_vector( ren))
      call flux_vector(iiPel, ppN4n, ppbacn, -eN4n*ren*insw_vector(-ren))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Dissolved Phosphorus dynamics
      ! Direct uptake of phosphate if P excretion rate is negative (rep < 0)
      ! This rate is assumed to occur with a timescale of p_ruep=1 day
      ! and controlled with a Michaelis-Menten function
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rep  =  (qpcPBA(bac,:) - p_qpcPBA(bac))*bacc*p_ruep(bac)
      call flux_vector(iiPel, ppbacp, ppN1p,       rep*insw_vector( rep))
      call flux_vector(iiPel, ppN1p, ppbacp, -eN1p*rep*insw_vector(-rep))

    case ( BACT2 ) ! Vichi et al. 2004

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium remineralization  (Eq. 28 Vichi et al. 2004, note that there 
      ! is a bug in the paper as there should be no division by B1c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      huln = (ruR6n + ruR1n) - p_qncPBA(bac)*run
      ren  = huln*insw_vector(huln)
      call flux_vector(iiPel, ppbacn, ppN4n, ren)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Nitrogen uptake  (Eq. 29 Vichi et al. 2004, there is a bug
      ! here as well, the min should be a max as the numbers are negative!)
      ! from Ammonium and Nitrate when N is not balanced (huln<0)
      ! (nitrate uptake with ammonium inhibition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rumn3 = p_qun(bac)*N3n(:)*bacc*(ONE-eN4n)
      rumn4 = p_qun(bac)*N4n(:)*bacc
      rumn  = rumn3 + rumn4
      ren   = max(-rumn,huln)*insw_vector(-huln)
      call flux_vector(iiPel, ppN4n, ppbacn, -ren*rumn4/rumn)
      call flux_vector(iiPel, ppN3n, ppbacn, -ren*rumn3/rumn)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Phosphate remineralization  (Eq. 28 Vichi et al. 2004, note that there 
      ! is an error in the paper as there should be no division by B1c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      hulp = (ruR6p + ruR1p) - p_qpcPBA(bac)*run
      rep  = hulp*insw_vector(hulp)
      call flux_vector(iiPel, ppbacp, ppN1p, rep)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Inorganic Phosphorus uptake  (Eq. 29 Vichi et al. 2004, there is a bug
      ! here as well, the min should be a max as the numbers are negative!)
      ! from dissolved phosphate whene P is not balanced (hulp<0)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rump = p_qup(bac)*N1p(:)*bacc
      rep  = max(-rump,hulp)*insw_vector(-hulp)
      call flux_vector(iiPel, ppN1p, ppbacp, -rep)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Excess carbon (also considering dissolved nutrient uptake ren and rep) 
      ! is released as R3c, no other excretion (reR2c=0)
      ! (eq. 30 Vichi et al. 2004, unfortunately there is another error in 
      ! the paper, the flux of dissolved nutrient is not written)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      r     = min(run, (ruR6n+ruR1n-ren)/p_qlnc(bac))
      reR3c = run - min(r, (ruR6p+ruR1p-rep)/p_qlpc(bac))
      reR3c = max(ZERO, reR3c)
      reR2c = ZERO

  end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion fluxes (only losses to R2 and R3)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector( iiPel, ppbacc,ppR2c, reR2c )
  call flux_vector( iiPel, ppbacc,ppR3c, reR3c )

  end subroutine PelBacDynamics

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
