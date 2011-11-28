#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MesoZoo
!
! DESCRIPTION
!   This submodel describes the carbon dynamics and associated
!    nutrient dynamics in carnivorous mesozooplankton (represented
!    by the state variable Z3) and in omnivorous zooplankton (in
!    the model known as Z4).
!    
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine MesoZooDynamics(zoo,  ppzooc, ppzoon, ppzoop)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: O2o, N1p, N4n, R6c, &
  ! R6p, R6n
  ! For the following Pelagic-group-states fluxes are &
  ! defined: PhytoPlankton, MicroZooPlankton, MesoZooPlankton
  ! The following Pelagic 1-d global boxvars  are used: ETW
  ! The following Pelagic 2-d global boxvars are used: qnPc, qpPc, qlPc, qsPc, &
  ! qn_mz, qp_mz, qnZc, qpZc
  ! The following groupmember vars are used: iiPhytoPlankton, &
  ! iiMicroZooPlankton, iiMesoZooPlankton
  ! The following constituent constants  are used: iiC, iiN, iiP, iiL
  ! The following 0-d global parameters are used: p_small
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN, ONE, ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: D3STATE, O2o, N1p, N4n, R6c, R6p, R2c, &
    R6n, PhytoPlankton, MicroZooPlankton, MesoZooPlankton
  use mem, ONLY: Depth, ppO2o, ppO3c, ppN1p, ppN4n, ppR6c, ppR6n, ppR6p, ppR6s, &
    ppPhytoPlankton, ppMicroZooPlankton, ppMesoZooPlankton, ETW, &
    qnPc, qpPc, qlPc, qsPc, qn_mz, qp_mz, qnZc, qpZc, iiPhytoPlankton, &
    iiMicroZooPlankton, iiMesoZooPlankton, iiC, iiN, iiP, iiL, iiS, NO_BOXES, &
    iiBen, iiPel, flux_vector,fixed_quota_flux_vector
#endif
#ifdef BFM_GOTM
  use mem, ONLY: jnetMeZc
#endif
  use mem_Param,  ONLY: p_small,check_fixed_quota
  use constants,ONLY: MIN_VAL_EXPFUN, MW_C
  use mem_MesoZoo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector, MM_power_vector

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
!   N. Broekhuizen and A.D. Bryant, ERSEM group
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer,dimension(NO_BOXES)  :: nut_lim
  integer, save :: first =0
  real(RLEN),allocatable,save,dimension(:) :: put_u,temp_p,temp_n,rumc,rugc,eo,  &
                                       et,rrs_c,rrs_n,rrs_p,rra_c,rra_n,rra_p,rut_c, &
                                       rut_n,rut_p,rd_c,rd_n,rd_p,sdo,rdo_c,  &
                                       rdo_n,rdo_p,ret_c,ret_n,ret_p,ru_c, &
                                       ru_n,ru_p,pu_e_n,pu_e_p,prI_R6,pe_R6c

  real(RLEN),allocatable,save,dimension(:) :: pe_N1p,pe_N4n,ruPIc,ruMIZc,ruMEZc,rq6c, &
                                       rq6n,rq6p,rrc,ren,rep,tfluxc,tfluxn,    &
                                       tfluxp,zooc
  real(RLEN),allocatable,save,dimension(:,:) :: rumPIc,rumMIZc,rumMEZc
  real(RLEN),allocatable,save,dimension(:) :: net,r
  integer :: AllocStatus, DeallocStatus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocate local memory
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (first==0) then
     first=1
     allocate(eo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eo"
     allocate(rumPIc(NO_BOXES,4),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumPIc"
       allocate(rumMIZc(NO_BOXES,2),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumMIZc"
       allocate(rumMEZc(NO_BOXES,2),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumMEZc"
       allocate(zooc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating zooc"
     allocate(tfluxp(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxp"
     allocate(tfluxn(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxn"
     allocate(tfluxc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating tfluxc"
     allocate(rep(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rep"
     allocate(ren(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ren"
     allocate(rrc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrc"
     allocate(rq6p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rq6p"
     allocate(rq6n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rq6n"
     allocate(rq6c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rq6c"
     allocate(ruMEZc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruMEZc"
     allocate(ruMIZc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruMIZc"
     allocate(ruPIc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ruPIc"
     allocate(pe_N4n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pe_N4n"
     allocate(pe_N1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pe_N1p"
     allocate(pe_R6c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pe_R6c"
     allocate(prI_R6(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating prI_R6"
     allocate(pu_e_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pu_e_p"
     allocate(pu_e_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating pu_e_n"
     allocate(ru_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ru_p"
     allocate(ru_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ru_n"
     allocate(ru_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ru_c"
     allocate(ret_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ret_p"
     allocate(ret_n(NO_BOXES),stat=AllocStatus)
      if (AllocStatus  /= 0) stop "error allocating ret_n"
     allocate(ret_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating ret_c"
     allocate(rdo_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rdo_p"
     allocate(rdo_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rdo_n"
     allocate(rdo_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rdo_c"
     allocate(sdo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating sdo"
     allocate(rd_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rd_p"
     allocate(rd_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rd_n"
     allocate(rd_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rd_c"
     allocate(rut_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rut_p"
     allocate(rut_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rut_n"
     allocate(rut_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rut_c"
     allocate(rra_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rra_p"
     allocate(rra_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rra_n"
     allocate(rra_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rra_c"
     allocate(rrs_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrs_p"
     allocate(rrs_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrs_n"
     allocate(rrs_c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rrs_c"
     allocate(et(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating et"
     allocate(put_u(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating put_u"
     allocate(temp_p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating temp_p"
     allocate(temp_n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating temp_n"
     allocate(rumc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rumc"
     allocate(rugc(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rugc"
     allocate(net(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating net"
     allocate(r(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating r"
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  zooc = D3STATE(ppzooc,:)

  tfluxc=ZERO
  tfluxn=ZERO
  tfluxp=ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !Physiological temperature and oxygen response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eo  =   MM_power_vector(  max(p_small,O2o(:)),  p_clO2o(zoo),3)
  et  =   eTq_vector(  ETW(:),  p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! rua: food uptake generalized. Addition of new FG become much more simpler!


  rumc  =   ZERO
  do i = 1 , ( iiPhytoPlankton)
#ifdef BFM_NS
    r = ONE-min(ONE,(R2c(:)*R2c(:))/40000.0_RLEN)
#else
    r = ONE
#endif

    rumPIc(:, i)  =   r* p_puPI(zoo,i)* PhytoPlankton(i,iiC)
    rumc  =   rumc+ rumPIc(:, i)
  end do


  do i = 1 , ( iiMicroZooPlankton)

    rumMIZc(:, i)  =   p_puMIZ(zoo,i)* MicroZooPlankton(i,iiC)
    rumc  =   rumc+ rumMIZc(:, i)
  end do


  do i = 1 , ( iiMesoZooPlankton)

    rumMEZc(:, i)  =   p_puMEZ(zoo,i)* MesoZooPlankton(i,iiC)
    rumc  =   rumc+ rumMEZc(:, i)
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  =   eo *et* p_sum(zoo)* MM_vector(  p_vum(zoo)* rumc,  p_sum(zoo))* zooc
  put_u  =   rugc/ ( 1.0D-80 + rumc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptakes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rut_c  =   ZERO
  rut_n  =   ZERO
  rut_p  =   ZERO
  do i = 1 , ( iiPhytoPlankton)

    ruPIc  =   put_u* rumPIc(:, i)
    call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                          ppPhytoPlankton(i,iiC),ppzooc, ruPIc ,tfluxc)
    call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                ppPhytoPlankton(i,iiN),ppzoon, ruPIc* qnPc(i,:),tfluxn )
    call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
              ppPhytoPlankton(i,iiP),ppzoop, ruPIc* qpPc(i,:),tfluxp )
    rut_c  =   rut_c+ ruPIc
    rut_n  =   rut_n+ ruPIc* qnPc(i,:)
    rut_p  =   rut_p+ ruPIc* qpPc(i,:)
    ! Chl is transferred to the sink
    call flux_vector( iiPel, ppPhytoPlankton(i,iiL), &
               ppPhytoPlankton(i,iiL),-( ruPIc* qlPc(i,:)) )
    ! PIs is directly transferred to R6s
    if ( ppPhytoPlankton(i,iiS) .gt. 0 ) & 
    call flux_vector( iiPel, ppPhytoPlankton(i,iiS), &
               ppR6s,+(ruPIc* qsPc(i,:)) )
    
  end do


  do i = 1 , ( iiMicroZooPlankton)

    ruMIZc  =   put_u* rumMIZc(:, i)
    call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
             ppMicroZooPlankton(i,iiC),ppzooc, ruMIZc,tfluxc )
    call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
           ppMicroZooPlankton(i,iiN),ppzoon, ruMIZc*qn_mz(i,:) ,tfluxn)
    call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
           ppMicroZooPlankton(i,iiP),ppzoop, ruMIZc* qp_mz(i,:) ,tfluxp)
    rut_c  =   rut_c+ ruMIZc
    rut_n  =   rut_n+ ruMIZc* qn_mz(i,:)
    rut_p  =   rut_p+ ruMIZc* qp_mz(i,:)
  end do


  do i = 1 , ( iiMesoZooPlankton)

    ruMEZc  =   put_u* rumMEZc(:, i)
    ! intra-group predation is not computed
    if ( i/= zoo) then
      call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                       ppMesoZooPlankton(i,iiC),ppzooc, ruMEZc, tfluxc )
      call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
            ppMesoZooPlankton(i,iiN),ppzoon, ruMEZc* qnZc(i,:) , tfluxn)
      call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
            ppMesoZooPlankton(i,iiP),ppzoop, ruMEZc* qpZc(i,:) , tfluxp)
    end if

    rut_c  =   rut_c+ ruMEZc
    rut_n  =   rut_n+ ruMEZc* qnZc(i,:)
    rut_p  =   rut_p+ ruMEZc* qpZc(i,:)
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Proportion of ingested food respired by zoo = prIR6/Z4R6
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  prI_R6  =   ONE- p_puI_u(zoo)- p_peI_R6(zoo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assimilated material
  ! Respectively Carbon, Nitrogen and Phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ru_c  =   p_puI_u(zoo)* rut_c
  ru_n  =  ( p_puI_u(zoo)+ prI_R6)* rut_n
  ru_p  =  ( p_puI_u(zoo)+ prI_R6)* rut_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! P:C and N:C ratios in assimilate
  ! Nitrogen & Phosphorus:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pu_e_n  =   ru_n/( p_small+ ru_c)
  pu_e_p  =   ru_p/( p_small+ ru_c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Determine whether C, P or N is the Limiting Nutrient. Variable nut_lim
  ! holds the kind of nutrient limitation (1, 2, 3).
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  nut_lim  =   1

  temp_p  =   pu_e_p/ qpZc(zoo,:)
  temp_n  =   pu_e_n/ qnZc(zoo,:)



    WHERE ( (( temp_p< temp_n) .OR.( abs(temp_p- temp_n)< p_small)) )
      where ( pu_e_p< qpZc(zoo,:))
        nut_lim  =   2
      end where

    ELSEWHERE
      where ( pu_e_n< qnZc(zoo,:))
        nut_lim  =   3
      end where

  END WHERE


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration and basal metabolism
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rra_c  =   prI_R6* rut_c
  rra_n  =   ZERO
  rra_p  =   ZERO

  rrs_c  =   p_srs(zoo)* et*eo* zooc
  rrs_n  =   p_srs(zoo)* et*eo* zooc * qnZc(zoo,:)
  rrs_p  =   p_srs(zoo)* et*eo* zooc * qpZc(zoo,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Defecation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ret_c  =   p_peI_R6(zoo)* rut_c
  ret_n  =   p_peI_R6(zoo)* rut_n
  ret_p  =   p_peI_R6(zoo)* rut_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Natural mortality + low oxygen mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rd_c  =   (p_sd(zoo)+ p_srs(zoo)*(ONE-eo))*et* zooc
  rd_n  =   (p_sd(zoo)+ p_srs(zoo)*(ONE-eo))* zooc* et * qnZc(zoo,:)
  rd_p  =   (p_sd(zoo)+ p_srs(zoo)*(ONE-eo))* zooc * et* qpZc(zoo,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Density dependent mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sdo  =   p_sdo(zoo)* (zooc)**(p_sds(zoo))
  rdo_c  =   sdo* zooc
  rdo_n  =   sdo* zooc* qnZc(zoo,:)
  rdo_p  =   sdo* zooc* qpZc(zoo,:)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Eliminate excess of non-limiting nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


    WHERE (( nut_lim)==1)
      pe_R6c  =  ZERO
      pe_N1p = max( ZERO, ( ONE- p_peI_R6(zoo))* rut_p- p_qpc(zoo)* &
        ru_c)
      pe_N1p  =   pe_N1p/( p_small+ rut_p)
      pe_N4n = max( ZERO, ( ONE- p_peI_R6(zoo))* rut_n- p_qnc(zoo)* &
        ru_c)
      pe_N4n  =   pe_N4n/( p_small+ rut_n)

    ELSEWHERE (( nut_lim)==2)
      pe_N1p  =   ZERO
      pe_R6c  =  max(ZERO,( p_qpc(zoo)* ru_c)-( ONE- p_peI_R6(zoo))* rut_p)
      pe_R6c  =   pe_R6c/( p_small+ p_qpc(zoo)* rut_c)
      pe_N4n = max( ZERO, ( ONE- p_peI_R6(zoo))* rut_n- p_qnc(zoo)*( &
        ru_c- pe_R6c* rut_c))
      pe_N4n  =   pe_N4n/( p_small+ rut_n)

    ELSEWHERE (( nut_lim)==3)
      pe_N4n  =   ZERO
      pe_R6c  = max(ZERO, ( p_qnc(zoo)* ru_c)-( ONE- p_peI_R6(zoo))* rut_n)
      pe_R6c  =   pe_R6c/( p_small+ p_qnc(zoo)* rut_c)
      pe_N1p = max( ZERO, ( ONE- p_peI_R6(zoo))* rut_p- p_qpc(zoo)*( &
        ru_c- pe_R6c* rut_c))
      pe_N1p  =   pe_N1p/( p_small+ rut_p)

  END WHERE
  ! End of select(nut_lim)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes for eliminated excess nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rq6c  =   rd_c+ ret_c+ rdo_c+ pe_R6c* rut_c
  rq6p  =   rd_p+ ret_p+ rdo_p
  rq6n  =   rd_n+ ret_n+ rdo_n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration: activity + basal metabolism
  ! Excretion: activity + basal metabolism + excess non-limiting nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrc  =   rra_c+ rrs_c
  ren  =   rra_n+ rrs_n+ pe_N4n* rut_n
  rep  =   rra_p+ rrs_p+ pe_N1p* rut_p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! flow statements
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrc/ MW_C) )
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                             ppzooc,ppO3c, rrc,tfluxc )
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
                             ppzoop,ppN1p, rep ,tfluxp)
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                             ppzoon,ppN4n, ren ,tfluxn)

  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                             ppzooc,ppR6c, rq6c ,tfluxc )
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
                             ppzoop,ppR6p, rq6p ,tfluxp)
  call fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                             ppzoon,ppR6n, rq6n ,tfluxn)

  ren=tfluxC*p_qnc(zoo)
  call fixed_quota_flux_vector( check_fixed_quota,-iiN,0,0,0,ren,tfluxN)
  rep=tfluxC*p_qpc(zoo)
  call fixed_quota_flux_vector( check_fixed_quota,-iiP,0,0,0,rep,tfluxP)

#ifdef BFM_GOTM
  net=rut_c-rra_c-ret_c-pe_R6c* rut_c
  jnetMeZc(1)=jnetMeZc(1)+sum(Depth(:)*net)
#endif
  end subroutine MesoZooDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
