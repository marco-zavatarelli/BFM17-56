#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: FilterFeeder
!
! DESCRIPTION
!   This process describes the carbon dynamics and associated
!   nutrient dynamics in benthic organism Y3 (suspension feeders)
!   Y3 is handled separately because it also feeds from the water column.
!
!
! !INTERFACE
  subroutine FilterFeederDynamics
!
#ifdef INCLUDE_BEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
  use constants, ONLY:MW_C
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: Y3c, Y3n, Y3p, Q6c, Q6n, Q6p, Q6s, G2o, K4n, K1p, D6m, D7m, D8m, D9m, &
    D1m
  use mem, ONLY: ppY3c, ppY3n, ppY3p, ppQ6c, ppQ6n, ppQ6p, ppQ6s, ppG2o, ppK4n,O2o_Ben, &
    ppK1p, ppD6m, ppD7m, ppD8m, ppD9m, ppD1m, rrBTo, reBTn, reBTp, jbotR6c, jbotR6n, &
    jbotR6p, jbotR6s, jPIY3c, jZIY3c, jRIY3c, jRIY3n, jRIY3p, jRIY3s, ETW_Ben, &
    iiPhytoPlankton, PI_Benc, PI_Benn, PI_Benp, PI_Bens, sediPPY_Ben, sediR6_Ben, & 
    ZI_Fc, RI_Fc, ZI_Fn, ZI_Fp, RI_Fn, RI_Fp, RI_Fs, ppG3c, &
    NO_BOXES_XY, Depth_ben, iiBen, iiPel, flux_vector, jbotO2o,jbotN1p,jbotN4n
#ifdef INCLUDE_BENCO2
  use mem, ONLY: jbotO3c
#endif
 use mem,  ONLY: Source_D2_vector_ben
#endif
  use mem_Param,  ONLY: p_d_tot,p_pe_R1c, p_pe_R1n, p_pe_R1p,p_small
  use mem_FilterFeeder

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, eramp_vector, &
  ! MM_vector, PartQ_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, eramp_vector, MM_vector, MM_power_vector, insw_vector, PartQ_vector

!
! !AUTHORS
!   W. Ebenhoeh and C. Kohlmeier 
!
!
! !REVISION_HISTORY
!   !
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  real(RLEN) :: clu
  real(RLEN),dimension(NO_BOXES_XY)  :: corr
  real(RLEN),dimension(NO_BOXES_XY)  :: fdepth
  real(RLEN),dimension(NO_BOXES_XY)  :: clm
  real(RLEN),dimension(NO_BOXES_XY)  :: cmm
  real(RLEN),dimension(NO_BOXES_XY)  :: et
  real(RLEN),dimension(NO_BOXES_XY)  :: eO
  real(RLEN),dimension(NO_BOXES_XY)  :: eNC
  real(RLEN),dimension(NO_BOXES_XY)  :: ePC
  real(RLEN),dimension(NO_BOXES_XY)  :: foodpm2
  real(RLEN),dimension(NO_BOXES_XY)  :: food
  real(RLEN),dimension(iiPhytoPlankton,NO_BOXES_XY)  :: food_PIc
  real(RLEN),dimension(NO_BOXES_XY)  :: food_PT
  real(RLEN),dimension(NO_BOXES_XY)  :: food_ZI
  real(RLEN),dimension(NO_BOXES_XY)  :: food_RI
  real(RLEN),dimension(NO_BOXES_XY)  :: food_Q6
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_c
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_n
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_p
  real(RLEN),dimension(NO_BOXES_XY)  :: eF
  real(RLEN),dimension(NO_BOXES_XY)  :: sgu
  real(RLEN),dimension(NO_BOXES_XY)  :: rgu
  real(RLEN),dimension(NO_BOXES_XY)  :: snuPI
  real(RLEN),dimension(NO_BOXES_XY)  :: snuZI
  real(RLEN),dimension(NO_BOXES_XY)  :: snuQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uPI
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uZI
  real(RLEN),dimension(NO_BOXES_XY)  :: se_uQ6
  real(RLEN),dimension(NO_BOXES_XY)  :: choice
  real(RLEN),dimension(NO_BOXES_XY)  :: rtY3c
  real(RLEN),dimension(NO_BOXES_XY)  :: rtY3n
  real(RLEN),dimension(NO_BOXES_XY)  :: rtY3p
  real(RLEN),dimension(NO_BOXES_XY)  :: rrc
  real(RLEN),dimension(NO_BOXES_XY)  :: sm
  real(RLEN),dimension(NO_BOXES_XY)  :: ren
  real(RLEN),dimension(NO_BOXES_XY)  :: rep
  real(RLEN),dimension(NO_BOXES_XY)  :: retR6c
  real(RLEN),dimension(NO_BOXES_XY)  :: retR6n
  real(RLEN),dimension(NO_BOXES_XY)  :: retR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: retQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: retQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: retQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIc
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIn
  real(RLEN),dimension(NO_BOXES_XY)  :: rePIp
  real(RLEN),dimension(NO_BOXES_XY)  :: reZIc
  real(RLEN),dimension(NO_BOXES_XY)  :: reZIn
  real(RLEN),dimension(NO_BOXES_XY)  :: reZIp
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6c
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6n
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: reR6s
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: reQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIc
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIn
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIp
  real(RLEN),dimension(NO_BOXES_XY)  :: ruPIs
  real(RLEN),dimension(NO_BOXES_XY)  :: ruZIc
  real(RLEN),dimension(NO_BOXES_XY)  :: ruZIn
  real(RLEN),dimension(NO_BOXES_XY)  :: ruZIp
  real(RLEN),dimension(NO_BOXES_XY)  :: RTc
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6p
  real(RLEN),dimension(NO_BOXES_XY)  :: ruR6s
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: su
  real(RLEN),dimension(NO_BOXES_XY)  :: r
  real(RLEN),dimension(NO_BOXES_XY)  :: puf
  real(RLEN),dimension(NO_BOXES_XY)  :: fsat ! filtering saturation : at high feed levels less filtering 
                                             ! is necessairy 
  real(RLEN),dimension(NO_BOXES_XY)  :: netto
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW_Ben(:),  p_q10)

  eo  =   MM_power_vector(  max(p_small,O2o_Ben(:)),  p_clO2o,3)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food Cfluxes!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  clu=p_clu+p_small;
  if ( sw_uptake == 1 ) clu=p_clu/p_dwat;
  food  =   p_small

  ! For phytoplankton:

  food_PT=ZERO
  do i=1,iiPhytoPlankton
     r =  PI_Benc(i,:) * MM_vector(  PI_Benc(i,:),  clu)
     call CorrectConcNearBed(Depth_Ben(:), sediPPY_Ben(i,:), p_height, & 
                                    p_max, p_vum*et*Y3c(:), corr)
     food_PIc(i,:)=r*corr*p_PI
     food_PT(:)  =   food_PT(:)+ food_PIc(i,:)
  enddo
  food  =   food  + food_PT(:)

  ! For microzooplankton:

  food_ZI  =   p_ZI * ZI_Fc(:) * MM_vector(  ZI_Fc(:),  clu)
  food  =   food+ food_ZI


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus, (if eaten) first calculate available amount
  ! and add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  r=   RI_Fc(:)* MM_vector(  RI_Fc(:),  clu)
  call CorrectConcNearBed(Depth_Ben(:), sediR6_Ben(:), p_height, & 
                                    p_max, p_vum*et*Y3c(:), corr)
  RTc=r*corr
  food_RI=RTc*p_R6
  food  =   food+ food_RI

  !

  select case (sw_uptake)
   case(1)
    ! This uptake procedure was developed for the one/two layer orginal ERSEM
    ! model where the layer above the sediment could have depths upto a few
    ! hunder meters. p_dwat is in this case the layer depth seen by the
    ! filterfeeders. o_dwat is used as an imporatent calibration parameter.

    ! In the orginal model the sedimentation of detritus (R6) was equal to the
    ! the sinking rate. By doing this we implicetly assumed that this rate was
    ! a gros seimentation rate. Therefor in the orignal setup filter took
    ! also food from the benthic system.

    fdepth=p_dwat
    foodpm2 =food*fdepth
    clm  =   p_clm
    cmm  =   p_cm
    availQ6_c  =   Q6c(:)* PartQ_vector(  D6m(:),  clm,  cmm,  p_d_tot)
    availQ6_n  =   Q6n(:)* PartQ_vector(  D7m(:),  clm,  cmm,  p_d_tot)
    availQ6_p  =   Q6p(:)* PartQ_vector(  D8m(:),  clm,  cmm,  p_d_tot)

    food_Q6  =   p_puQ6* availQ6_c* MM_vector(  availQ6_c,  clu)
    foodpm2  =   foodpm2+ food_Q6

    cmm  =  ( p_clm+ p_cm)* 0.5D+00

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Correct for too much food:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     eF  =   MM_vector(  foodpm2,  p_chu)

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Correction of growth rate for environmental factors:
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     ! The minimal uptake rate is equal to rest respiration.
     ! With filtering the filterfeeder provide himself also with oxygen.
     rgu  =max( p_su* eO* eF,p_srr)* Y3c(:)* et
    
     fsat=ONE;
     rrc = max(eo * p_sra, p_srr)* Y3c(:)* et

   case(2)
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Alternative food uptake as in zooplankton:
     !  using the modifed Holling response equation which tkae in account
     !  the maximum growth rate and the volume filtered.
     !
     !  It is assumed that the detritus sedimentation is defined as a netto ptocess
     !  ( p_bursel << P_sediR6). Therefor it assumed that filterfeeders do noet eat Q6.  
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     cmm = ZERO;

     fdepth=Depth_Ben(:)
     su  =  et* eO*  p_su* MM_vector(  p_vum* food,  p_su)
     fsat=min(ONE,su/(et*eo*p_vum*food))
     rgu= su *Y3c(:)
     rrc = max(eo * p_sra*fsat, p_srr)* Y3c(:)* et
     foodpm2 =food*fdepth
   case(3)
     fdepth=Depth_Ben(:)
     su  =  p_su* MM_vector(  p_vum* food,  p_su)
     netto= (ONE-(p_pueQ6*food_RI+p_puePI*food_PT +p_pueZI*food_ZI)/food ) * (ONE-p_pur);
     su=su * insw_vector(netto * su-p_sra);
     fsat=min(ONE,su/(p_vum*food));
     rgu= et* eO*  su *Y3c(:)
     rrc = max(eo * p_sra*fsat, p_srr)* Y3c(:)* et
     foodpm2 =food*fdepth
   case(4)
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Alternative food uptake as in zooplankton:
     !  using the modifed Holling response equation which tkae in account
     !  the maximum growth rate and the volume filtered.
     !  Further is assumed that the filterfeeder (nearly) stop filtering as soon as  
     !  the costs for filtering  are lower than the  profit
     !  For this we solve the next equation in which r is the unknown:
     !    (left side == profit , right side=costs)
     !    r* p_su* MM_vector(  p_vum* food,  r* p_su)* Y3c(:)*netto = p_sra *r 
     !  If r > ONE : there is enough food to grow
     !  if r < ONE : there is balance between costs and profit if  r*p_puf*sgu is larger than
     !  the rest respiration.
     !
     !  It is assumed that the detritus sedimentation is defined as a netto ptocess
     !  ( p_bursel << P_sediR6). Therefor it assumed that filterfeeders do noet eat Q6.  
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     cmm = ZERO;

     fdepth=Depth_Ben(:)
     puf = p_sra/p_su;

     netto= (ONE-(p_pueQ6*food_RI+p_puePI*food_PT +p_pueZI*food_ZI)/food ) * (ONE-p_pur);
     r= netto /puf - p_su/(p_vum *food);

     r=min(ONE,max(1.0e-6_RLEN,r));

     ! Calculate relative uptake
     su= ONE/( ONE/(p_small+ r* p_vum *food * et *eO ) + ONE/(p_small+ p_su *et *eO  ))  ;
     ! The minimal uptake rate is equal to rest respiration.
     ! With filtering the filterfeeder provide himself also with oxygen.
     rgu= su *Y3c(:)

     ! filtering saturation ( high at low , low at hight food)
     fsat=su/(p_small+ et*eo*p_vum*food);
     ! Calculate cost of enregy based on realized rate of uptake.
     rrc = max(eo * p_su*puf*fsat, p_srr)* Y3c(:)* et
     
     foodpm2 =food*fdepth
  end select

  ! Relative growth rate corrected for actual amount of food:

  sgu  =   rgu/ foodpm2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Net uptake:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  snuPI  =   sgu*( ONE- p_puePI)
  snuZI  =   sgu*( ONE- p_pueZI)
  snuQ6  =   sgu*( ONE- p_pueQ6)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Execreted part:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  se_uPI  =   sgu- snuPI
  se_uZI  =   sgu- snuZI
  se_uQ6  =   sgu- snuQ6

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of uptake rate:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  eNC=(ONE-p_pe_R1n)/(ONE-p_pe_R1c)
  ePC=(ONE-p_pe_R1p)/(ONE-p_pe_R1c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic Phytoplankton:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruPIc  =   ZERO
  ruPIn  =   ZERO
  ruPIp  =   ZERO
  ruPIs  =   ZERO

  rePIc  =   ZERO
  rePIn  =   ZERO
  rePIp  =   ZERO

  do i=1,iiPhytoPlankton
    choice=food_PIc(i,:)* fdepth/(p_small + PI_Benc(i,:))
    jPIY3c(i,:) =       PI_Benc(i,:)* sgu* choice
    ruPIc  = ruPIc  +   PI_Benc(i,:)* sgu* choice
    ruPIn  = ruPIn  +   PI_Benn(i,:)* sgu* choice
    ruPIp  = ruPIp  +   PI_Benp(i,:)* sgu* choice
    ruPIs  = ruPIs  +   PI_Bens(i,:)* sgu* choice

    rePIc  = rePIc  +   PI_Benc(i,:)* se_uPI* choice
    rePIn  = rePIn  +   PI_Benn(i,:)* se_uPI* choice*eNC 
    rePIp  = rePIp  +   PI_Benp(i,:)* se_uPI* choice*ePC 
  enddo

  call flux_vector( iiBen, ppY3c,ppY3c, ruPIc )
  call flux_vector( iiBen, ppY3n,ppY3n, ruPIn )
  call flux_vector( iiBen, ppY3p,ppY3p, ruPIp )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic MicroZooplankton:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  choice  =   p_ZI* MM_vector(  ZI_Fc(:),  clu)* fdepth

  ruZIc  =   ZI_Fc(:)* sgu* choice
  ruZIn  =   ZI_Fn(:)* sgu* choice
  ruZIp  =   ZI_Fp(:)* sgu* choice

  call flux_vector( iiBen, ppY3c,ppY3c, ruZIc )
  call flux_vector( iiBen, ppY3n,ppY3n, ruZIn )
  call flux_vector( iiBen, ppY3p,ppY3p, ruZIp )

  ! flux definitions from Z -> Y3 are found in BentoPelCoup
  jZIY3c(:)  =   ruZIc

  reZIc  =   ZI_Fc(:)* se_uZI* choice
  reZIn  =   ZI_Fn(:)* se_uZI* choice* eNC
  reZIp  =   ZI_Fp(:)* se_uZI* choice* ePC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic Detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  choice  =   food_RI * fdepth/(p_small+RI_Fc(:))

  ruR6c  =   RI_Fc(:)* sgu* choice
  ruR6n  =   RI_Fn(:)* sgu* choice
  ruR6p  =   RI_Fp(:)* sgu* choice
  ruR6s  =   RI_Fs(:)* sgu* choice

  call flux_vector( iiBen, ppY3c,ppY3c, ruR6c )
  call flux_vector( iiBen, ppY3n,ppY3n, ruR6n )
  call flux_vector( iiBen, ppY3p,ppY3p, ruR6p )

  reR6c  =   RI_Fc(:)* se_uQ6* choice
  reR6n  =   RI_Fn(:)* se_uQ6* choice *eNC
  reR6p  =   RI_Fp(:)* se_uQ6* choice *ePC


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Benthic Detritus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( sw_uptake == 1) then
    choice  =   p_puQ6* MM_vector(  availQ6_c,  clu)

    ruQ6c  =   sgu* choice* availQ6_c
    ruQ6n  =   sgu* choice* availQ6_n
    ruQ6p  =   sgu* choice* availQ6_p

    call flux_vector( iiBen, ppQ6c,ppY3c, ruQ6c )
    call flux_vector( iiBen, ppQ6n,ppY3n, ruQ6n )
    call flux_vector( iiBen, ppQ6p,ppY3p, ruQ6p )

    reQ6c  =   se_uQ6* availQ6_c* choice
    reQ6n  =   se_uQ6* availQ6_n* choice *eNC
    reQ6p  =   se_uQ6* availQ6_p* choice *ePC

  else
    reQ6c=ZERO
    reQ6n=ZERO
    reQ6p=ZERO
    ruQ6c=ZERO
    ruQ6n=ZERO
    ruQ6p=ZERO
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Book keeping
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rtY3c  =   ruPIc+ ruZIc+ ruR6c+ ruQ6c
  rtY3n  =   ruPIn+ ruZIn+ ruR6n+ ruQ6n
  rtY3p  =   ruPIp+ ruZIp+ ruR6p+ ruQ6p

  retR6c  =   rePIc+ reZIc+ reR6c
  retR6n  =   rePIn+ reZIn+ reR6n
  retR6p  =   rePIp+ reZIp+ reR6p

  retQ6c  =   reQ6c
  retQ6n  =   reQ6n
  retQ6p  =   reQ6p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrc  =   rrc+ p_pur*( foodpm2* sgu- retR6c- retQ6c)
  
  !jnetY3c(:)=rtY3c(:)- p_pur*( foodpm2* sgu- retR6c- retQ6c)-retQ6c-retR6c


  call flux_vector( iiBen, ppY3c,ppG3c, rrc*(ONE-p_pePel) )
  call flux_vector(iiBen, ppG2o,ppG2o,-( rrc/ MW_C)* ( ONE-p_pePel))
  call flux_vector(iiBen, ppY3c,ppY3c, -rrc*p_pePel )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- =-=-=-=-=-=-=-=-=-=-=-=-=-

  sm  =   p_sd* et  +p_sd2 * Y3c(:)

  reQ6c  =   Y3c(:)* sm
  reQ6n  =   Y3n(:)* sm
  reQ6p  =   Y3p(:)* sm

  retQ6c  =   retQ6c+ reQ6c
  retQ6n  =   retQ6n+ reQ6n
  retQ6p  =   retQ6p+ reQ6p

  ! in case of a negative value of one of the following values there is a &
  ! situation
  ! of startvation and very low biomass values. Check on quota in the food is &
  ! out of order

  rtY3c  =   max(  ZERO,  rtY3c -retR6c-retQ6c-rrc)
  rtY3n  =   max(  ZERO,  rtY3n -retR6n-retQ6n)
  rtY3p  =   max(  ZERO,  rtY3p -retR6p-retQ6p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of nutrient release and correction of C:N:P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ren  =   rtY3n- rtY3c* p_qn 
  rep  =   rtY3p- rtY3c* p_qp 

  r=ZERO

  where ( ren< ZERO)

    reQ6c  =  - ren/ p_qn
    r      =   r + reQ6c
    rtY3c  =   rtY3c- reQ6c

    ren  =   rtY3n- rtY3c* p_qn
    rep  =   rtY3p- rtY3c* p_qp

  end where


  where ( rep< ZERO)

    reQ6c  =  - rep/ p_qp
    r      =   r   + reQ6c
    rtY3c  =   rtY3c- reQ6c

    ren  =   rtY3n- rtY3c* p_qn
    rep  =   rtY3p- rtY3c* p_qp

  end where
 
  retQ6c= retQ6c +r * retQ6c/(p_small + retQ6c + retR6c);
  retR6c= retR6c +r * retR6c/(p_small + retQ6c + retR6c);


  ren = max( ZERO, ren+ Y3n(:) -p_qn* Y3c(:))
  rep = max( ZERO, rep+ Y3p(:) -p_qp* Y3c(:))

  call flux_vector( iiBen, ppY3n,ppY3n, -ren * p_pePel)
  call flux_vector( iiBen, ppY3p,ppY3p, -rep * p_pePel)
  call flux_vector( iiBen, ppY3n,ppK4n,  ren * (ONE-p_pePel))
  call flux_vector( iiBen, ppY3p,ppK1p,  rep * (ONE-p_pePel))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Add respiration and excretion to the benthic totals
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrBTo(:)  =   rrBTo(:)+ rrc/ MW_C * ( ONE-p_pePel)
  reBTn(:)  =   reBTn(:)+ ren * ( ONE-p_pePel)
  reBTp(:)  =   reBTp(:)+ rep * ( ONE-p_pePel)

#ifdef INCLUDE_BENCO2
  jbotO3c(:)=jbotO3c(:)+rrc*p_pePel
#endif
  jbotO2o(:)=jbotO2o(:)-rrc/MW_C *p_pePel
  jbotN4n(:)=jbotN4n(:)+ren *p_pePel
  jbotN1p(:)=jbotN1p(:)+rep *p_pePel

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux from Suspension feeders to Q6:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppY3c,ppQ6c, retQ6c )
  call flux_vector( iiBen, ppY3n,ppQ6n, retQ6n )
  call flux_vector( iiBen, ppY3p,ppQ6p, retQ6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to uptake of Benthic detritus in distribution of
  ! state variables (Dx.m is an undetermined source)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector(iiBen, ppD6m,ppD6m,( cmm- D6m(:))*( retQ6c- ruQ6c)/ Q6c(:))
  call flux_vector(iiBen, ppD7m,ppD7m,( cmm- D7m(:))*( retQ6n- ruQ6n)/ Q6n(:))
  call flux_vector(iiBen, ppD8m,ppD8m,( cmm- D8m(:))*( retQ6p- ruQ6p)/ Q6p(:))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total flux from Suspension feeders to R6:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiBen, ppY3c,ppY3c,- retR6c )
  call flux_vector( iiBen, ppY3n,ppY3n,- retR6n )
  call flux_vector( iiBen, ppY3p,ppY3p,- retR6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate NET flux from R6 to Suspension feeders :
  ! (can be negative!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  jRIY3c(:)  =   ruR6c  - retR6c
  jRIY3n(:)  =   ruR6n  - retR6n
  jRIY3p(:)  =   ruR6p  - retR6p
  ! The ruR6s which is uptaken by filter feeders is directly relased back 
  ! to R6: net food flux  from Y3 to/from R6 is 0
  jRIY3s(:)  =  ZERO      - ruPIs

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! excretion of food orginating from the Pelagic realm is also a
  ! sedimentation from from pelagic to benthic, and thus is 
  ! added to the total benthic boundary flux
  ! jbot< 0 : flux out of the system
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  jbotR6c(:)  =   jbotR6c(:)- retR6c* (ONE-p_pR6Pel)
  jbotR6n(:)  =   jbotR6n(:)- retR6n* (ONE-p_pR6Pel)
  jbotR6p(:)  =   jbotR6p(:)- retR6p* (ONE-p_pR6Pel)

  ! The silicate is directly transferred to Q6.s
  ! the ruPis which is put back in R6 is however sedimentating:
  jbotR6s(:)  =   jbotR6s(:)- ruPIs- ruR6s -reR6s


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! pseudo faeces production
  ! This production lead only to a flux to the sediment!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( sw_uptake /= 1) then
    
    r =  min( p_Rps * et * eo * p_vum * fsat* Y3c(:) *RTc,0.25*RI_Fc(:))
    r =  r/(p_small + RI_Fc(:))
 

    reR6c=  max(ZERO,r * RI_Fc(:) -ruR6c)
    reR6n=  max(ZERO,r * RI_Fn(:) -ruR6n)
    reR6p=  max(ZERO,r * RI_Fp(:) -ruR6p)
    reR6s=  max(ZERO,r * RI_Fs(:) -ruR6s)

    jbotR6c(:)  =   jbotR6c(:)- reR6c
    jbotR6n(:)  =   jbotR6n(:)- reR6n
    jbotR6p(:)  =   jbotR6p(:)- reR6p

  endif
#endif

  end subroutine FilterFeederDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
