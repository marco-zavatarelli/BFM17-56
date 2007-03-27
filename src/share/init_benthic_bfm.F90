#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise benthic BFM variables
!
! !INTERFACE:
   subroutine init_benthic_bfm(namlst,fname,unit,setup)
!
! !DESCRIPTION:
!  Initialize benthic variables from namelist. 
!  Benthic nutrients are initialized either from namelist or
!  by assuming equlibrium conditions with water column concentrations
!
! !USES:
#ifndef NOT_STANDALONE
   use api_bfm
   use global_mem
#endif
#ifdef BFM_GOTM
   use bio_var
   use bio_bfm
#endif
   use mem
   use mem_BenBac, ONLY: p_qnc,p_qpc
   use mem_Param,  ONLY: CalcBenthicFlag, p_small, &
                        CalcBenOrganisms, CalcBenBacteria
   use string_functions, ONLY: getseq_number,empty
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: namlst
   character(len=*), intent(in)        :: fname
   integer,          intent(in)        :: unit
   integer,          intent(in)        :: setup
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi & Piet Ruardij
!
! !LOCAL VARIABLES:
   integer   :: icontrol,i,j,iiLastElement,n
   REALTYPE  :: Y1c0, Y2c0, Y3c0, Y4c0, Y5c0,       &
                Q1c0, q11c0, Q6c0, K3n0, G4n0,      &
                H1c0, H2c0, K1p0, k11p0, k21p0,     &
                K4n0, k14n0, k24n0, K6r0,K5s0,      &
                D1m0, D2m0, D6m0, D7m0, D8m0, D9m0, &
                G2o0, p_qpQIc,p_qnQIc,p_qsQIc
  namelist /bfm_ben_init_nml/  calc_init_bennut_states,    &
                           p_qpQIc,p_qnQIc,p_qsQIc,          &
                           Y1c0, Y2c0, Y3c0, Y4c0, Y5c0,     &
                           Q1c0, q11c0, Q6c0, K3n0, G4n0,    &
                           H1c0, H2c0, K1p0, k11p0, k21p0,   &
                           K4n0, k14n0, k24n0, K6r0,K5s0,    &
                           D1m0, D2m0, D6m0, D7m0, D8m0, D9m0, G2o0
   interface
      subroutine init_cnps(c,n,p,s,l,nc,pc,sc,lc)
         REALTYPE,dimension(:),intent(in)             :: c
         REALTYPE,intent(in),optional                 :: nc,pc,sc,lc
         REALTYPE,dimension(:),intent(inout),optional :: n,p,s,l
      end subroutine init_cnps
   end interface
! COPYING
!
!   Copyright (C) 2006 P. Ruardij and Marcello Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_benthic_bfm: done also if pelagic setup only'
   !---------------------------------------------
   ! Give reasonable initial values
   ! Overwritten by namelist parameters
   !---------------------------------------------

   !---------------------------------------------
   ! Benthic variables
   !---------------------------------------------
   Y1c0  = _ONE_
   Y2c0  = _ONE_
   Y3c0  = _ONE_
   Y4c0  = _ONE_
   Y5c0  = _ONE_
   Q1c0  = _ONE_
   Q11c0 = _ONE_
   Q6c0  = _ONE_
   H1c0  = _ONE_
   H2c0  = _ONE_
   K1p0  = _ONE_
   K11p0 = _ONE_
   K21p0 = _ONE_
   K3n0  = _ONE_
   G4n0  = _ONE_
   K4n0  = _ONE_
   K14n0 = _ONE_
   K24n0 = _ONE_
   K6r0  = _ONE_
   K5s0  = _ONE_
   D1m0  = _ONE_
   D2m0  = _ONE_
   D6m0  = _ONE_
   D7m0  = _ONE_
   D8m0  = _ONE_
   D9m0  = _ONE_
   G2o0  = _ONE_
   p_qpQIc = -_ONE_
   p_qnQIc = -_ONE_
   p_qsQIc = -_ONE_

   !---------------------------------------------
   ! Open and read the namelist
   !---------------------------------------------
   icontrol=0
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bfm_ben_init_nml,err=101)
   close(namlst)
   icontrol=1
98 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(fname)
     LEVEL3 'The initial values of the BFM variables are set to ONE'
     LEVEL3 'If thats not what you want you have to supply ',trim(fname)
   end if

   !---------------------------------------------
   ! Check benthic model
   !---------------------------------------------
   select case (CalcBenthicFlag)
     case (0)
        LEVEL3 "Benthic model is: not used"
     case (1)
        LEVEL3 "Benthic model is: simple nutrient return"
     case (2)
        LEVEL3 "Benthic model is: benthos + intermediate nutrient return"
     case (3)
        LEVEL3 "Benthic model is: benthos + Ruardij & Van Raaphorst"
   end select

   !---------------------------------------------
   ! Write defined variables to stdout
   !---------------------------------------------
   if (bio_setup >= 2) then
      LEVEL3 'Benthic variables:'
      do n=stBenStateS,stBenStateE
        LEVEL4 trim(var_names(n)),'  ',trim(var_units(n)), &
           '  ',trim(var_long(n))
      end do
   end if

   !---------------------------------------------
   ! Initialise benthic state variables
   !---------------------------------------------
   !MAV: need to always give initial non-zero values
   ! because there are still part of the
   ! benthic system which are computed when setup=1
      Y1c  = Y1c0
      Y2c  = Y2c0
      Y3c  = Y3c0
      Y4c  = Y4c0
      Y5c  = Y5c0
      Q1c  = Q1c0
      Q11c = Q11c0
      Q6c  = Q6c0
      H1c  = H1c0
      H2c  = H2c0
      K1p  = K1p0
      K11p = K11p0
      K21p = K21p0
      K3n  = K3n0
      G4n  = G4n0
      K4n  = K4n0
      K14n = K14n0
      K24n = K24n0
      K6r  = K6r0
      K5s  = K5s0
      D1m  = D1m0
      D2m  = D2m0
      D6m  = D6m0
      D7m  = D7m0
      D8m  = D8m0
      D9m  = D9m0
      G2o  = G2o0

      !---------------------------------------------
      ! Initialise organisms' internal components
      ! with Redfield
      !---------------------------------------------
      call init_cnps(c=Y1c,n=Y1n,p=Y1p)
      call init_cnps(c=Y2c,n=Y2n,p=Y2p)
      call init_cnps(c=Y3c,n=Y3n,p=Y3p)
      call init_cnps(c=Y4c,n=Y4n,p=Y4p)
      call init_cnps(c=Y5c,n=Y5n,p=Y5p)
      call init_cnps(c=H1c,n=H1n,p=H1p,nc=p_qnc(iiH1), &
           pc=p_qpc(iiH1))
      call init_cnps(c=H2c,n=H2n,p=H2p,nc=p_qnc(iiH2), &
           pc=p_qpc(iiH2))

      !---------------------------------------------
      ! Initialise detritus' components  with Redfield
      !---------------------------------------------
      call init_cnps(c=Q1c,n=Q1n,p=Q1p,nc=p_qnQIc,pc=p_qpQIc)
      call init_cnps(c=Q11c,n=Q11n,p=Q11p,nc=p_qnQIc,pc=p_qpQIc)
      call init_cnps(c=Q6c,n=Q6n,p=Q6p,s=Q6s,nc=p_qnQIc,pc=p_qpQIc,sc=p_qsQIc)

      !---------------------------------------------
      ! Initialise benthic nutrients from water column
      ! conditions
      !---------------------------------------------
#ifndef BFM_GETM
      if ( calc_init_bennut_states == 1) then
           LEVEL4 "Benthic nutrient variables are initialised by assuming"
           LEVEL4 "equilibrium between input (nutrient regeneratation) and output"
           LEVEL4 "(flux to water column and definitive loss processes)"
           call InitBenthicNutrientDynamics
      end if
#endif

   !---------------------------------------------
   ! Zeroing of the switched off state variables
   !---------------------------------------------
   do j = 1,iiBenOrganisms
      iiLastElement = iiP
      if (.NOT.CalcBenOrganisms(j)) then
         do i = iiC,iiLastElement
            D2STATE(ppBenOrganisms(j,i),:) = p_small
         end do
      end if
   end do
   do j = 1,iiBenBacteria
      iiLastElement = iiP
      if (.NOT.CalcBenBacteria(j)) then
         do i = iiC,iiLastElement
            D2STATE(ppBenBacteria(j,i),:) = p_small
         end do
      end if
   end do

   return

101  FATAL 'I could not read bfm_ben_init_nml'
    stop 'init_benthic_bfm'

   end subroutine init_benthic_bfm
!EOC

