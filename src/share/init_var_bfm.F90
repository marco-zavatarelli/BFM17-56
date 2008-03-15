#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise BFM variables
!
! !INTERFACE:
   subroutine init_var_bfm(namlst,fname,unit,setup)
!
! !DESCRIPTION:
!  Allocate BFM variables and give initial values of
!  parameters and state variables
!  Only pelagic variables are initialized here.
!  Benthic variables are done in a special routine init_benthic_bfm
!
! !USES:
#ifndef NOT_STANDALONE
   use api_bfm
   use global_mem
#endif
#ifdef BFM_GOTM
   use bio_var
   use bio_bfm
   use global_mem, ONLY: RLEN,ZERO,ONE
#endif
   use mem
   use mem_Phyto, ONLY: p_qnRc,p_qpRc,p_qsRc
   use constants, ONLY: HOURS_PER_DAY
   use mem_Param, ONLY: CalcPelagicFlag,CalcBenthicFlag,p_small,p_qchlc, &
                        CalcPhytoPlankton,CalcMicroZooPlankton,          &
                        CalcMesoZooPlankton,CalcPelChemistry
#ifdef INCLUDE_BEN
   use mem_Param, ONLY: CalcBenOrganisms, CalcBenBacteria, CalcBacteria
   use mem_BenBac, ONLY: p_qnc,p_qpc
#endif
   use mem_Param, ONLY: AssignPelBenFluxesInBFMFlag
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
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer              :: icontrol,i,j,iiLastElement,n
   integer,parameter    :: NSAVE=100  ! Maximum no variables which can be saved
   character(len=64),dimension(NSAVE):: var_save
   character(len=64),dimension(NSAVE):: ave_save
   REALTYPE  :: N1p0,N3n0,N4n0,N5s0,N6r0,  &
                P1c0,P2c0,P3c0,P4c0,Z3c0,  &
                Z4c0,Z5c0,Z6c0,B1c0,R1c0,  &
                R2c0,R6c0,R7c0,O2o0,O4n0,  &
#ifdef INCLUDE_PELCO2
                O3c0,O3h0,                 &
#endif
                P1l0,P2l0,P3l0,P4l0,       &
                P1n0,P2n0,P3n0,P4n0,       &
                P1p0,P2p0,P3p0,P4p0,P1s0

   namelist /bfm_init_nml/ surface_flux_method,       &
                           n_surface_fluxes,          &
                           N1p0,N3n0,N4n0,N5s0,N6r0,  &
                           P1c0,P2c0,P3c0,P4c0,Z3c0,  &
                           Z4c0,Z5c0,Z6c0,B1c0,R1c0,  &
                           R2c0,R6c0,R7c0,O2o0,O4n0,  &
#ifdef INCLUDE_PELCO2
                           O3c0,O3h0,                 &
#endif
                           P1l0,P2l0,P3l0,P4l0,       &
                           P1n0,P2n0,P3n0,P4n0,       &
                           P1p0,P2p0,P3p0,P4p0,P1s0

   namelist /bfm_save_nml/ var_save, ave_save

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

   LEVEL2 'init_var_bfm'
   !---------------------------------------------
   ! Give zero initial values
   ! Overwritten by namelist parameters
   !---------------------------------------------
   surface_flux_method = -1
   n_surface_fluxes = 1

   !---------------------------------------------
   ! Pelagic variables
   !---------------------------------------------
   N1p0 = _ONE_
   N3n0 = _ONE_
   N4n0 = _ONE_
   N5s0 = _ONE_
   N6r0 = _ONE_
   O2o0 = 280.0_RLEN
#ifdef INCLUDE_PELCO2
   O3c0 = 24785.0_RLEN ! conversion from GLODAP 2303 umol/kg to mg/m3
   O3h0 = 2303.0_RLEN  ! GLODAP surface average (umol/kg)
#endif
   O4n0 = _ONE_
   P1c0 = _ZERO_
   P2c0 = _ZERO_
   P3c0 = _ZERO_
   P4c0 = _ZERO_
   P1l0 = _ZERO_
   P2l0 = _ZERO_
   P3l0 = _ZERO_
   P4l0 = _ZERO_
   P1n0 = _ZERO_
   P2n0 = _ZERO_
   P3n0 = _ZERO_
   P4n0 = _ZERO_
   P1p0 = _ZERO_
   P2p0 = _ZERO_
   P3p0 = _ZERO_
   P4p0 = _ZERO_
   P1s0 = _ZERO_
   Z3c0 = _ZERO_
   Z4c0 = _ZERO_
   Z5c0 = _ZERO_
   Z6c0 = _ZERO_
   B1c0 = _ZERO_
   R1c0 = _ZERO_
   R2c0 = _ZERO_
   R6c0 = _ZERO_
   R7c0 = _ZERO_

   !---------------------------------------------
   ! Open and read the namelist
   !---------------------------------------------
   icontrol=0
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bfm_init_nml,err=99)
   var_save=""
   ave_save=""
   var_ave=.false.
   read(namlst,nml=bfm_save_nml,err=100)
   close(namlst)
   icontrol=1
98 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(fname)
     LEVEL3 'The initial values of the BFM variables are set to ONE'
     LEVEL3 'If thats not what you want you have to supply ',trim(fname)
   end if

   !---------------------------------------------
   ! Check variable to be saved and
   ! set the corresponding flag value in var_ids
   !---------------------------------------------
   do i=1,NSAVE
      if (.NOT.empty(var_save(i))) then
            j=getseq_number(var_save(i),var_names,stBenFluxE,.TRUE.)
            if ( j > 0 ) var_ids(j)=-1
      end if
      if ( .NOT.empty(var_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(var_save(i)),' does not exist!'
      end if
   end do
   do i=1,NSAVE
      if (.NOT.empty(ave_save(i))) then
         j=getseq_number(ave_save(i),var_names,stBenFluxE,.TRUE.)
         if ( .NOT.empty(ave_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(ave_save(i)),' does not exist!'
         else if ( var_ids(j) <0 ) then
            STDERR 'Warning: Variable ',trim(ave_save(i)), &
               ' is already selected for output in var_save'
         else if ( j > 0 ) then
            var_ids(j)=-1
            var_ave(j)=.true.
         end if
      end if
   end do

#ifdef BFM_GOTM
   !---------------------------------------------
   ! Create pointers
   !---------------------------------------------
    call pointers_gotm_bfm()
#endif

   !---------------------------------------------
   ! Initialize BFM parameters
   !---------------------------------------------
   call Initialize

   !---------------------------------------------
   ! Initially set the number of sun hours
   ! equal to the number of hours in a day.
   !---------------------------------------------
   SUNQ = HOURS_PER_DAY

   !---------------------------------------------
   ! Initialise pelagic state variables
   ! also if using a benthic-only setup
   ! (for boundary conditions)
   !---------------------------------------------
      N1p = N1p0
      N3n = N3n0
      N4n = N4n0
      N5s = N5s0
      N6r = N6r0
      O2o = O2o0
#ifdef INCLUDE_PELCO2
      O3c = O3c0
      O3h = O3h0
#endif
      O4n = O4n0
      P1c = P1c0
      P1n = P1n0
      P1p = P1p0
      P1l = P1l0
      P1s = P1s0
      P2c = P2c0
      P2n = P2n0
      P2p = P2p0
      P2l = P2l0
      P3c = P3c0
      P3n = P3n0
      P3p = P3p0
      P3l = P3l0
      P4c = P4c0
      P4n = P4n0
      P4p = P4p0
      P4l = P4l0
      Z3c = Z3c0
      Z4c = Z4c0
      Z5c = Z5c0
      Z6c = Z6c0
      B1c = B1c0
      R1c = R1c0
      R2c = R2c0
      R6c = R6c0
      R7c = R7c0

      !---------------------------------------------
      ! Initialise other internal components
      ! with Redfield
      !---------------------------------------------
      do i = 1 , ( iiPhytoPlankton)
         if (ppPhytoPlankton(i,iiS)>0) then
            call init_cnps(c=PhytoPlankton(i,iiC),              &
                           n=D3STATE(ppPhytoPlankton(i,iiN),:), &
                           p=D3STATE(ppPhytoPlankton(i,iiP),:), &
                           s=D3STATE(ppPhytoPlankton(i,iiS),:), &
                           l=D3STATE(ppPhytoPlankton(i,iiL),:), &
                           lc=p_qchlc(i), nc=p_qnRc(i),         &
                           pc=p_qpRc(i),  sc=p_qsRc(i))
         else
            call init_cnps(c=PhytoPlankton(i,iiC),              &
                           n=D3STATE(ppPhytoPlankton(i,iiN),:), &
                           p=D3STATE(ppPhytoPlankton(i,iiP),:), &
                           l=D3STATE(ppPhytoPlankton(i,iiL),:), &
                           lc=p_qchlc(i), nc=p_qnRc(i), pc=p_qpRc(i))
         end if
      end do
      call init_cnps(c=B1c,n=B1n,p=B1p)
      call init_cnps(c=R1c,n=R1n,p=R1p)
      call init_cnps(c=R6c,n=R6n,p=R6p,s=R6s)
      ! Initialise zooplankton components checking for fixed-quota
      do i = 1 , ( iiMicroZooPlankton)
         if ( (ppMicroZooPlankton(i,iiP)>0) .and. (ppMicroZooPlankton(i,iiN)>0) ) &
            call init_cnps(c=MicroZooPlankton(i,iiC),  &
                           n=D3STATE(ppMicroZooPlankton(i,iiN),:), &
                           p=D3STATE(ppMicroZooPlankton(i,iiP),:))
      end do
      do i = 1 , ( iiMesoZooPlankton)
         if ( (ppMesoZooPlankton(i,iiP) > 0) .and. (ppMesoZooPlankton(i,iiN)>0) ) &
            call init_cnps(c=MesoZooPlankton(i,iiC),  &
                           n=D3STATE(ppMesoZooPlankton(i,iiN),:), &
                           p=D3STATE(ppMesoZooPlankton(i,iiP),:))
      end do

   !---------------------------------------------
   ! Check setup settings
   ! and finalize initialization
   !---------------------------------------------
   select case (setup)
      case (0)
      case (1) ! Pelagic only
         LEVEL2 "Pelagic-only setup (bio_setup=1), Switching off the benthic system"
         CalcBenthicFlag = 0
#ifndef BFM_GOTM
         ! force computation of bottom fluxes in the BFM
         AssignPelBenFluxesInBFMFlag = .TRUE.
#endif
#ifdef INCLUDE_BEN
      case (2) ! Benthic only
         LEVEL2 "Benthic-only setup (bio_setup=2), Switching off the pelagic system"
         CalcPelagicFlag = .FALSE.
         CalcPhytoPlankton=.FALSE.
         CalcBacteria=.FALSE.
         CalcMesoZooPlankton=.FALSE.
         CalcMicroZooPlankton=.FALSE.
      case (3) ! Pelagic-Benthic coupling
         LEVEL2 "Pelagic-Benthic coupled setup (bio_setup=3)"
         if (CalcBenthicFlag == 0) &
            LEVEL3 'Warning, benthic system is switched off!'
         if (.NOT.CalcPelagicFlag) &
            LEVEL3 'Warning, pelagic system is switched off!'
#endif
   end select

   !---------------------------------------------
   ! Write defined variables to stdout
   !---------------------------------------------
   if (setup /= 2) then
      LEVEL3 'Pelagic variables:'
      do n=stPelStateS,stPelStateE
         LEVEL4 trim(var_names(n)),'  ',trim(var_units(n)), &
           '  ',trim(var_long(n))
      end do
   endif

   !---------------------------------------------
   ! Zeroing of the switched off state variables
   !---------------------------------------------
   do j = 1,iiPhytoPlankton
      iiLastElement = iiL
      if (.NOT.CalcPhytoPlankton(j)) then
         if (j==iiP1) iiLastElement=iiS
         do i = iiC,iiLastElement
            D3STATE(ppPhytoPlankton(j,i),:) = p_small
            D3STATETYPE(ppPhytoPlankton(j,i)) = NOTRANSPORT
         end do
      end if
   end do
   do j = 1,iiMesoZooPlankton
      if ( (ppMesoZooPlankton(j,iiP)>0) .and. &
           (ppMesoZooPlankton(j,iiN)>0) ) then
         iiLastElement = iiP
      else
         iiLastElement = iiC
      end if
      if (.NOT.CalcMesoZooPlankton(j)) then
         do i = iiC,iiLastElement
            D3STATE(ppMesoZooPlankton(j,i),:) = p_small
            D3STATETYPE(ppMesoZooPlankton(j,i)) = NOTRANSPORT
         end do
      end if
   end do
   do j = 1,iiMicroZooPlankton
      if ( (ppMicroZooPlankton(j,iiP)>0) .and. &
           (ppMicroZooPlankton(j,iiN)>0) ) then
         iiLastElement = iiP
      else
         iiLastElement = iiC
      end if
      if (.NOT.CalcMicroZooPlankton(j)) then
         do i = iiC,iiLastElement
            D3STATE(ppMicroZooPlankton(j,i),:) = p_small
            D3STATETYPE(ppMicroZooPlankton(j,i)) = NOTRANSPORT
         end do
      end if
   end do
#ifdef INCLUDE_BEN
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
   if (.NOT.CalcBacteria) then
      B1c = p_small; B1n = p_small; B1p = p_small;
      D3STATETYPE(ppB1c) = NOTRANSPORT
      D3STATETYPE(ppB1n) = NOTRANSPORT
      D3STATETYPE(ppB1p) = NOTRANSPORT
   end if
#endif

   return

99  FATAL 'I could not read bfm_init_nml'
    stop 'init_var_bfm'
100 FATAL 'I could not read bfm_save_nml'
    stop 'init_var_bfm'

   end subroutine init_var_bfm
!EOC

