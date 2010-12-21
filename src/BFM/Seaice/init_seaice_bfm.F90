#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise seaice BFM variables
!
! !INTERFACE:
   subroutine init_seaice_bfm(namlst,fname,unit,setup)
!
! !DESCRIPTION:
!  Initialize seaice variables from namelist. 
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
   use mem_Seaicealgae, ONLY: p_qnRc,p_qpRc,p_qsRc,p_qchlcSI
   use mem_Param,  ONLY:  p_small, &
                          CalcSeaiceBacteria
   use mem_Param,  ONLY: CalcSeaiceAlgae, CalcSeaiceZoo
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
   integer   :: icontrol,i,j,iiLastElement

   REALTYPE  :: S1c0, S1n0, S1p0, S1s0, S1l0, S2c0, S2n0, S2p0, S2l0, &
                U1c0, U1n0, U1p0, U6c0, U6n0, U6p0, U6s0, &
                X1c0, X1n0, X1p0,        &
                T1c0, T1n0, T1p0,        &
                I1p0, I3n0, I4n0, I5s0, F2o0, F3c0

  namelist /bfm_seaice_init_nml/  S1c0, S1n0, S1p0, S1s0, S2c0, S2n0, S2p0, &
                U1c0, U1n0, U1p0, U6c0, U6n0, U6p0, U6s0, &
                X1c0, X1n0, X1p0,        &
                T1c0, T1n0, T1p0,        &
                I1p0, I3n0, I4n0, I5s0, F2o0, F3c0

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
   LEVEL2 'init_seaice_bfm: done also if pelagic setup only'
   !---------------------------------------------
   ! Give reasonable initial values
   ! Overwritten by namelist parameters
   !---------------------------------------------

   !---------------------------------------------
   ! Seaice variables
   !---------------------------------------------
   S1c0  = _ONE_
   S1n0  = _ZERO_
   S1p0  = _ZERO_
   S1s0  = _ZERO_
   S1l0  = _ZERO_
   S2c0  = _ONE_
   S2n0  = _ZERO_
   S2p0  = _ZERO_
   S2l0  = _ZERO_
   U1c0  = _ONE_
   U1n0  = _ZERO_
   U1p0  = _ZERO_
   U6c0  = _ONE_
   U6n0  = _ZERO_
   U6p0  = _ZERO_
   U6s0  = _ZERO_
   X1c0  = _ONE_
   X1n0  = _ZERO_
   X1p0  = _ZERO_
   T1c0  = _ONE_
   T1n0  = _ZERO_
   T1p0  = _ZERO_
   I1p0  = _ONE_
   I3n0  = _ONE_
   I4n0  = _ONE_
   I5s0  = _ONE_
   F2o0  = _ONE_
   F3c0  = _ONE_

   !---------------------------------------------
   ! Open and read the namelist
   !---------------------------------------------
   icontrol=0
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bfm_seaice_init_nml,err=101)
   close(namlst)
   icontrol=1
98 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(fname)
     LEVEL3 'The initial values of the BFM variables are set to ONE'
     LEVEL3 'If thats not what you want you have to supply ',trim(fname)
   end if

   !---------------------------------------------
   ! Initialise seaice state variables
   !---------------------------------------------
   !MAV: need to always give initial non-zero values
   ! because there are still part of the
   ! seaice system which are computed when setup=1
   S1c  = S1c0  
   S1n  = S1n0
   S1p  = S1p0
   S1s  = S1s0
   S1l  = S1l0
   S2c  = S2c0
   S2n  = S2n0
   S2p  = S2p0
   S2l  = S2l0
   U1c  = U1c0
   U1n  = U1n0
   U1p  = U1p0
   U6c  = U6c0
   U6n  = U6n0
   U6p  = U6p0
   U6s  = U6s0
   X1c  = X1c0 
   X1n  = X1n0 
   X1p  = X1p0
   T1c  = T1c0 
   T1n  = T1n0 
   T1p  = T1p0
   I1p  = I1p0
   I3n  = I3n0
   I4n  = I4n0
   I5s  = I5s0
   F2o  = F2o0
   F3c  = F3c0

   !---------------------------------------------
   ! Initialise organisms' internal components
   ! with Redfield
   !---------------------------------------------
   do i = 1,iiSeaiceAlgae
         if (ppSeaiceAlgae(i,iiS)>0) then
            call init_cnps(c=SeaiceAlgae(i,iiC),              &
                           n=D2STATE(ppSeaiceAlgae(i,iiN),:), &
                           p=D2STATE(ppSeaiceAlgae(i,iiP),:), &
                           s=D2STATE(ppSeaiceAlgae(i,iiS),:), &
                           l=D2STATE(ppSeaiceAlgae(i,iiL),:), &
                           lc=p_qchlcSI(i), nc=p_qnRc(i),         &
                           pc=p_qpRc(i),  sc=p_qsRc(i))
         else
            call init_cnps(c=SeaiceAlgae(i,iiC),              &
                           n=D2STATE(ppSeaiceAlgae(i,iiN),:), &
                           p=D2STATE(ppSeaiceAlgae(i,iiP),:), &
                           l=D2STATE(ppSeaiceAlgae(i,iiL),:), &
                           lc=p_qchlcSI(i), nc=p_qnRc(i), pc=p_qpRc(i))
         end if
   end do

   ! Initialise bacteria components 
   do i = 1,iiSeaiceBacteria
      call init_cnps(c=SeaiceBacteria(i,iiC),  &
                     n=D2STATE(ppSeaiceBacteria(i,iiN),:), & 
                     p=D2STATE(ppSeaiceBacteria(i,iiP),:))
   end do

   ! Initialise zooplankton components checking for fixed-quota
   do i = 1 , ( iiSeaiceZoo)
      if ( (ppSeaiceZoo(i,iiP)>0) .and. (ppSeaiceZoo(i,iiN)>0)) &
            call init_cnps(c=SeaiceZoo(i,iiC),  &
                           n=D2STATE(ppSeaiceZoo(i,iiN),:), &
                           p=D2STATE(ppSeaiceZoo(i,iiP),:))
   end do

   !---------------------------------------------
   ! Initialise detritus' components  
   !---------------------------------------------
   do i = 1,iiSeaiceDetritus
      if (ppSeaiceDetritus(i,iiS)>0) then
      call init_cnps(c=SeaiceDetritus(i,iiC),     &
            n=D2STATE(ppSeaiceDetritus(i,iiN),:), &
            p=D2STATE(ppSeaiceDetritus(i,iiP),:), &
            s=D2STATE(ppSeaiceDetritus(i,iiS),:))
      else
         call init_cnps(c=SeaiceDetritus(i,iiC),     &
            n=D2STATE(ppSeaiceDetritus(i,iiN),:), &
            p=D2STATE(ppSeaiceDetritus(i,iiP),:))
      end if
   end do

   !---------------------------------------------
   ! Zeroing of the switched off state variables
   !---------------------------------------------
   do j = 1,iiSeaiceAlgae
      iiLastElement = iiL
      if (.NOT.CalcSeaiceAlgae(j)) then
         if (j==iiP1) iiLastElement=iiS
         do i = iiC,iiLastElement
            D2STATE(ppSeaiceAlgae(j,i),:) = p_small
         end do
      end if
   end do

   do j = 1,iiSeaiceZoo
      iiLastElement = iiP
      if (.NOT.CalcSeaiceZoo(j)) then
         do i = iiC,iiLastElement
            D2STATE(ppSeaiceZoo(j,i),:) = p_small
         end do
      end if
   end do

   do j = 1,iiSeaiceBacteria
        iiLastElement = iiP
        if (.NOT.CalcSeaiceBacteria(j)) then
            do i = iiC,iiLastElement
            D2STATE(ppSeaiceBacteria(j,i),:) = p_small
            end do
        end if
    end do

   return

101  FATAL 'I could not read bfm_seaice_init_nml'
    stop 'init_seaice_bfm'

   end subroutine init_seaice_bfm
!EOC

