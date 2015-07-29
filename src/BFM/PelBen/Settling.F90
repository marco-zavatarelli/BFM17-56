#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Settling
!
! DESCRIPTION
!   This process describes the dynamics of sedimentation and
!    deposition of phytoplankton (Pi) and detritus (R6)
!    in the benthic system.
!    A burial velocity is defined, which controls the magnitude
!    of the inflow rate of detritus from the pelagic form R6 to
!    the benthic form Q6.
!    The processes described here takes place only in the bottom layer
!
! !INTERFACE
  subroutine SettlingDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: R1c, R2c, R6c, R1n, R6n, R1p, R6p, R6s, PhytoPlankton, D3STATE
  use mem, ONLY: ppR1c, ppR6c, ppR1n, ppR6n, ppR1p, ppR6p, &
    ppR6s, ppPhytoPlankton, &
    NO_BOXES_XY, BoxNumberXY, Depth, jbotR6c, jbotR6n, jbotR6p, &
    jbotR6s, jbotR1c, jbotR1n, jbotR1p, sediPPY, sediR2, sediR6, iiPhytoPlankton, &
    iiC, iiN, iiP, iiL, iiS, iiBen, iiPel, PELBOTTOM, flux
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF,R6f,ppR6f,jbotR6f
#endif
#endif
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p
  use mem_Settling
  use mem_PelBac, ONLY: p_suhR1,p_sulR1,p_suR2,p_suR6, &
                        p_qncPBA,p_qpcPBA
#ifdef BFM_GOTM
  use bio_var, ONLY: BOTindices
#else
  use api_bfm, ONLY: BOTindices
#endif

!
!
! !AUTHORS
!   the BFM team
!
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Vectors used to change the scope of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
  real(RLEN), dimension(:), pointer  ::lcl_ppPhytoPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer  :: j
  integer  :: kbot
  real(RLEN)  :: sedi
  real(RLEN)  :: ruQIc
  real(RLEN)  :: ruQIn
  real(RLEN)  :: ruQIp
  real(RLEN)  :: ruQ1c
  real(RLEN)  :: ruQ1n
  real(RLEN)  :: ruQ1p
  real(RLEN)  :: ruQ6c
  real(RLEN)  :: ruQ6n
  real(RLEN)  :: ruQ6p
  real(RLEN)  :: ruQIl
  real(RLEN)  :: ruQ6s
#ifdef INCLUDE_PELFE
  real(RLEN)  :: ruQ6f
#endif
  real(RLEN)  :: s
  real(RLEN)  :: p

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Starting physical settling:
   ! phytoplankton settling are temporarely stored as R1 and R6 and
   ! subsequently exported (BentoPelCoup) and incorporated in 
   ! the benthic system (sedimentation)
   ! Organic material thus enters the sediment with a smaller velocity
   ! than the sedimentation velocity along the water column
   ! This burial rate can be larger than the physical burial rate in case
   ! some material are particularly sticky
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! loop over the number of bottom boxes
      do BoxNumberXY = 1,NO_BOXES_XY 
         kbot = BOTindices(BoxNumberXY)
         do i = 1 , ( iiPhytoPlankton)
            sedi  =   sediPPY(i,kbot)
            if ( sedi> ZERO ) then
               ! Phytoplankton carbon
               j=ppPhytoPlankton(i,iiC)
               lcl_PhytoPlankton => PhytoPlankton(i,iiC)
               ruQIc  =   sedi* lcl_PhytoPlankton(kbot)
               ruQ1c  =   p_pe_R1c* ruQIc
               ruQ6c  =   ruQIc- ruQ1c
               PELBOTTOM(j,BoxNumberXY)=  - ruQIc
               jbotR1c(BoxNumberXY)  =   jbotR1c(BoxNumberXY)- ruQ1c
               jbotR6c(BoxNumberXY)  =   jbotR6c(BoxNumberXY)- ruQ6c
               ! Phytoplankton nitrogen
               j=ppPhytoPlankton(i,iiN)
               lcl_PhytoPlankton => PhytoPlankton(i,iiN)
               ruQIn  =   sedi* lcl_PhytoPlankton(kbot)
               ruQ1n  =   p_pe_R1n* ruQIn
               ruQ6n  =   ruQIn- ruQ1n
               PELBOTTOM(j,BoxNumberXY)=  - ruQIn
               jbotR1n(BoxNumberXY)  =   jbotR1n(BoxNumberXY)- ruQ1n
               jbotR6n(BoxNumberXY)  =   jbotR6n(BoxNumberXY)- ruQ6n
               ! Phytoplankton phosphorus
               j=ppPhytoPlankton(i,iiP)
               lcl_PhytoPlankton => PhytoPlankton(i,iiP)
               ruQIp  =   sedi* lcl_PhytoPlankton(kbot)
               ruQ1p  =   p_pe_R1p* ruQIp
               ruQ6p  =   ruQIp- ruQ1p
               PELBOTTOM(j,BoxNumberXY)= - ruQIp
               jbotR1p(BoxNumberXY) =jbotR1p(BoxNumberXY)- ruQ1p
               jbotR6p(BoxNumberXY) =jbotR6p(BoxNumberXY)- ruQ6p
               ! Phytoplankton chlorophyll (stored but not used)
               j=ppPhytoPlankton(i,iiL)
               lcl_PhytoPlankton => PhytoPlankton(i,iiL)
               ruQIl  =   sedi* lcl_PhytoPlankton(kbot)
               PELBOTTOM(j,BoxNumberXY)=  - ruQIl
               ! Phytoplankton silica
               j=ppPhytoPlankton(i,iiS)
               if ( j> 0) then
                  lcl_PhytoPlankton => PhytoPlankton(i,iiS)
                  ruQ6s  =   sedi* lcl_PhytoPlankton(kbot)
                  PELBOTTOM(j,BoxNumberXY)  = - ruQ6s
                  jbotR6s(BoxNumberXY)  =   jbotR6s(BoxNumberXY)- ruQ6s
               end if

            else
                PELBOTTOM(ppPhytoPlankton(i,iiC),BoxNumberXY) = ZERO
                PELBOTTOM(ppPhytoPlankton(i,iiN),BoxNumberXY) = ZERO
                PELBOTTOM(ppPhytoPlankton(i,iiP),BoxNumberXY) = ZERO
                PELBOTTOM(ppPhytoPlankton(i,iiL),BoxNumberXY) = ZERO
                j=ppPhytoPlankton(i,iiS)
                if ( j>0) PELBOTTOM(j,BoxNumberXY) = ZERO
            end if
         end do !over phytoplankton groups

         ! R2 into Q1 and Q6
         ! NOTE: ALL DETRITUS FLUXES TO THE SEDIMENT ARE DIRECTED TO Q6 VIA R6 
         if ( sediR2(kbot) > ZERO) then
           ! Calculate how the intermediate degradable R2 has to be redistrubuted 
           ! between R1 and R6 in such a way that the degradability is the same
           ! Calculate first actual degradability of LOC ( dependent of quotum NC,PC)
           ! and compare with bacterial preferencial quota (first bacteria group)
           p = min(ONE, R1n(kbot)/(R1c(kbot)* p_qncPBA(1)), &
                                    R1p(kbot)/(R1c(kbot)*p_qpcPBA(1)))
           ! Calculate actual degradability of R1
           s= p_suhR1(1)*p-p_sulR1(1)* (ONE-p)
           ! Calculate distribution factor for R2 between R1 and R6
           p=(p_suR2(1)-p_suR6(1))/(s-p_suR6(1))
           jbotR1c(BoxNumberXY)=jbotR1c(BoxNumberXY)- &
                                            p*sediR2(kbot) * R2c(kbot)
           jbotR6c(BoxNumberXY)=jbotR6c(BoxNumberXY)-&
                                    (ONE-p)*sediR2(kbot) * R2c(kbot)
         end if

         ! R6 into Q6 through burial
         ! NOTE: ALL DETRITUS FLUXES TO THE SEDIMENT ARE DIRECTED TO Q6 VIA R6 
         ruQ6c  =   sediR6(kbot)* R6c(kbot)
         ruQ6n  =   sediR6(kbot)* R6n(kbot)
         ruQ6p  =   sediR6(kbot)* R6p(kbot)
         ruQ6s  =   sediR6(kbot)* R6s(kbot)
         jbotR6c(BoxNumberXY)  =   jbotR6c(BoxNumberXY)- ruQ6c
         jbotR6n(BoxNumberXY)  =   jbotR6n(BoxNumberXY)- ruQ6n
         jbotR6p(BoxNumberXY)  =   jbotR6p(BoxNumberXY)- ruQ6p
         jbotR6s(BoxNumberXY)  =   jbotR6s(BoxNumberXY)- ruQ6s
#ifdef INCLUDE_PELFE
         ruQ6f  =   sediR6(kbot)* R6f(kbot)
         jbotR6f(BoxNumberXY)  =   jbotR6f(BoxNumberXY)- ruQ6f
#endif
      end do

  end subroutine SettlingDynamics

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
