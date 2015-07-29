#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BentoPelCoup
!
! DESCRIPTION
!   This is the coupling interface between benthic and pelagic systems.
!   This routine is called after the computation of benthic processes
!   and assign boundary fluxes
!
! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij and M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
! !INTERFACE
  subroutine BentoPelCoupDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: R6c, R6n, R6p, R6s, O2o, N1p, N3n, N4n, N5s, N6r, &
    R1c, R1n, R1p, PhytoPlankton, ppPhytoPlankton,iiPhytoPlankton,  &
    MicroZooPlankton, ppMicroZooPlankton,iiMicroZooPlankton
  use mem, ONLY: ppR6c, ppR6n, ppR6p, ppR6s, ppO2o, ppN1p, &
    ppN3n, ppN4n, ppN5s, ppN6r, ppR1c, ppR1n, ppR1p, &
    NO_BOXES_XY, ERHO, BoxNumberXY, Depth, &
    jbotO2o, jbotN1p, jbotN3n, jbotN4n, jbotN5s, jbotN6r, jbotR6c, jbotR6n, &
    jbotR6p, jbotR6s, jbotR1c, jbotR1n, jbotR1p, PELBOTTOM, &
    iiP1, iiC, iiN, iiP, iiL, iiS, iiBen, iiPel, flux
#ifdef INCLUDE_BEN
  use mem, ONLY: jPIY3c, jZIY3c, ZI_Fc, jRIY3c, jRIY3n, jRIY3p, jRIY3s
#endif
#ifdef INCLUDE_PELFE
  use mem, ONLY: iiF,R6f,ppR6f,jbotR6f,R1f,ppR1f,jbotR1f
#endif
#if defined INCLUDE_PELCO2 && defined INCLUDE_BENCO2
  use mem, ONLY: ppO3h, ppO3c,jbotO3c,jbotO3h
#endif
#endif

#ifdef BFM_GOTM
 use bio_var, ONLY: BOTindices
#else
 use api_bfm, ONLY: BOTindices
#endif
 use constants, ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
 use mem_Param, ONLY: CalcBenthicFlag
 use mem_Param, ONLY: AssignPelBenFluxesInBFMFlag, p_small

!  
!
! !AUTHORS
!   Piet Ruardij  
!
!
! !REVISION_HISTORY
!   Created at Mon Apr 19 00:08:12 CEST 2004
!
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij & M. Vichi
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
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
  real(RLEN), dimension(:), pointer  ::lcl_MicroZooPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer  :: j
  integer  :: kbot
  real(RLEN)  :: uptake
  real(RLEN)  :: Pc,Zc

#ifdef INCLUDE_BEN
  ! fluxes from pelagic to benthic organisms only when Benthic model
  ! is more than BENTHIC_RETURN
  if  (CalcBenthicFlag > BENTHIC_RETURN) then
   do BoxNumberXY = 1,NO_BOXES_XY
      kbot = BOTindices(BoxNumberXY)
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Add Fluxes to Filter-feeder from Pelagic for
      ! all phytoplankton subgroups and constituents
      ! Note: this part is vectorized
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      do i = 1 , ( iiPhytoPlankton)
         lcl_PhytoPlankton => PhytoPlankton(i,iiC)
         Pc  =   lcl_PhytoPlankton(kbot)
         ! there is a precautionary check on the amount of phytoplankton
         ! this may generate a small mass loss, but is negligible
         if ( Pc > p_small) then
            j = ppPhytoPlankton(i,iiC)
            uptake  =   jPIY3c(i,BoxNumberXY)
            PELBOTTOM(j,BoxNumberXY)  =  PELBOTTOM(j,BoxNumberXY)  -uptake
            j = ppPhytoPlankton(i,iiN)
            lcl_PhytoPlankton => PhytoPlankton(i,iiN)
            PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                              -uptake* lcl_PhytoPlankton(kbot)/ Pc
            j = ppPhytoPlankton(i,iiP)
            lcl_PhytoPlankton => PhytoPlankton(i,iiP)
            PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                              -uptake* lcl_PhytoPlankton(kbot)/ Pc
            j = ppPhytoPlankton(i,iiL)
            lcl_PhytoPlankton => PhytoPlankton(i,iiL)
            PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                              -uptake* lcl_PhytoPlankton(kbot)/ Pc
            j = ppPhytoPlankton(i,iiS)
            if ( j> 0) then
               lcl_PhytoPlankton => PhytoPlankton(i,iiS)
               PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                              -uptake* lcl_PhytoPlankton(kbot)/ Pc
            end if
#ifdef INCLUDE_PELFE
            j = ppPhytoPlankton(i,iiF)
            if ( j> 0) then
               lcl_PhytoPlankton => PhytoPlankton(i,iiF)
               PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                              -uptake* lcl_PhytoPlankton(kbot)/ Pc
            end if
#endif
         end if
      end do

      ! Microzooplankton
      ! In contrast with phytoplankton there is only one flux to the sediment!
      ! Only grazing. (see settling.F90)
      do i = 1 , ( iiMicroZooPlankton)
         lcl_MicroZooPlankton => MicroZooPlankton(i,iiC)
         Zc  =   lcl_MicroZooPlankton(kbot)
         if ( Zc > p_small) then
            j = ppMicroZooPlankton(i,iiC)
            uptake  =   jZIY3c(BoxNumberXY) * Zc/ ZI_Fc(BoxNumberXY)
            PELBOTTOM(j,BoxNumberXY)  =   -uptake
            j = ppMicroZooPlankton(i,iiN)
            if ( j> 0) then
              lcl_MicroZooPlankton => MicroZooPlankton(i,iiN)
              PELBOTTOM(j,BoxNumberXY) =  -uptake* lcl_MicroZooPlankton(kbot)/ &
                                          ZI_fc(BoxNumberXY)
            endif
            j = ppMicroZooPlankton(i,iiP)
            if ( j> 0) then
              lcl_MicroZooPlankton => MicroZooPlankton(i,iiP)
              PELBOTTOM(j,BoxNumberXY) =  -uptake* lcl_MicroZooPlankton(kbot)/ &
                                          ZI_fc(BoxNumberXY)
            endif
         else
            PELBOTTOM(ppMicroZooPlankton(i,iiC),:)  = ZERO
            j = ppMicroZooPlankton(i,iiN)
            if ( j> 0) PELBOTTOM(j,BoxNumberXY)  =  ZERO
            j = ppMicroZooPlankton(i,iiP)
            if ( j> 0) PELBOTTOM(j,BoxNumberXY)  =  ZERO
         end if
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Net detritus Fluxes to Benthic from Pelagic by Y3
      ! net flux= uptake - excretion of food : flux may be negative!
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux(kbot, iiPel, ppR6c, ppR6c, -( jRIY3c(BoxNumberXY)/ &
         Depth(kbot)) )
      call flux(kbot, iiPel, ppR6n, ppR6n, -( jRIY3n(BoxNumberXY)/ &
         Depth(kbot)) )
      call flux(kbot, iiPel, ppR6p, ppR6p, -( jRIY3p(BoxNumberXY)/ &
         Depth(kbot)) )
      call flux(kbot, iiPel, ppR6s, ppR6s, -( jRIY3s(BoxNumberXY)/ &
         Depth(kbot)) )
   end do ! loop over NO_BOXES_XY
  end if ! benthic model includes benthos
#endif

   if ( AssignPelBenFluxesInBFMFlag) then
   ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! All Fluxes to Benthic from Pelagic defined for the
   ! Pelagic State variables (done only if BFM computes them)
   ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      do BoxNumberXY = 1,NO_BOXES_XY
         kbot = BOTindices(BoxNumberXY)
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Nutrient Fluxes to Benthic from Pelagic
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call flux(kbot, iiPel, ppO2o, ppO2o, jbotO2o(BoxNumberXY)/ &
          Depth(kbot) )
        call flux(kbot, iiPel, ppN1p, ppN1p, jbotN1p(BoxNumberXY)/ &
          Depth(kbot) )
        call flux(kbot, iiPel, ppN3n, ppN3n, jbotN3n(BoxNumberXY)/ &
          Depth(kbot) )
        call flux(kbot, iiPel, ppN4n, ppN4n, jbotN4n(BoxNumberXY)/ &
          Depth(kbot) )
        call flux(kbot, iiPel, ppN5s, ppN5s, jbotN5s(BoxNumberXY)/ &
          Depth(kbot) )
        call flux(kbot, iiPel, ppN6r, ppN6r, jbotN6r(BoxNumberXY)/ &
          Depth(kbot) )
#if defined INCLUDE_PELCO2 && defined INCLUDE_BENCO2
        call flux(kbot, iiPel, ppO3c, ppO3c, jbotO3c(BoxNumberXY)/ &
          Depth(kbot) )
        ! convert the units from mmol/m3/d to umol/kg/d
        call flux(kbot, iiPel, ppO3h, ppO3h, jbotO3h(BoxNumberXY)/&
          Depth(kbot))
#endif

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! PhytoPlankton Fluxes to Benthic from Pelagic (mostly by Y3 uptake)
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        do i = 1 , ( iiPhytoPlankton)
              j = ppPhytoPlankton(i,iiC)
              call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(kbot) )
              j = ppPhytoPlankton(i,iiN)
              call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(kbot) )
              j = ppPhytoPlankton(i,iiP)
              call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(kbot) )
              j = ppPhytoPlankton(i,iiL)
              call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(kbot) )
              j = ppPhytoPlankton(i,iiS)
              if ( j> 0) then
                call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                  Depth(kbot) )
              end if
#ifdef INCLUDE_PELFE
              j = ppPhytoPlankton(i,iiF)
              if ( j> 0) then
                call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                  Depth(kbot) )
              end if
#endif
        end do

        do i = 1 , ( iiMicroZooPlankton)
              j = ppMicroZooPlankton(i,iiC)
              call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(kbot) )
              j = ppMicroZooPlankton(i,iiN)
              if ( j> 0) then
                call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                  Depth(kbot) )
              endif
              j = ppMicroZooPlankton(i,iiP)
              if ( j> 0) then
                call flux(kbot, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                  Depth(kbot) )
              endif
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Total Sedimentation flux to Benthic from Pelagic defined
        ! for the detritus state variables.
        ! (See sedimentation for definition of the fluxes of benthic &
        ! variables)
        ! !!!!!!! ALL DETRITUS FLUXES TO THE SEDIMENT ARE DIRECTED VIA R6 TO Q6 !!!!!!!!
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call flux(kbot, iiPel, ppR6c, ppR6c, ( jbotR6c(BoxNumberXY)/ &
          Depth(kbot)) )
        call flux(kbot, iiPel, ppR6n, ppR6n, ( jbotR6n(BoxNumberXY)/ &
          Depth(kbot)) )
        call flux(kbot, iiPel, ppR6p, ppR6p, ( jbotR6p(BoxNumberXY)/ &
          Depth(kbot)) )
        call flux(kbot, iiPel, ppR6s, ppR6s, ( jbotR6s(BoxNumberXY)/ &
          Depth(kbot)) )
        call flux(kbot, iiPel, ppR1c, ppR1c, ( jbotR1c(BoxNumberXY)/ &
          Depth(kbot)) )
        call flux(kbot, iiPel, ppR1n, ppR1n, ( jbotR1n(BoxNumberXY)/ &
          Depth(kbot)) )
        call flux(kbot, iiPel, ppR1p, ppR1p, ( jbotR1p(BoxNumberXY)/ &
          Depth(kbot)) )
#ifdef INCLUDE_PELFE
        call flux(kbot, iiPel, ppR1f, ppR1f, ( jbotR1f(BoxNumberXY)/ &
          Depth(kbot)) )
        call flux(kbot, iiPel, ppR6f, ppR6f, ( jbotR6f(BoxNumberXY)/ &
          Depth(kbot)) )
#endif
      end do ! loop over NO_BOXES_XY
   end if ! AssignPelBenFluxesInBFMFlag

  end subroutine  BentoPelCoupDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
