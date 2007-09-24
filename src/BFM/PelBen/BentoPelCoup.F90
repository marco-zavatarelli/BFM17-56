#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BentoPelCoup
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BentoPelCoupDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: R6c, R6n, R6p, R6s, &
  ! O2o, N1p, N3n, N4n, N5s, N6r, R1c, R1n, R1p
  ! For the following Pelagic-group-states fluxes are defined: PhytoPlankton
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY
  ! The following Pelagic 1-d global boxvars  are used: Depth
  ! The following Benthic 1-d global boxvars are used: jPIY3c, &
  ! PIc, jRIY3c, jRIY3n, jRIY3p, jRIY3s, jbotO2o, jbotN1p, jbotN3n, jbotN4n, &
  ! jbotN5s, jbotN6r, jbotR6c, jbotR6n, jbotR6p, jbotR6s, jbotR1c, jbotR1n, jbotR1p
  ! The following Benthic 2-d global boxvars are modified : PELBOTTOM
  ! The following groupmember vars  are used: iiPhytoPlankton, iiP1
  ! The following constituent constants  are used: iiC, iiN, iiP, iiL, iiS
  ! The following 0-d global parameters are used: &
  ! AssignPelBenFluxesInBFMFlag
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem, ONLY: R6c, R6n, R6p, R6s, O2o, N1p, N3n, N4n, N5s, N6r, &
    R1c, R1n, R1p, PhytoPlankton, D2STATE
  use mem, ONLY: ppR6c, ppR6n, ppR6p, ppR6s, ppO2o, ppN1p, &
    ppN3n, ppN4n, ppN5s, ppN6r, ppR1c, ppR1n, ppR1p, ppPhytoPlankton, BoxNumberZ, &
    NO_BOXES_Z, NO_BOXES_XY, BoxNumber, &
    BoxNumberXY, Depth, jPIY3c, PIc, jRIY3c, jRIY3n, jRIY3p, jRIY3s, jbotO2o, &
    jbotN1p, jbotN3n, jbotN4n, jbotN5s, jbotN6r, jbotR6c, jbotR6n, jbotR6p, jbotR6s, &
    jbotR1c, jbotR1n, jbotR1p, PELBOTTOM, &
    iiPhytoPlankton, iiP1, iiC, iiN, iiP, iiL, iiS, iiBen, iiPel, flux
  use mem_Param,  ONLY: AssignPelBenFluxesInBFMFlag, p_small


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
  real(RLEN), dimension(:), pointer  ::lcl_ppPhytoPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer  :: j
  real(RLEN)  :: uptake
  real(RLEN)  :: Pc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  BoxNumberZ = NO_BOXES_Z
    DO BoxNumberXY=1,NO_BOXES_XY
      BoxNumber=D3toD1(BoxNumberXY,BoxNumberZ)


      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate Phyto Fluxes to Filterfeeder from Pelagic for
      ! all phyt types/constituents
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        if ( jPIY3c(BoxNumberXY)> 0) then
          do i = 1 , ( iiPhytoPlankton)

          lcl_PhytoPlankton => PhytoPlankton(i,iiC)
          Pc  =   lcl_PhytoPlankton(BoxNumber)
          if ( Pc> p_small) then
            j = ppPhytoPlankton(i,iiC)
            uptake  =   jPIY3c(BoxNumberXY)* Pc/ PIc(BoxNumberXY)
            PELBOTTOM(j,BoxNumberXY)  =  PELBOTTOM(j,BoxNumberXY)  -uptake
            j = ppPhytoPlankton(i,iiN)
            lcl_PhytoPlankton => PhytoPlankton(i,iiN)
            PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                                      -uptake* lcl_PhytoPlankton(BoxNumber)/ Pc
            j = ppPhytoPlankton(i,iiP)
            lcl_PhytoPlankton => PhytoPlankton(i,iiP)
            PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                                      -uptake* lcl_PhytoPlankton(BoxNumber)/ Pc
            j = ppPhytoPlankton(i,iiL)
            lcl_PhytoPlankton => PhytoPlankton(i,iiL)
            PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                                      -uptake* lcl_PhytoPlankton(BoxNumber)/ Pc
            j = ppPhytoPlankton(i,iiS)
            if ( j> 0) then
              lcl_PhytoPlankton => PhytoPlankton(i,iiS)
              PELBOTTOM(j,BoxNumberXY) =  PELBOTTOM(j,BoxNumberXY)  &
                              -uptake* lcl_PhytoPlankton(BoxNumber)/ Pc
            end if

          end if

        end do

      end if


      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Net detritus Fluxes to Benthic from Pelagic by Y3
      ! net flux= uptake - excretion of food : flux may ber negative!
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux(BoxNumber, iiPel, ppR6c, ppR6c, -( jRIY3c(BoxNumberXY)/ &
        Depth(BoxNumber)) )
      call flux(BoxNumber, iiPel, ppR6n, ppR6n, -( jRIY3n(BoxNumberXY)/ &
        Depth(BoxNumber)) )
      call flux(BoxNumber, iiPel, ppR6p, ppR6p, -( jRIY3p(BoxNumberXY)/ &
        Depth(BoxNumber)) )
      call flux(BoxNumber, iiPel, ppR6s, ppR6s, -( jRIY3s(BoxNumberXY)/ &
        Depth(BoxNumber)) )

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! All Fluxes to Benthic from Pelagic defined for the
      ! Pelagic State variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( AssignPelBenFluxesInBFMFlag) then

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Nutrient Fluxes to Benthic from Pelagic
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call flux(BoxNumber, iiPel, ppO2o, ppO2o, jbotO2o(BoxNumberXY)/ &
          Depth(BoxNumber) )
        call flux(BoxNumber, iiPel, ppN1p, ppN1p, jbotN1p(BoxNumberXY)/ &
          Depth(BoxNumber) )
        call flux(BoxNumber, iiPel, ppN3n, ppN3n, jbotN3n(BoxNumberXY)/ &
          Depth(BoxNumber) )
        call flux(BoxNumber, iiPel, ppN4n, ppN4n, jbotN4n(BoxNumberXY)/ &
          Depth(BoxNumber) )
        call flux(BoxNumber, iiPel, ppN5s, ppN5s, jbotN5s(BoxNumberXY)/ &
          Depth(BoxNumber) )
        call flux(BoxNumber, iiPel, ppN6r, ppN6r, jbotN6r(BoxNumberXY)/ &
          Depth(BoxNumber) )

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! PhytoPlankton Fluxes to Benthic from Pelagic by Y3
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        do i = 1 , ( iiPhytoPlankton)

              j = ppPhytoPlankton(i,iiC)
              call flux(BoxNumber, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(BoxNumber) )
              j = ppPhytoPlankton(i,iiN)
              call flux(BoxNumber, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(BoxNumber) )
              j = ppPhytoPlankton(i,iiP)
              call flux(BoxNumber, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(BoxNumber) )
              j = ppPhytoPlankton(i,iiL)
              call flux(BoxNumber, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                Depth(BoxNumber) )
              if ( i== iiP1) then
                !No Y3.s defined, all silicate uptaken is moved into sink
                j = ppPhytoPlankton(i,iiS)
                call flux(BoxNumber, iiPel, j, j, PELBOTTOM(j,BoxNumberXY)/ &
                  Depth(BoxNumber) )
              end if
        end do

        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Total Sedimentation flux to Benthic from Pelagic defined
        ! for the pealgic state variables.
        ! (See sedimentation for definition of the fluxes of benthic &
        ! variables)
        ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call flux(BoxNumber, iiPel, ppR6c, ppR6c, ( jbotR6c(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR6n, ppR6n, ( jbotR6n(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR6p, ppR6p, ( jbotR6p(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR6s, ppR6s, ( jbotR6s(BoxNumberXY)/ &
          Depth(BoxNumber)) )

        call flux(BoxNumber, iiPel, ppR1c, ppR1c, ( jbotR1c(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR1n, ppR1n, ( jbotR1n(BoxNumberXY)/ &
          Depth(BoxNumber)) )
        call flux(BoxNumber, iiPel, ppR1p, ppR1p, ( jbotR1p(BoxNumberXY)/ &
          Depth(BoxNumber)) )
      end if


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
