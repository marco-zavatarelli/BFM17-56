#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenDenitriDepth
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenDenitriDepthDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D2m
  ! The following Benthic-states are used (NOT in fluxes): K3n,K4n,D1m
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, BoxNumberY, NO_BOXES_Y, &
  ! BoxNumber, BoxNumberXY, dummy, InitializeModel
  ! The following Benthic 1-d global boxvars are modified : shiftD1m,shiftD2m
  ! The following Benthic 1-d global boxvars  are used: KNO3, KNH4
  ! The following 0-d global parameters are used: p_d_tot, p_clD1D2m
  ! The following global constants are used: RLEN
  ! The following constants are used: EQUATION, STANDARD, GET, LABDA_1, &
  ! ONE_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: D2m, D1m, D2STATE
  use mem, ONLY: K3n,K4n,ppD2m, ppD1m, NO_BOXES_XY, BoxNumberXY, dummy, &
    InitializeModel, shiftD1m, shiftD2m, KNO3, KNH4, iiBen, iiPel, flux
  use mem,ONLY: jbotN3n,jbotN4n,N3n_Ben,N4n_Ben,K14n,K24n,D6m,D7m
  use constants,  ONLY: EQUATION, STANDARD, GET, LABDA_1, ONE_PER_DAY
  use mem_Param,  ONLY: p_d_tot, p_clD1D2m
  use mem_BenDenitriDepth


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:CalculateFromSet, GetInfoFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface,   ONLY: CalculateFromSet, GetInfoFromSet
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
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
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: control
  real(RLEN)  :: pmM3n
  real(RLEN)  :: ushiftD2m
  real(RLEN)  :: D2mNew
  real(RLEN)  :: M3n_D1m

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: PrintSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    do BoxNumberXY=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate concentration of nitrate in porewater in M3n:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M3n_D1m = CalculateFromSet( KNO3(BoxNumberXY), EQUATION, &
        STANDARD, D1m(BoxNumberXY), dummy)

      select case ( M3n_D1m<= 0.0D+00)

        case( .TRUE. )

          ! Do nothing give a warning
          ! Let if D1m increase let D2m increase too to avoid zero thickness of 
          ! denitrifiaction layer
          if ( InitializeModel ==0 ) then
            write(LOGUNIT,'(''D1m='',F12.3,'' D2m='',F12.3)') D1m(BoxNumberXY), D2m(BoxNumberXY)
            write(LOGUNIT,'(''D6m='',F12.3,'' D7m='',F12.3)') D6m(BoxNumberXY), D7m(BoxNumberXY)
            write(LOGUNIT,'(''K3n='',F12.3,'' K4n='',F12.3)') K3n(BoxNumberXY), K4n(BoxNumberXY)
            write(LOGUNIT,'(''K14n='',F12.3,'' K24n='',F12.3)') K14n(BoxNumberXY), K24n(BoxNumberXY)
            write(LOGUNIT,'(''fluxN3='',F12.3,'' fluxK4n='',F12.3)') jbotN3n(BoxNumberXY), jbotN4n(BoxNumberXY)
            write(LOGUNIT,'(''N3n='',F12.3,'' N4n='',F12.3)') N3n_Ben(BoxNumberXY), N4n_Ben(BoxNumberXY)
            control  =   PrintSet(  KNH4(BoxNumberXY),"concentration nitrate on D1m < 0")
            control  =   PrintSet(  KNO3(BoxNumberXY),"concentration nitrate on D1m < 0")
         endif

          D2mnew =D2m(BoxnumberXY)+shiftD1m(boxNumberXY)
          shiftD2m(BoxNumberXY) = max( min( p_d_tot- &
            p_clD1D2m, D2mNew), D1m(BoxNumberXY)+ p_clD1D2m)- D2m(BoxNumberXY)

        case( .FALSE. )

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Calculate fraction with which new depth of D2.m is calculated:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          pmM3n  =   max(  p_pmM3n* M3n_D1m,  p_clM3n_D2)/ (1.0D-80 +M3n_D1m)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! According solution nitrate concentrate decreases
          ! exponentially. Use exponent to calculate depth at which
          ! concentration is pmM3n * M3n_D1m. Use this new depth as
          ! uncorrected new denitrification depth
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          D2mNew = min( p_d_tot-p_clD1D2m, ( log( pmM3n)/ &
            GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_1, 21) &
                                                     + D1m(BoxNumberXY)))

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! 1. Calculate uncorrected shift of D2.m
          ! 2. limit shift incase of D2mnew moves in the direction of D1m
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          shiftD2m(BoxNumberXY) = max( min( p_d_tot- &
            p_clD1D2m, D2mNew), D1m(BoxNumberXY)+ p_clD1D2m)- D2m(BoxNumberXY)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Correct by damping the change of D2m in case large changes:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          if ( InitializeModel== 0) then
            shiftD2m(BoxNumberXY) = shiftD2m(BoxNumberXY)/ &
              ONE_PER_DAY* (D2m(BoxNumberXY)/( &
              D2m(BoxNumberXY)+ abs(shiftD2m(BoxNumberXY))))**(p_xdamping)


             call flux(BoxNumberXY, iiBen, ppD2m, ppD2m, shiftD2m(BoxNumberXY) )

          end if


      end select


  end DO

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
