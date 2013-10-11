#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: D2m, D1m
  use mem, ONLY: K3n,K4n,K6r,ppD2m, ppD1m,    &
    NO_BOXES_XY,    BoxNumberXY_ben, &
    InitializeModel, shiftD1m, shiftD2m, KNO3, KNH4, iiBen, iiPel, flux
  use mem,ONLY: jbotN3n,jbotN4n,N3n_Ben,N4n_Ben,K14n,K24n,D6m,D7m
#endif
  use constants,  ONLY: EQUATION, STANDARD, GET, LABDA_1, ONE_PER_DAY,MASS,INTEGRAL,LABDA_2
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
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
  real(RLEN)  :: Dxm
  real(RLEN)  :: Dzm
  real(RLEN)  :: pmM3n
  real(RLEN)  :: ushiftD2m
  real(RLEN)  :: D2mNew
  real(RLEN)  :: M3n_D1m
  real(RLEN)  :: M3n_Dxm
  real(RLEN)  :: sK3G4n
  real(RLEN)  :: lambda
  real(RLEN)  :: dummy

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: PrintSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate concentration of nitrate in porewater in M3n:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M3n_D1m = CalculateFromSet( KNO3(BoxNumberXY_ben), EQUATION, &
        STANDARD, D1m(BoxNumberXY_ben), dummy)
      Dxm=D1m(BoxNumberXY_ben)+D2m(BoxNumberXY_ben) * 0.5
      M3n_Dxm= CalculateFromSet( KNO3(BoxNumberXY_ben), EQUATION, &
        STANDARD, Dxm, dummy)

      select case ( M3n_D1m<= 0.0D+00)

        case( .TRUE. )

          ! Do nothing give a warning
          ! Let if D1m increase let D2m increase too to avoid zero thickness of 
          ! denitrifiaction layer
          if ( InitializeModel ==0 ) then
            write(LOGUNIT,'(''D1m='',F12.3,'' D2m='',F12.3)') &
                                     D1m(BoxNumberXY_ben), D2m(BoxNumberXY_ben)
            write(LOGUNIT,'(''D6m='',F12.3,'' D7m='',F12.3)') &
                                     D6m(BoxNumberXY_ben), D7m(BoxNumberXY_ben)
            write(LOGUNIT,'(''K3n='',F12.3,'' K4n='',F12.3)') &
                                     K3n(BoxNumberXY_ben), K4n(BoxNumberXY_ben)
            write(LOGUNIT,'(''K14n='',F12.3,'' K24n='',F12.3)') &
                                     K14n(BoxNumberXY_ben), K24n(BoxNumberXY_ben)
            write(LOGUNIT,'(''fluxN3='',F12.3,'' fluxK4n='',F12.3)') &
                                     jbotN3n(BoxNumberXY_ben), jbotN4n(BoxNumberXY_ben)
            write(LOGUNIT,'(''N3n='',F12.3,'' N4n='',F12.3)') &
                                     N3n_Ben(BoxNumberXY_ben), N4n_Ben(BoxNumberXY_ben)
            write(LOGUNIT,'(''K6r='',F12.3,'' lambda K3n='',F12.3)') &
                                     K6r(BoxNumberXY_ben),GetInfoFromSet( KNO3(BoxNumberXY_ben), GET, LABDA_1, 21)  
            control  =   PrintSet(  KNH4(BoxNumberXY_ben),"concentration nitrate on D1m < 0")
            control  =   PrintSet(  KNO3(BoxNumberXY_ben),"concentration nitrate on D1m < 0")
         endif

          D2mnew =D2m(BoxNumberXY_ben)
          shiftD2m(BoxNumberXY_ben) = max( min( p_d_tot- &
            p_clD1D2m, D2mNew), D1m(BoxNumberXY_ben)+ p_clD1D2m)- D2m(BoxNumberXY_ben)

        case( .FALSE. )

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! According solution nitrate concentrate decreases
          ! exponentially. Calculate depth at where below the 
          ! denitrification /m2 is  below a (small) fixed number. 
          ! Use this new depth as the (uncorrected) new denitrification depth
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          sK3G4n= GetInfoFromSet( KNO3(BoxNumberXY_ben), GET, LABDA_2, 31,dummy,dummy) 
          lambda= GetInfoFromSet( KNO3(BoxNumberXY_ben), GET, LABDA_1, 31,dummy,dummy) 
          D2mnew= log(-p_jlK3G4n*lambda/(sK3G4n* M3n_Dxm))/(2.0*lambda)
          Dzm=p_d_tot-p_clD1D2m
          D2mNew = min( Dzm, ( D2mnew + Dxm))

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! 1. Calculate uncorrected shift of D2.m
          ! 2. limit shift incase of D2mnew moves in the direction of D1m
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          shiftD2m(BoxNumberXY_ben) = max( D2mNew,D1m(BoxNumberXY_ben)+ p_clD1D2m) &
                                                              - D2m(BoxNumberXY_ben)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Correct by damping the change of D2m in case large changes:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          if ( InitializeModel== 0) then
            shiftD2m(BoxNumberXY_ben) = shiftD2m(BoxNumberXY_ben)/ &
              ONE_PER_DAY* (D2m(BoxNumberXY_ben)/( &
              D2m(BoxNumberXY_ben)+ abs(shiftD2m(BoxNumberXY_ben))))**(p_xdamping)  &
               * Dzm/(Dzm+D2mnew)


             call flux(BoxNumberXY_ben, iiBen, ppD2m, ppD2m, shiftD2m(BoxNumberXY_ben) )

          end if


      end select

  end do

  end subroutine BenDenitriDepthDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
