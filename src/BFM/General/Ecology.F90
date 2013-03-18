#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Ecology
!
! DESCRIPTION
!   !	This submodel calls all other submodels
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine EcologyDynamics
!
! !USES:
  ! The following 0-d global parameters are used: CalcPelagicFlag, &
  ! CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem,  ONLY: iiBen, iiPel, iiReset, flux
  use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
  use mem_Param,  ONLY: CalcPelagicFlag, CalcBenthicFlag, CalcConservationFlag
#ifdef INCLUDE_SEAICE
  use mem_Param,  ONLY: CalcSeaiceFlag
#endif

!  
!
! !AUTHORS
!   Marcello Vichi & Piet Ruardij
!
!
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
#ifdef DEBUG
  call  flux(1,iiReset,1,1,0.0)
#endif

#ifdef INCLUDE_SEAICE
  if ( CalcSeaiceFlag) then

    call SeaiceSystemDynamics

  end if
#endif

  if ( CalcPelagicFlag) then

    call PelagicSystemDynamics

  end if
#ifdef INCLUDE_BEN
  if ( CalcBenthicFlag > 0 ) then

         call SettlingDynamics
  
       select case ( CalcBenthicFlag)

         case ( BENTHIC_RETURN )  ! Simple benthic return
           call BenthicReturn1Dynamics

         case ( BENTHIC_BIO )  ! Intermediate benthic return
           call PelForcingForBenDynamics
           call BenthicSystemDynamics
           call BenthicNutrient2Dynamics

         case ( BENTHIC_FULL )  ! Full benthic nutrients
           call PelForcingForBenDynamics
           call BenthicSystemDynamics
           call BenthicNutrient3Dynamics

       end select

       call ControlBennutBuffersDynamics

       call BentoPelCoupDynamics

       call SedimentationDynamics

  else 
       ! This case considers an inactive benthic system 
       ! (the benthic arrays are defined but not used)
       ! only the net sink at the bottom is computed
       call SettlingDynamics
       call BentoPelCoupDynamics

  endif
#else
  ! only the net sink at the bottom is computed
  call SettlingDynamics
  call BentoPelCoupDynamics

#endif

  if (CalcConservationFlag) &
     call CheckMassConservationDynamics

  end subroutine EcologyDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
