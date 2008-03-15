#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenGlobal
!
! DESCRIPTION
!  Initialise global variables in the benthic system 
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenGlobalDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D6m, D7m, D8m, D9m
  ! The following Benthic-states are used (NOT in fluxes): Y2c, Y5c, Y1c, Y4c
  ! The following Benthic 1-d global boxvars are modified : turenh
  ! The following Benthic 1-d global boxvars got a value: irrenh, rrBTo, &
  ! reBTn, reBTp, rrATo, reATn, reATp
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem,  ONLY: D6m, D7m, D8m, D9m, Y2c, Y5c, Y1c, Y4c
#endif
  use mem, ONLY: ppD6m, ppD7m, ppD8m, ppD9m, ppY2c, ppY5c, ppY1c, &
    ppY4c, turenh, irrenh, rrBTo, rrATo, reBTn, reBTp, reATn, reATp, &
    jbotO2o,jbotN1p,jbotN3n,jbotN4n,jbotN5s, jbotN6r, &
    ETW_Ben, NO_BOXES_XY, iiBen, iiPel, flux_vector
#ifdef INCLUDE_BENCO2
  use mem, ONLY:  jbotO3h,jbotO3c
#endif


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: MM_vector

!  
!
! !AUTHORS
!   P.Ruardij
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! In the betnhic processes all respirations and excretions
  ! are added to the rr???? and re??? variables.
  ! There rates are input to the Benthic Nutrient model
  ! first these variables are initialized:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rrBTo(:)  =   ZERO  ! mgO2/m2 # Total Benthic oxic respiration
  reBTn(:)  =   ZERO  ! mmN/m2  # Total Benthic oxic N mineralization
  reBTp(:)  =   ZERO  ! mmP/m2  # Total Benthic oxic P mineralization
  rrATo(:)  =   ZERO  ! mgO2/m2 # Total Benthic anoxic respiration
  reATn(:)  =   ZERO  ! mmN/m2  # Total Benthic anoxic N mineralization
  reATp(:)  =   ZERO  ! mmP/m2  # Total Benthic anoxic P mineralization

 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! All nutrient flux rates are initalized here because nutrient are given
  ! back by Filterfeeders and by the nutrient regeration model 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   
#ifdef INCLUDE_BENCO2
  jbotO3h(:)=ZERO
  jbotO3c(:)=ZERO
#endif
  jbotO2o(:)=ZERO
  jbotN1p(:)=ZERO
  jbotN3n(:)=ZERO
  jbotN4n(:)=ZERO
  jbotN5s(:)=ZERO
  jbotN6r(:)=ZERO

  end subroutine BenGlobalDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
