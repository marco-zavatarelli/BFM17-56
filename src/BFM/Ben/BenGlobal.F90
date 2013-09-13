#include "DEBUG.h"
#include "INCLUDE.h"
#ifdef INCLUDE_BEN
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: D6m, D7m, D8m, D9m, Y2c, Y5c, Y1c, Y4c
  use mem, ONLY: ppD6m, ppD7m, ppD8m, ppD9m, ppY2c, ppY5c, ppY1c, &
    ppY4c, turenh, irrenh, rrBTo, rrATo, reBTn, reBTp, reATn, reATp, &
    jbotO2o,jbotN1p,jbotN3n,jbotN4n,jbotN5s, jbotN6r, &
    ETW_Ben, iiBen, iiPel, flux_vector
#ifdef INCLUDE_BENCO2
  use mem, ONLY:  jbotO3h,jbotO3c
#endif
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
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
#endif
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
