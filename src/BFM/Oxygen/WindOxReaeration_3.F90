#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: WindOxReaeration_3
!
! DESCRIPTION
!   Model describes reaeration between air and water column.
!       as forced by temperature and wind.
!
!       The equation and correlation used in this routine
!       are found in the 
!		R. Wanninkhof (1992), Relationship between windspeed and gas
!		exchange over the oecean
!               J. GeoPhys. Res. 97, 7373-7382
!
! !INTERFACE
  subroutine OxygenReaerationDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: O2o, D3STATE, EICE
  use mem, ONLY: ppO2o, NO_BOXES_XY, &
    BoxNumberXY, EWIND, ETW, cxoO2, Depth, &
    jsurO2o, iiBen, iiPel, flux
#endif
  use mem_Param,  ONLY: AssignAirPelFluxesInBFMFlag
  use mem_WindOxReaeration_3
#ifdef BFM_GOTM
  use bio_var, ONLY: SRFindices
#else
  use api_bfm, ONLY: SRFindices
#endif
  use global_interface,   ONLY: CalcSchmidtNumberOx

!
! !AUTHORS
!   11 March 1998 Original version by P. Ruardij
!	              JWB 1999/03/25 Corrected k 
!
! !REVISION_HISTORY
!   
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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
  real(RLEN)  :: reacon
  real(RLEN)  :: p_schmidt
  integer     :: ksur

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY
     ksur = SRFindices(BoxNumberXY)

      !
      ! Reference Schimdt number for CO2 (=reference) of 660.0
      !
      p_schmidt  =   CalcSchmidtNumberOx(  ETW(ksur))/ 660.0D+00

      !
      ! Calculate wind dependency:
      !`
      reacon  =   k* (EWIND(BoxNumberXY))**(2.0D+00)/ sqrt(  p_schmidt)

      jsurO2o(BoxNumberXY)  = jsurO2o(BoxNumberXY) +  &
                              (ONE-EICE(BoxNumberXY))*reacon*( cxoO2(ksur)- O2o(ksur))

      if ( AssignAirPelFluxesInBFMFlag) then
        call flux(ksur, iiPel, ppO2o, ppO2o, jsurO2o(BoxNumberXY)/ &
          Depth(ksur) )
      end if

  end do

  end subroutine OxygenReaerationDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
