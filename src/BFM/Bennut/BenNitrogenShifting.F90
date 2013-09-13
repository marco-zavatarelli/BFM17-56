#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenNitrogenShifting
!
! DESCRIPTION
!   Description of shifting of dissolved N (amm and nitrate) between
!       layers
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenNitrogenShiftingDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,ONE
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: K14n, K4n, K24n, K3n, D1m, D7m, D2m
  use mem, ONLY: ppK14n, ppK4n, ppK24n, ppK3n, ppD1m, ppG4n, ppD7m, &
    ppD2m,    NO_BOXES_XY,   &
     BoxNumberXY_ben, LocalDelta, shiftD1m, KNH4, reATn, shiftD2m, &
    KNO3, jK34K24n, jK13K3n, iiBen, iiPel, flux
#endif
  use constants,  ONLY: SHIFT, LAYER1, DERIVATIVE, RFLUX, LAYER2
  use mem_Param,  ONLY: p_poro, p_clDxm, p_d_tot
  use mem_BenthicNutrient3, ONLY:p_max_shift_change,p_max_state_change

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface,   ONLY: CalculateFromSet


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:insw, IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: insw, IntegralExp
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
!
! !REVISION_HISTORY
!   April 15, 1994 by EGM Embsen and P Ruardij:
!               Created a new version of the this process
!               so that it can be used with OpenSESAME.
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
  real(RLEN)  :: Dnew, dummy
  real(RLEN)  :: shiftmass
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: jK14K4n
  real(RLEN)  :: jK24K14n
  real(RLEN)  :: alpha
  real(RLEN)  :: r

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY_ben=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium Fluxes at the oxic/denitrification boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dnew  =   D1m(BoxNumberXY_ben)+ LocalDelta* shiftD1m(BoxNumberXY_ben)

      ! Calculate mass shifted in upwards direction:
      shiftmass = CalculateFromSet( KNH4(BoxNumberXY_ben), SHIFT, LAYER1, &
        D1m(BoxNumberXY_ben), Dnew)/ LocalDelta

      jK14K4n = CalculateFromSet( KNH4(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
        D1m(BoxNumberXY_ben), ZERO)+ shiftmass

      call LimitShift(jK14K4n,K4n(BoxNumberXY_ben),K14n(BoxNumberXY_ben),p_max_shift_change)
      call flux(BoxNumberXY_ben, iiBen, ppK14n, ppK4n,   jK14K4n* insw(  jK14K4n) )
      call flux(BoxNumberXY_ben, iiBen, ppK4n, ppK14n, - jK14K4n* insw( -jK14K4n) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! All the nutrient mineralization source term in the anoxic layer
      ! has been added to K14.n in BenBacDynamics
      ! However in the model this layer is subdivided and hence a partition
      ! flux is here calculated according to the exponential distribution.
      !
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D7.m is the average penetration depth for N-detritus
      !                 +
      ! Anoxic Mineralization at D1.m, using the exponential distribution
      !                 +
      !          Anoxic Mineralization at D2.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   ONE/ max(  p_clDxm,  D7m(BoxNumberXY_ben))
      zuD1 = max( 1.E-20_RLEN, reATn(BoxNumberXY_ben))/ p_poro(BoxNumberXY_ben)/ &
                        IntegralExp( -alpha, p_d_tot- D1m(BoxNumberXY_ben))
      zuD2  =   zuD1* exp( - alpha*( D2m(BoxNumberXY_ben)- D1m(BoxNumberXY_ben)))

      jK24K14n = - zuD2* p_poro(BoxNumberXY_ben)* IntegralExp( - alpha, &
        p_d_tot- D2m(BoxNumberXY_ben))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium Fluxes at the denitrification/anoxic boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       
      Dnew  =   D2m(BoxNumberXY_ben)+ LocalDelta* shiftD2m(BoxNumberXY_ben)

      ! Calculate mass shifted in upwards direction:
      shiftmass = CalculateFromSet( KNH4(BoxNumberXY_ben), SHIFT, LAYER2, &
        D2m(BoxNumberXY_ben), Dnew)/ LocalDelta  &
      + CalculateFromSet( KNH4(BoxNumberXY_ben), DERIVATIVE, &
        RFLUX, D2m(BoxNumberXY_ben), ZERO)

      call LimitShift(shiftmass,K14n(BoxNumberXY_ben),K24n(BoxNumberXY_ben),p_max_shift_change)
      jK24K14n=jK24K14n + shiftmass
       
    
      call flux(BoxNumberXY_ben, iiBen, ppK24n, ppK14n, jK24K14n* insw( jK24K14n) )
      call flux(BoxNumberXY_ben, iiBen, ppK14n, ppK24n,-jK24K14n* insw(-jK24K14n) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Fluxes at the lower boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Ammonium:
      jK34K24n(BoxNumberXY_ben)  = CalculateFromSet( KNH4(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
        p_d_tot, ZERO)
      call flux(BoxNumberXY_ben, iiBen, ppK24n, ppK24n, jK34K24n(BoxNumberXY_ben) )


      ! Nitrate:
      shiftmass = CalculateFromSet( KNO3(BoxNumberXY_ben), SHIFT, LAYER2, &
        D2m(BoxNumberXY_ben), Dnew)/ LocalDelta
      shiftmass=shiftmass * insw(shiftmass *shiftD2m(BoxNumberXY_ben))
      jK13K3n(BoxNumberXY_ben)  = CalculateFromSet( KNO3(BoxNumberXY_ben), DERIVATIVE, RFLUX, &
        D2m(BoxNumberXY_ben), dummy)+ shiftmass


      call LimitChange(1,jK13K3n(BoxNumberXY_ben),K3n(BoxNumberXY_ben),p_max_shift_change)
      call flux(BoxNumberXY_ben, iiBen, ppK3n, ppK3n,  jK13K3n(BoxNumberXY_ben) ) 


  end do

  end subroutine BenNitrogenShiftingDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
