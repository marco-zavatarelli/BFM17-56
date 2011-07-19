#include "DEBUG.h"
#include "INCLUDE.h"
#ifdef INCLUDE_BEN
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CorrectConcNearBed
!
! DESCRIPTION
!
! !INTERFACE
      subroutine CorrectConcNearBed(depthlayer, sedi, fto, & 
                                    p_max, wf, correction)
! !USES:
     use constants, only: RLEN, SEC_PER_DAY
#ifdef BFM_GOTM
     use turbulence,  ONLY: kappa
#endif
#ifdef NOPOINTERS
     use mem
#else
     use mem,         ONLY: ETAUB,NO_BOXES_XY
#endif

     IMPLICIT NONE
! !INPUT PARAMETERS:
     real(RLEN), intent(IN),dimension(NO_BOXES_XY)   ::DepthLayer
     real(RLEN), intent(IN),dimension(NO_BOXES_XY)   ::Sedi
     real(RLEN), intent(IN)                          ::fto
     real(RLEN), intent(IN)                          ::p_max
     real(RLEN), intent(IN),dimension(NO_BOXES_XY)   ::wf          ! volumefiltered*Y3c (m/d)
! !OUTPUT PARAMETERS:
     real(RLEN),dimension(NO_BOXES_XY),intent(OUT)   ::correction

!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
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

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Local Variables
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      real(RLEN),dimension(NO_BOXES_XY)                ::b
      real(RLEN),dimension(NO_BOXES_XY)                ::f
      real(RLEN),dimension(NO_BOXES_XY)                ::r
      real(RLEN),parameter                             ::p_small=1.0D-10
#ifndef BFM_GOTM
      real(RLEN),parameter                             ::kappa=1.0D-5
#endif
    
      f=  min(fto,DepthLayer)
      b = min(100.0D+00,Sedi/(p_small+SEC_PER_DAY*kappa*ETAUB(:))) 

      where (b < 1.0D+00)
        r=1.0/(1.0D+00-b)*(0.5*DepthLayer)**b
        correction=r*(f**(1.0D+00-b) -p_small**(1.0D+00-b))/(f-p_small)
      elsewhere
        correction=1.0/p_small
      endwhere


      return
      end subroutine CorrectConcNearBed
#endif
!EOC
!-----------------------------------------------------------------------

