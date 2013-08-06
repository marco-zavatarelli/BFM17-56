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
     use mem_Param, only: p_small
     use global_mem, only: ONE,ZERO
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
     real(RLEN), intent(IN)                              ::fto
     real(RLEN), intent(IN)                              ::p_max
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
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
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
#ifndef BFM_GOTM
      real(RLEN),parameter                             ::kappa=1.0E-5_RLEN
#endif
    
      f=  min(fto,DepthLayer)
      b = min(100.0_RLEN,Sedi/(p_small+SEC_PER_DAY*kappa*ETAUB(:))) 

      where (b < ONE)
        r=ONE/(ONE-b)*(0.5*DepthLayer)**b
        correction=r*(f**(ONE-b) -p_small**(ONE-b))/(f-p_small)
      elsewhere
        correction=ONE/p_small
      endwhere


      return
      end subroutine CorrectConcNearBed

#endif
!EOC
!-----------------------------------------------------------------------

