#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   RecalcPenetrationDepth.f90
!
! FILE
!   RecalcPenetrationDepth.f90
!
! DESCRIPTION
!   
!	Recompute the penetration depth of detritus to conserve mass
!       in the lower layers. 
!       The new thickness is calculated Newton's rule of approximation:
!           x(i+1)=x(i)-f(x)/f'(x for f(x)==0
!  
! !INTERFACE
    subroutine RecalcPenetrationDepth(D1m, Dxm, input, mass,newDxm )
!
! !AUTHORS
!   Piet Ruardij   
!
#ifdef INCLUDE_BEN
! !USES:
     use global_mem, ONLY:RLEN,ONE
     use mem,ONLY:NO_BOXES_XY
     IMPLICIT  NONE
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2007 P. Ruardij, the BFM team
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
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOC
        real(RLEN),intent(IN)     :: D1m(NO_BOXES_XY)
        real(RLEN),intent(IN)     :: Dxm(NO_BOXES_XY)
        real(RLEN),intent(IN)     :: input(NO_BOXES_XY)
        real(RLEN),intent(IN)     :: mass(NO_BOXES_XY)
        real(RLEN),intent(OUT)    :: newDxm(NO_BOXES_XY)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
        real(RLEN),parameter      :: EPS=1.0e-5_RLEN
        real(RLEN)                :: alpha
        real(RLEN)                :: old
        real(RLEN)                :: fx
        real(RLEN)                :: dfx
        real(RLEN)                :: c
        real(RLEN)                :: newalpha
        integer                   :: i
      
        do i = 1,NO_BOXES_XY ! loop over no. of grid points
           alpha=ONE/Dxm(i)
           ! calculate total mass present below D1m
           old = mass(i) * exp(-alpha *D1m(i))/alpha 
           c=ONE
           ! start iteration with old alpha
           newalpha = alpha
           do while ( c > EPS) 
              ! calc mass below D1m with new input
              fx = (mass(i) + input(i)) * exp(-newalpha *D1m(i))/newalpha
              ! determine first derivative
              dfx  = fx * ( -D1m(i) - ONE/newalpha ) 
              ! keep old value
              c =newalpha
              ! calc new alpha for (fx-old)==0
              newalpha = newalpha - (fx-old) /dfx
              ! use c to calculate iteration precision
              c = abs(c -newalpha)/c
           end do
           newDxm(i) = ONE/newalpha
        end do

        return
#endif
    end subroutine  RecalcPenetrationDepth
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
