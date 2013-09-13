#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ControlBennutBuffers
!
! DESCRIPTION
!   This process forces the system to calculate a flux from
!       dissolved pelagic nutrient to the particualte nutreit pool
!       in the sediment. The reason for this is that if one of the
!       benthic pools become zero the partial differential equations 
!       used to calclate the diagenict processes and vertical fluxes
!       are not valid any more.
!       Therfor this process maximize input to benthos in
!       cases when nearly no particulate detritus enters the benthic system.
!       This is done by transferring dissolved nutrients N3.n,N4.n,N5.s
!       in the layer adjacent to the sediment into
!       benthic particulate nutrients (Q6.n,Q6.p,Q6.s).
!       A fixed fraction per day will be transferred to the benthos :
!        a. Minimized to  a fixed maximum value
!        b. Only when the settling rate is lower than the fixed maximum value
!        c. Only when het particulate nutrient buffer is really low.
!           This is determined to check if the quotient of the nutrient buffer
!           and the  fixed fraction is below a parameter, standing for the
!           specific refill rate.
!       The values of fixed fraction and the maximum value for the 
!       different nutrients are calculated according to Redfield
!       The user decides which nutrients should be handled
!          
! !INTERFACE
  subroutine ControlBennutBuffersDynamics
!
#ifdef INCLUDE_BEN
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO
  use mem_Param, ONLY: p_small
  use mem,  ONLY: Q6n, D9m, D7m, Q6p, D8m, Q6s, D2STATE_BEN
  use mem, ONLY: ppQ6n, ppD9m, ppD7m, ppQ6p, ppD8m, ppQ6s, &
    NO_BOXES_XY, BoxNumberXY_ben, jbotN4n, jbotN1p, jbotN5s, N4n_Ben, Depth_Ben, &
    jbotR6n, N1p_Ben, jbotR6p, N5s_Ben, jbotR6s, iiBen, iiPel, flux
  use mem_ControlBennutBuffers


!
! !AUTHORS
!   Piet Ruardij
!
!
! !REVISION_HISTORY
!   Created at Fri Jan 27 08:47:37 CET 2006
!
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij & M. Vichi
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
  real(RLEN)  :: r
  real(RLEN)  :: rate

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    do BoxNumberXY_ben=1,NO_BOXES_XY

      !----N:----------------------------------------------------------------
      if ( p_control_n) then
        rate = min( p_rn, p_s* N4n_Ben(BoxNumberXY_ben)* &
          Depth_Ben(BoxNumberXY_ben))
        if ( -jbotR6n(BoxNumberXY_ben)< p_rn .AND. Q6n(BoxNumberXY_ben)/ rate< &
          p_refill_time) then
          rate  =   max(  ZERO,  rate+ jbotR6n(BoxNumberXY_ben))
          call flux(BoxNumberXY_ben, iiBen, ppQ6n, ppQ6n, -(- rate) )
          call flux(BoxNumberXY_ben, iiBen, ppD9m, ppD9m, ( ZERO- &
            D7m(BoxNumberXY_ben))* rate/( p_small+ Q6n(BoxNumberXY_ben)) )
          jbotN4n(BoxNumberXY_ben)  =   jbotN4n(BoxNumberXY_ben)- rate
        end if

      end if

      !----P:----------------------------------------------------------------
      if ( p_control_p) then
        r  =   p_rn/ 16.0_RLEN  ! Correct according Redfield.
        rate  =   min(  r,  p_s* N1p_Ben(BoxNumberXY_ben)* Depth_Ben(BoxNumberXY_ben))
        if ( -jbotR6p(BoxNumberXY_ben)< r .AND. Q6p(BoxNumberXY_ben)/ rate< p_refill_time) &
          then
          rate  =   max(  ZERO,  rate+ jbotR6p(BoxNumberXY_ben))
          call flux(BoxNumberXY_ben, iiBen, ppQ6p, ppQ6p, -(- rate) )
          call flux(BoxNumberXY_ben, iiBen, ppD9m, ppD9m, ( ZERO- &
            D8m(BoxNumberXY_ben))* rate/( p_small+ Q6p(BoxNumberXY_ben)) )
          jbotN1p(BoxNumberXY_ben)  =   jbotN1p(BoxNumberXY_ben)- rate
        end if

      end if

      !----Si:---------------------------------------------------------------
      if ( p_control_s) then
        r  =   p_rn* 22.0_RLEN/ 16.0_RLEN  ! Correct according Readfield.
        rate  =   min(  r,  p_s* N5s_Ben(BoxNumberXY_ben)* Depth_Ben(BoxNumberXY_ben))
        if ( -jbotR6s(BoxNumberXY_ben)< r .AND. Q6s(BoxNumberXY_ben)/ rate< p_refill_time) &
          then
          rate  =   max(  ZERO,  rate+ jbotR6s(BoxNumberXY_ben))
          call flux(BoxNumberXY_ben, iiBen, ppQ6s, ppQ6s, -(- rate) )
          call flux(BoxNumberXY_ben, iiBen, ppD9m, ppD9m, ( ZERO- &
            D9m(BoxNumberXY_ben))* rate/( p_small+ Q6s(BoxNumberXY_ben)) )
          jbotN5s(BoxNumberXY_ben)  =   jbotN5s(BoxNumberXY_ben)- rate
        end if

      end if

  end do

#endif
  end subroutine ControlBennutBuffersDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
