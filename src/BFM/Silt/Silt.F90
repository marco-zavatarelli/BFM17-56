#include "DEBUG.h"
#include "INCLUDE.h"
#include"cppdefs.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Silt
!
! DESCRIPTION
!
!
! !INTERFACE
  subroutine SiltDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,PI
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,  ONLY: R9x, ppR9x,ERHO,EWIND
#endif
  use mem, ONLY: Depth,ETW, NO_BOXES, iiBen, iiPel, flux_vector,InitializeModel
  use mem_Param,  ONLY: p_small, p_poro
  use global_interface,   ONLY: eTq
  use mem_Silt, ONLY: siltmethod


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:MM_vector, eTq_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!  
!
! !AUTHORS
!   Original version by P. Ruardij and  J. van der Molen
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



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     REAL(RLEN)                        :: delta
     REAL(RLEN)                        :: psilt
     REAL(RLEN)                        :: tdepth
     REAL(RLEN),dimension(NO_BOXES)    :: new_R9x
     REAL(RLEN),dimension(NO_BOXES)    :: rate
     REAL(RLEN)                        :: U, Hs, Tz, tau, uw, rate_1, R9x_new
     REAL(RLEN), PARAMETER             :: g=9.81, fw=0.1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     real(RLEN), external  :: GetDelta


     delta=GetDelta()
     psilt=(p_poro(1) - 0.38662 )/ 0.00415
     tdepth=sum(Depth)
     select case (siltmethod)
       case (1)
        STDERR , 'silt=', R9x(1:2)
        delta=GetDelta()
        psilt=max(0.01,(p_poro(1) - 0.38662 )/ 0.00415)
        new_R9x(:) = 10000.0 * min(1.0,10.0/tdepth)* max(0.5,psilt)/7.0 &
                      /eTq(  ETW(1), 2.0D+00) 
        if ( new_R9x(1) < 0.0 ) then
         STDERR, 'NewSilt=', new_R9x(1:2)
        endif
        select case ( InitializeModel ) 
         case(0) 
           rate=(new_R9x(:)-R9x(:))/delta
          case(1)
           rate=new_R9x(:)/delta;
        end select
       case (2)
         !JM simple wave generation: vRijn p 331
         U=0.7*(EWIND(1)**1.2)
         Hs=max(1.0,0.243*(U**2)/g)
         Hs=min(Hs,0.4*tdepth)
         Tz=max(5.0,8.14*U/g)
         uw=PI*Hs/(Tz*sinh(2*PI/Tz*sqrt(tdepth/g)))  !JM linear wave theory
         tau=0.5*fw*ERHO(1)*(uw**2)                     !JM 0.5*fw*rho*u**2
!        R9x_new = 0.5 * max(0.5,psilt) * tau
         R9x_new = 100.0 * max(0.3,psilt/100.0) * tau
         if ( InitializeModel == 0) then
           rate_1=(R9x_new-R9x(1))/delta
!          rate_1=max(rate_1,-R9x(1)*0.003)            !JM time lag in settling
           rate_1=max(rate_1,-R9x(1)*0.75)            !JM time lag in settling
         else
           rate_1=R9x_new/delta
         endif
         rate(:)=rate_1
     end select

      select case (InitializeModel)         !JM always 0 because InitBenthicNutrient happens first
        case (0)
          call flux_vector(iiPel,ppR9x,ppR9x,rate);
        case (1)
         R9x(:)=rate*delta
      end select

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
