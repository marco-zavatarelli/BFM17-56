#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:
!
! !INTERFACE:
   subroutine trcbfm
!
! !DESCRIPTION:
!  BFM stepping inside NEMO
!  This routine computes the BFM trends, integrates the benthic 
!  model with a simple Euler Forward and the pelagic model if 
!  a "no transport" simulation is prescribed.
!
! !USES:
   use global_mem, only: RLEN
   use constants,  only: SEC_PER_DAY
   use mem, only: D3STATE,D3SOURCE,NO_D3_BOX_STATES, &
                  D3STATETYPE,NO_BOXES
#ifndef ONESOURCE
   use mem, only: D3SINK
#endif 
#ifdef INCLUDE_BEN
   use mem, only: D2STATE,D2SOURCE,NO_D2_BOX_STATES, &
                  D2STATETYPE
#ifndef ONESOURCE
   use mem, only: D2SINK
#endif 
#endif
   use mem_param, only: CalcTransportFlag, CalcBenthicFlag, &
                        CalcPelagicFlag
   ! NEMO
   use oce_trc          ! ocean dynamics and active tracers variables
   use trc              ! ocean passive tracers variables

   implicit none
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): M. Vichi
!
! !LOCAL VARIABLES:
   integer               :: j
   real(RLEN)            :: delt
   real(RLEN)            :: tmp(NO_D3_BOX_STATES,NO_BOXES)
   real(RLEN)            :: tmp2(NO_BOXES)
!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Biological timestep 
   !---------------------------------------------
   delt  = rdt*FLOAT(ndttrc)

   !---------------------------------------------
   ! Compute external forcing functions
   !---------------------------------------------
   call envforcing_bfm

   !---------------------------------------------
   ! Compute Biogeochemical trends
   !---------------------------------------------
   call EcologyDynamics

   !---------------------------------------------
   ! Basic ODE solver for benthic variables or for 
   ! 0D test simulations when coupled with NEMO
   !---------------------------------------------
   if (CalcPelagicFlag .AND. .NOT.CalcTransportFlag) then
      do j=1,NO_D3_BOX_STATES
         if (D3STATETYPE(j).ge.0) then
#ifdef ONESOURCE
            D3State(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:),1)
#else
            D3State(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         end if
      end do
   end if
#ifdef INCLUDE_BEN
   if (CalcBenthicFlag /= 0) then
      do j=1,NO_D2_BOX_STATES
         if (D2STATETYPE(j).ge.0) then
#ifdef ONESOURCE
            D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:),1)
#else
            D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
#endif
         end if
      end do
   end if
#endif

   return
   end subroutine trcbfm

!EOC

