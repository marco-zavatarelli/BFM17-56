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
   use mem, only: D3STATE,D3SOURCE,D3SINK,NO_D3_BOX_STATES, &
                  D2STATE,D2SOURCE,D2SINK,NO_D2_BOX_STATES, &
                  D3STATETYPE,D2STATETYPE
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
!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Biological timestep (in days)
   !---------------------------------------------
   delt  = rdt*FLOAT(ndttrc)/SEC_PER_DAY

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
         if (D3STATETYPE(j).ge.0) &
            D3STATE(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
      end do
   end if
   if (CalcBenthicFlag /= 0) then
      do j=1,NO_D2_BOX_STATES
         if (D2STATETYPE(j).ge.0) &
            D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
      end do
   end if

   return
   end subroutine trcbfm

!EOC

