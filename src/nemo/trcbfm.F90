#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:
!
! !INTERFACE:
   subroutine trc_bfm( kt )
!
! !DESCRIPTION:
!  BFM stepping inside NEMO
!  This routine computes the BFM trends, integrates the benthic 
!  model with a simple Euler Forward and the pelagic model if 
!  a "no transport" simulation is prescribed.
!
! !USES:
   use global_mem, only: RLEN,LOGUNIT
   use constants,  only: SEC_PER_DAY
   use mem, only: D3STATE,D3SOURCE,NO_D3_BOX_STATES, &
                  D3STATETYPE,NO_BOXES
#ifdef EXPLICIT_SINK
   use mem, only: D3SINK
#endif 
#ifdef INCLUDE_BEN
   use mem, only: D2STATE_BEN,D2SOURCE_BEN,NO_D2_BOX_STATES_BEN, &
                  D2STATETYPE_BEN
#ifdef EXPLICIT_SINK
   use mem, only: D2SINK_BEN
#endif 
#endif
   use mem_param, only: CalcTransportFlag, CalcBenthicFlag, &
                        CalcPelagicFlag
   use time,      only: bfmtime
   ! NEMO
   use oce_trc          ! ocean dynamics and active tracers variables
   use trc              ! ocean passive tracers variables
   USE prtctl          ! Print control                    (prt_ctl routine)
   
   implicit none
!
! !INPUT PARAMETERS:
!
   integer, intent( IN ) ::  kt  ! ocean time-step index
   
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): M. Vichi
!
! !LOCAL VARIABLES:
   integer               :: j
   integer,save          :: first=0
   real(RLEN)            :: delt
   real(RLEN), allocatable,save,dimension(:,:) :: tmp
   real(RLEN), allocatable,save, dimension(:)  :: tmp2
   integer               :: AllocStatus,DeAllocStatus



!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Biological timestep 
   !---------------------------------------------
   delt  = rdt*real(nn_dttrc,RLEN)
   
   !---------------------------------------------
   ! BFM internal time
   !--------------------------------------------- 
   bfmtime%stepnow  = bfmtime%stepnow + nn_dttrc
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
#ifndef EXPLICIT_SINK
            D3State(j,:) = D3STATE(j,:) + delt*D3SOURCE(j,:)
#else
            D3State(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         end if
      end do
   end if
#ifdef INCLUDE_BEN
   if (CalcBenthicFlag /= 0) then
      do j=1,NO_D2_BOX_STATES_BEN
         if (D2STATETYPE_BEN(j).ge.0) then
#ifndef EXPLICIT_SINK
            D2STATE_BEN(j,:) = D2STATE_BEN(j,:) + delt*D2SOURCE_BEN(j,:)
#else
            D2STATE_BEN(j,:) = D2STATE_BEN(j,:) + delt*sum(D2SOURCE_BEN(j,:,:)-D2SINK_BEN(j,:,:),1)
#endif
         end if
      end do
   end if
#endif

   return
   end subroutine trc_bfm

!EOC

