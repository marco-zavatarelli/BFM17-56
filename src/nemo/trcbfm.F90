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
!  This routine computes the BFM trends
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s):
!  01Jan2000: Ver. ?.?.? (
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Compute external forcing functions
   !---------------------------------------------
   call envforcing_bfm

   !---------------------------------------------
   ! Compute Biogeochemical trends
   !---------------------------------------------
   call EcologyDynamics

   ! MAV: add a 1st order ODE solver for benthic in the future

   return
   end subroutine trcbfm

!EOC

