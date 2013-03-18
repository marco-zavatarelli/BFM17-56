#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Euler-forward time-integration
!
! !INTERFACE
   SUBROUTINE integrationEfw
!
! !DESCRIPTION
!  Euler-forward integration with time step adjustment
! !USES
   use global_mem, ONLY:RLEN
   use mem,  ONLY: NO_D3_BOX_STATES,NO_BOXES,D3SOURCE,D3STATE, &
                   D3STATETYPE
#ifndef ONESOURCE
   use mem,  ONLY: D3SINK
#endif
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
   use mem,  ONLY: NO_D2_BOX_STATES,D2SOURCE,D2STATE,NO_BOXES_XY, &
                   D2STATETYPE
#ifndef ONESOURCE
   use mem,  ONLY: D2SINK
#endif
#endif
   use standalone
   use api_bfm
   use time, ONLY: update_time
   implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO)
!
!
! COPYING
!
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
! !LOCAL VARIABLES:
   real(RLEN),parameter      :: eps=0.
   real(RLEN)                :: min3D,min2D
   integer                   :: i,j,ll
   integer,dimension(2,2)    :: blccc
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'integration efw: starting delt = ',delt
#endif
   bbccc3D=D3STATE
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
   bbccc2D=D2STATE
#endif
   TLOOP : DO
   ! Integration step:
      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
#ifdef D1SOURCE
            D3STATE(j,:) = D3STATE(j,:) + delt * D3SOURCE(j,:)
#else
            D3STATE(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:),1)
#endif
#else
            D3STATE(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         END IF
      END DO
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      DO j=1,NO_D2_BOX_STATES
         IF (D2STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
            D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:),1)
#else
            D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
#endif
         END IF
      END DO
#endif
      nmin=nmin+nstep 
   !  Check for negative concentrations
      min3D=minval(D3STATE)
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      min2D=minval(D2STATE)
      IF(min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
#else
      IF(min3D.lt.eps) THEN ! cut timestep
#endif
         IF (nstep.eq.1) THEN
            LEVEL1 'Necessary Time Step too small! Exiting...'
            blccc(:,1)=minloc(D3STATE)
#ifdef ONESOURCE
#ifdef D1SOURCE
            bbccc3D = D3SOURCE(:,:)
#else
            bbccc3D = sum(D3SOURCE(:,:,:),2)
#endif
#else
            bbccc3D = sum(D3SOURCE(:,:,:)-D3SINK(:,:,:),2)
#endif
            LEVEL1 'Pelagic Variable:',trim(var_names(stPelStateS+blccc(1,1)-1))
            LEVEL1 'Value: ',D3STATE(blccc(1,1),blccc(2,1)),' Rate: ', &
                        bbccc3D(blccc(1,1),blccc(2,1))
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
            blccc(:,2)=minloc(D2STATE)
            bbccc2D = sum(D2SOURCE(:,:,:)-D2SINK(:,:,:),2)
            LEVEL1 'Benthic Variable:',trim(var_names(stBenStateS+blccc(1,2)-1))
            LEVEL1 'Value: ',D2STATE(blccc(1,2),blccc(2,2)),' Rate: ', &
                           bbccc2D(blccc(1,2),blccc(2,2))
            LEVEL1 'EXIT at  time ',timesec
#endif
            STOP 'integration-efw'
         END IF
         nstep=nstep/2
         nmin=0
         D3STATE=bbccc3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
         D2STATE=bbccc2D
#endif
         dtm1=delt
         delt=nstep*mindelt
         timesec=ntime*maxdelt
         LEVEL1 'Time Step cut! delt= ',delt,' nstep= ',nstep
      ELSE
#ifdef DEBUG
      LEVEL2 'Internal time step= ',nmin
#endif
         IF(nmin.eq.nmaxdelt) EXIT TLOOP
         timesec=timesec+delt
#ifdef DEBUG
      LEVEL2 'Internal time= ',timesec
#endif
      END IF
      ! Recalculate Sources:
      call ResetFluxes
      call envforcing_bfm
      call EcologyDynamics
   END DO TLOOP
   nstep=nmaxdelt
   nmin=0
   ntime=ntime+1
   call update_time(ntime)
   dtm1=delt
   delt=nstep*mindelt
   timesec=delt*ntime
#ifdef DEBUG
   LEVEL1 'ntime: ',ntime
   LEVEL1 'Integration time: ',timesec
#endif

   END SUBROUTINE integrationEfw
!EOC
!-----------------------------------------------------------------------
