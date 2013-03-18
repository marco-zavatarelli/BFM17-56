#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Runge-Kutta 2nd-order time-integration
!
! !INTERFACE
   SUBROUTINE integrationRK2
!
! !DESCRIPTION
!  Runge-Kutta 2nd-order integration with time step adjustment
! !USES
   use global_mem, ONLY:RLEN,ONE
   use mem, ONLY: NO_D3_BOX_STATES,NO_BOXES,D3SOURCE,D3STATE, &
                  D3STATETYPE
#ifndef ONESOURCE
   use mem, ONLY: D3SINK
#endif
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
   use mem, ONLY: NO_D2_BOX_STATES,D2SOURCE,D2STATE,NO_BOXES_XY, &
                  D2STATETYPE
#ifndef ONESOURCE
   use mem, ONLY: D2SINK
#endif
#endif
   use standalone
   use api_bfm
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
   real(RLEN),parameter    :: eps=0.0_RLEN
   real(RLEN)              :: min3D,min2D
   integer                 :: i,j,ll
   integer,dimension(2,2)  :: blccc
!
! !EOP
!-----------------------------------------------------------------------
!BOC
   min3D =  ONE
   min2D =  ONE
#ifdef DEBUG
   LEVEL1 'integration RK: starting delt = ',delt
#endif
   bbccc3D=D3STATE
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
   bbccc2D=D2STATE
#endif
   TLOOP : DO
   ! Integration step:
#ifdef ONESOURCE
#ifdef D1SOURCE
      bccc3D=D3SOURCE
#else
      bccc3D=sum(D3SOURCE,2)
#endif
#else
      bccc3D=sum(D3SOURCE-D3SINK,2)
#endif
      ccc_tmp3D=D3STATE
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
#ifdef ONESOURCE
      bccc2D=sum(D2SOURCE,2)
#else
      bccc2D=sum(D2SOURCE-D2SINK,2)
#endif
      ccc_tmp2D=D2STATE
#endif
      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
#ifdef D1SOURCE
            D3STATE(j,:) = ccc_tmp3D(j,:) + delt*D3SOURCE(j,:)
#else
            D3STATE(j,:) = ccc_tmp3D(j,:) + delt*sum(D3SOURCE(j,:,:),1)
#endif
#else
            D3STATE(j,:) = ccc_tmp3D(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         END IF
      END DO
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      DO j=1,NO_D2_BOX_STATES
         IF (D2STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
            D2STATE(j,:) = ccc_tmp2D(j,:) + delt*sum(D2SOURCE(j,:,:),1)
#else
            D2STATE(j,:) = ccc_tmp2D(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
#endif
         END IF
      END DO
#endif
      nmin=nmin+nstep 
      ! Check for negative concentrations
      min3D=minval(D3STATE)
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      min2D=minval(D2STATE)
#endif
      IF(min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
         IF (nstep.eq.1) THEN
            LEVEL1 'Necessary Time Step too small! Exiting...'
            blccc(:,1)=minloc(D3STATE)
            LEVEL1 var_names(stPelStateS+blccc(1,1)-1)
            LEVEL1 ccc_tmp3D(blccc(1,1),blccc(2,1)), &
                           bccc3D(blccc(1,1),blccc(2,1))
            D3STATE=bbccc3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
            blccc(:,2)=minloc(D2STATE)
            LEVEL1 var_names(stBenStateS+blccc(1,2)-1)
            LEVEL1 ccc_tmp2D(blccc(1,2),blccc(2,2)), &
                           bccc2D(blccc(1,2),blccc(2,2))
            D2STATE=bbccc2D
#endif
            LEVEL1 'EXIT at time: ',timesec
            STOP
         END IF
         nstep=nstep/2
         nmin=0
         D3STATE=bbccc3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
         D2STATE=bbccc2D
#endif
         dtm1=maxdelt
         delt=nstep*mindelt
         timesec=ntime*maxdelt
#ifdef DEBUG
         LEVEL2 'Time Step cut! delt= ',delt/2.,' nstep= ',nstep
#endif
         ! Recalculate Sources:
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
      ELSE
#ifdef DEBUG
         LEVEL2 'Internal time step= ',nmin
#endif
         timesec=timesec+delt
#ifdef DEBUG
         LEVEL2 'Internal time= ',timesec
#endif   
         ! Recalculate sources:
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
         DO j=1,NO_D3_BOX_STATES
            IF (D3STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
#ifdef D1SOURCE
               D3STATE(j,:) = ccc_tmp3D(j,:) +.5*delt*(D3SOURCE(j,:)+bccc3D(j,:)) 
#else
               D3STATE(j,:) = ccc_tmp3D(j,:) + &
                  .5*delt*(sum(D3SOURCE(j,:,:),1)+bccc3D(j,:))
#endif
#else
               D3STATE(j,:) = ccc_tmp3D(j,:) + &
                  .5*delt*(sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)+bccc3D(j,:))
#endif
            END IF
         END DO
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
         DO j=1,NO_D2_BOX_STATES
            IF (D2STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
               D2STATE(j,:) = ccc_tmp2D(j,:) + &
                  .5*delt*(sum(D2SOURCE(j,:,:),1)+bccc2D(j,:))
#else
               D2STATE(j,:) = ccc_tmp2D(j,:) + &
                  .5*delt*(sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)+bccc2D(j,:))
#endif
            END IF
         END DO
#endif

         min3D=minval(D3STATE)
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
         min2D=minval(D2STATE)
         IF (min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
#else
         IF (min3D.lt.eps) THEN ! cut timestep
#endif
            IF (nstep.eq.1) THEN
               LEVEL1 'Necessary Time Step too small! Exiting...'
               blccc(:,1)=minloc(D3STATE)
               LEVEL1 var_names(stPelStateS+blccc(1,1)-1)
               LEVEL1 ccc_tmp3D(blccc(1,1),blccc(2,1)), &
                              bccc3D(blccc(1,1),blccc(2,1))
               D3STATE=bbccc3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
               blccc(:,2)=minloc(D2STATE)
               LEVEL1 var_names(stBenStateS+blccc(1,2)-1)
               LEVEL1 ccc_tmp2D(blccc(1,2),blccc(2,2)), &
                              bccc2D(blccc(1,2),blccc(2,2))
               D2STATE=bbccc2D
#endif
               LEVEL1 'EXIT at time: ',timesec
               STOP 'integration-RK2'
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
            LEVEL1 'Time Step cut at RK2! delt= ',delt,' nstep= ',nstep
         ELSE
            IF (nmin.eq.nmaxdelt) EXIT TLOOP
         ENDIF
         ! Recalculate Sources:
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
      END IF
   END DO TLOOP
   nstep=nmaxdelt
   nmin=0
   ntime=ntime+1
   delt=nstep*mindelt
   timesec=delt*ntime
#ifdef DEBUG
   LEVEL1 'ntime: ',ntime
!  LEVEL1 'Integration time: ',time
#endif

   END SUBROUTINE integrationRK2
!EOC
!-----------------------------------------------------------------------
