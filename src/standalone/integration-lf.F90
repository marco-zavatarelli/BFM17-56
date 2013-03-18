#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Leap Frog time-integration
!
! !INTERFACE
  subroutine integrationLf
!
! !DESCRIPTION
!  Leap Frog time-integration with time-step adjustment
!
! !USES
   use global_mem, ONLY:RLEN
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
   real(RLEN),parameter    :: eps=0.
   real(RLEN),parameter    :: ass=.05
   real(RLEN)              :: min3D,min2D
   integer                 :: i,j,n
   integer,dimension(2,2)  :: blccc
   real(RLEN),dimension(NO_D3_BOX_STATES,NO_BOXES)    :: bc3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
   real(RLEN),dimension(NO_D2_BOX_STATES,NO_BOXES_XY) :: bc2D
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'integration-lf: starting delt = ',delt
#endif
   ! save initial states and additional variables
   bccc3D=D3STATE
   bc3D=bbccc3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
   bccc2D=D2STATE
   bc2D=bbccc2D
#endif

   TLOOP : DO
      ! Integration step:
      DO j=1,NO_D3_BOX_STATES
         IF(D3STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
#ifdef D1SOURCE
            ccc_tmp3D(j,:) = bbccc3D(j,:) + delt * D3SOURCE(j,:)
#else
            ccc_tmp3D(j,:) = bbccc3D(j,:) + delt*sum(D3SOURCE(j,:,:),1)
#endif
#else
            ccc_tmp3D(j,:) = bbccc3D(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         END IF
      END DO
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      DO j=1,NO_D2_BOX_STATES
         IF(D2STATETYPE(j).ge.0) THEN
#ifdef ONESOURCE
                  ccc_tmp2D(j,:) = bbccc2D(j,:) + delt*sum(D2SOURCE(j,:,:),1)
#else
                  ccc_tmp2D(j,:) = bbccc2D(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
#endif
         END IF
      END DO
#endif
      nmin=nmin+nstep 
      ! Check for negative concentrations
      min3D=minval(ccc_tmp3D)
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      min2D=minval(ccc_tmp2D)
      IF (min3D.lt.eps.OR.min2D.lt.eps) THEN ! cut timestep
#else
      IF (min3D.lt.eps) THEN ! cut timestep
#endif
         IF (nstep.eq.1) THEN
            LEVEL1 'Necessary Time Step too small! Exiting...'
               blccc(1,:)=minloc(ccc_tmp3D)
               LEVEL2 ccc_tmp3D(blccc(1,1),blccc(1,2)), &
                           bbccc3D(blccc(1,1),blccc(1,2))
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
               blccc(2,:)=minloc(ccc_tmp2D)
               LEVEL2 ccc_tmp2D(blccc(2,1),blccc(2,2)), &
                           bbccc2D(blccc(2,1),blccc(2,2))
#endif
               LEVEL2 'EXIT at time: ',timesec
            STOP
         END IF
         nstep=nstep/2
         nmin=0
         D3STATE=bccc3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
         D2STATE=bccc2D
#endif
         dtm1=.5*delt
         delt=2.*nstep*mindelt
         timesec=maxdelt*ntime
         n=nmaxdelt/nstep
      ELSE
         ! filtering and advancement:
         write(6,*) 'filtered'
         DO j=1,NO_D3_BOX_STATES
            IF (D3STATETYPE(j).ge.0) THEN
               bbccc3D(j,:) = D3STATE(j,:) + &
                     ass*(bbccc3D(j,:)-2.*D3STATE(j,:)+ccc_tmp3D(j,:))
            END IF
         END DO
         D3STATE=ccc_tmp3D
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
         DO j=1,NO_D2_BOX_STATES
            IF (D2STATETYPE(j).ge.0) THEN
               bbccc2D(j,:) = D2STATE(j,:) + &
                     ass*(bbccc2D(j,:)-2.*D2STATE(j,:)+ccc_tmp2D(j,:))
            END IF
         END DO
         D2STATE=ccc_tmp2D
#endif
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
      END IF
      IF(nmin.eq.nmaxdelt) EXIT TLOOP
      IF(nmin.eq.0) THEN
         ! 2nd order approximation of backward State from Taylor expansion:
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
         bbccc2D=bc2D/n**2+D2STATE*(1.-1./n**2)+ &
            sum((D2SOURCE-D2SINK),2)*.5*delt*(1./n**2-1./n)
#endif
#ifdef ONESOURCE
#ifdef D1SOURCE
         bbccc3D=bc3D/n**2+D3STATE*(1.-1./n**2) + D3SOURCE*.5*delt*(1./n**2-1./n)
#else
         bbccc3D=bc3D/n**2+D3STATE*(1.-1./n**2)+ &
            sum(D3SOURCE,2)*.5*delt*(1./n**2-1./n)
#endif
#else
         bbccc3D=bc3D/n**2+D3STATE*(1.-1./n**2)+ &
            sum((D3SOURCE-D3SINK),2)*.5*delt*(1./n**2-1./n)
#endif
#ifdef DEBUG
         LEVEL2 'Time Step cut! delt= ',delt/2.,' nstep= ',nstep
#endif
      ELSE
#ifdef DEBUG
         LEVEL2 'Internal time step= ',nmin
#endif
         timesec=timesec+nstep*mindelt
#ifdef DEBUG
         LEVEL2 'Internal time= ',timesec
#endif
      END IF
   END DO TLOOP

   IF (nstep.ne.nmaxdelt) THEN 
      ! filter and forwarding of central value if the time step was cut:
#ifdef DEBUG
      LEVEL2 'Full Step Filter after cutting'
#endif
      bbccc3D=bccc3D
      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
            bbccc3D(j,:) = bccc3d(j,:) + &
               ass*(bc3D(j,:)-2.*bccc3d(j,:)+D3STATE(j,:))
         END IF
      END DO
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      bbccc2D=bccc2D
      DO j=1,NO_D2_BOX_STATES
         IF (D2STATETYPE(j).ge.0) THEN
            bbccc2d(j,:) = bccc2d(j,:) + &
               ass*(bc2D(j,:)-2.*bccc2d(j,:)+D2STATE(j,:))
         END IF
      END DO
#endif
   ENDIF
   nstep=nmaxdelt
   nmin=0
   ntime=ntime+1
   dtm1=.5*delt
   delt=2*nstep*mindelt
   timesec=delt/2.*ntime
#ifdef DEBUG
   LEVEL1 'ntime: ',ntime
   LEVEL1 'Integration time: ',timesec
#endif

 END subroutine integrationLf
!EOC
!-----------------------------------------------------------------------
