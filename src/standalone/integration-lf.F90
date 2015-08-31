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
   use global_mem, ONLY:RLEN, ONE
   use mem, ONLY: NO_D3_BOX_STATES,NO_BOXES,D3SOURCE,D3STATE, &
                  D3STATETYPE
#ifdef EXPLICIT_SINK
   use mem, ONLY: D3SINK
#endif

#if defined INCLUDE_SEAICE
   use mem, ONLY: NO_D2_BOX_STATES_ICE, D2SOURCE_ICE, &
        NO_BOXES_XY, D2STATE_ICE, D2STATETYPE_ICE
#ifdef EXPLICIT_SINK
   use mem, ONLY: D2SINK_ICE
#endif
#endif

#if defined INCLUDE_BEN
   use mem, ONLY: NO_D2_BOX_STATES_BEN, D2SOURCE_BEN, &
        NO_BOXES_XY, D2STATE_BEN, D2STATETYPE_BEN
#ifdef EXPLICIT_SINK
   use mem, ONLY: D2SINK_BEN
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
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
   real(RLEN),parameter    :: ass=.05_RLEN
   real(RLEN)              :: min3D,min2D
   logical                 :: min
   integer                 :: i,j,n
   integer,dimension(2,2)  :: blccc
   real(RLEN),dimension(NO_D3_BOX_STATES,NO_BOXES)    :: bc3D
#if defined INCLUDE_SEAICE
   real(RLEN)                                                 :: min2D_ice
   integer,dimension(2,2)                                     :: blccc_ice
   real(RLEN),dimension(NO_D2_BOX_STATES_ICE,NO_BOXES_XY) :: bc2D_ice
#endif
#if defined INCLUDE_BEN
   real(RLEN)                                                 :: min2D_ben
   integer,dimension(2,2)                                     :: blccc_ben
   real(RLEN),dimension(NO_D2_BOX_STATES_BEN,NO_BOXES_XY) :: bc2D_ben
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

#if defined INCLUDE_SEAICE
   bccc2D_ice=D2STATE_ICE
   bc2D_ice=bbccc2D_ice
#endif

#if defined INCLUDE_BEN
   bccc2D_ben=D2STATE_BEN
   bc2D_ben=bbccc2D_ben
#endif

   TLOOP : DO
      ! Integration step:
      DO j=1,NO_D3_BOX_STATES
         IF(D3STATETYPE(j).ge.0) THEN
#ifndef EXPLICIT_SINK
            ccc_tmp3D(j,:) = bbccc3D(j,:) + delt * D3SOURCE(j,:)
#else
            ccc_tmp3D(j,:) = bbccc3D(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         END IF
      END DO

#if defined INCLUDE_SEAICE
      DO j=1,NO_D2_BOX_STATES_ICE
         IF(D2STATETYPE_ICE(j).ge.0) THEN
#ifndef EXPLICIT_SINK
                  ccc_tmp2D_ice(j,:) = bbccc2D_ice(j,:) + delt*D2SOURCE_ICE(j,:)
#else
                  ccc_tmp2D_ice(j,:) = bbccc2D_ice(j,:) + delt*sum(D2SOURCE_ICE(j,:,:)-D2SINK_ICE(j,:,:),1)
#endif
         END IF
      END DO
#endif

#if defined INCLUDE_BEN
      DO j=1,NO_D2_BOX_STATES_BEN
         IF(D2STATETYPE_BEN(j).ge.0) THEN
#ifndef EXPLICIT_SINK
                  ccc_tmp2D_ben(j,:) = bbccc2D_ben(j,:) + delt*D2SOURCE_BEN(j,:)
#else
                  ccc_tmp2D_ben(j,:) = bbccc2D_ben(j,:) + delt*sum(D2SOURCE_BEN(j,:,:)-D2SINK_BEN(j,:,:),1)
#endif
         END IF
      END DO
#endif

      nmin=nmin+nstep 
      ! Check for negative concentrations
      min3D=minval(ccc_tmp3D)
      min =  ( min3D.lt.eps )
#if defined INCLUDE_SEAICE
      min2D_ice=minval(ccc_tmp2D_ice)
      min = ( min .OR. ( min2D_ice .lt. eps ) )
#elif defined INCLUDE_BEN
      min2D_ben=minval(ccc_tmp2D_ben)
      min = ( min .OR. ( min2D_ben .lt. eps ) )
#endif
      IF ( min ) THEN ! cut timestep
         IF (nstep.eq.1) THEN
            LEVEL1 'Necessary Time Step too small! Exiting...'
               blccc(1,:)=minloc(ccc_tmp3D)
               LEVEL2 ccc_tmp3D(blccc(1,1),blccc(1,2)), &
                           bbccc3D(blccc(1,1),blccc(1,2))
#if defined INCLUDE_SEAICE
               blccc_ice(2,:)=minloc(ccc_tmp2D_ice)
               LEVEL2 ccc_tmp2D_ice(blccc_ice(2,1),blccc_ice(2,2)), &
                           bbccc2D_ice(blccc_ice(2,1),blccc_ice(2,2))
#endif
#if defined INCLUDE_BEN
               blccc_ben(2,:)=minloc(ccc_tmp2D_ben)
               LEVEL2 ccc_tmp2D_ben(blccc_ben(2,1),blccc_ben(2,2)), &
                           bbccc2D_ben(blccc_ben(2,1),blccc_ben(2,2))
#endif
               LEVEL2 'EXIT at time: ',timesec
            STOP
         END IF
         nstep=nstep/2_RLEN
         nmin=0
         D3STATE=bccc3D
#if defined INCLUDE_SEAICE
         D2STATE_ICE=bccc2D_ice
#endif
#if defined INCLUDE_BEN
         D2STATE_BEN=bccc2D_ben
#endif
         dtm1=.5_RLEN*delt
         delt=2._RLEN*nstep*mindelt
         timesec=maxdelt*ntime
         n=nmaxdelt/nstep
      ELSE
         ! filtering and advancement:
         write(6,*) 'filtered'
         DO j=1,NO_D3_BOX_STATES
            IF (D3STATETYPE(j).ge.0) THEN
               bbccc3D(j,:) = D3STATE(j,:) + &
                     ass*(bbccc3D(j,:)-2._RLEN*D3STATE(j,:)+ccc_tmp3D(j,:))
            END IF
         END DO
         D3STATE=ccc_tmp3D
#if defined INCLUDE_SEAICE
         DO j=1,NO_D2_BOX_STATES_ICE
            IF (D2STATETYPE_ICE(j).ge.0) THEN
               bbccc2D_ice(j,:) = D2STATE_ICE(j,:) + &
                     ass*(bbccc2D_ice(j,:)-2._RLEN*D2STATE_ICE(j,:)+ccc_tmp2D_ice(j,:))
            END IF
         END DO
         D2STATE_ICE=ccc_tmp2D_ice
#endif
#if defined INCLUDE_BEN
         DO j=1,NO_D2_BOX_STATES_BEN
            IF (D2STATETYPE_BEN(j).ge.0) THEN
               bbccc2D_ben(j,:) = D2STATE_BEN(j,:) + &
                     ass*(bbccc2D_ben(j,:)-2._RLEN*D2STATE_BEN(j,:)+ccc_tmp2D_ben(j,:))
            END IF
         END DO
         D2STATE_BEN=ccc_tmp2D_ben
#endif
         call ResetFluxes
         call envforcing_bfm
         call EcologyDynamics
      END IF
      IF(nmin.eq.nmaxdelt) EXIT TLOOP
      IF(nmin.eq.0) THEN
         ! 2nd order approximation of backward State from Taylor expansion:
#if defined INCLUDE_SEAICE
#ifndef EXPLICIT_SINK
         bbccc2D_ice=bc2D_ice/n**2+D2STATE_ICE*(ONE-ONE/n**2) + &
            D2SOURCE_ICE*.5_RLEN*delt*(ONE/n**2 - ONE/n)
#else
         bbccc2D_ice=bc2D_ice/n**2+D2STATE_ICE*(ONE-ONE/n**2)+ &
            sum((D2SOURCE_ICE-D2SINK_ICE),2)*.5_RLEN*delt*(ONE/n**2 - ONE/n)
#endif
#endif

#if defined INCLUDE_BEN
#ifndef EXPLICIT_SINK
         bbccc2D_ben=bc2D_ben/n**2+D2STATE_BEN*(ONE-ONE/n**2) + &
            D2SOURCE_BEN*.5_RLEN*delt*(ONE/n**2 - ONE/n)
#else
         bbccc2D_ben=bc2D_ben/n**2+D2STATE_BEN*(ONE-ONE/n**2)+ &
            sum((D2SOURCE_BEN-D2SINK_BEN),2)*.5_RLEN*delt*(ONE/n**2 - ONE/n)
#endif
#endif

#ifdef DEBUG
         LEVEL2 'Time Step cut! delt= ',delt/2._RLEN,' nstep= ',nstep
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
               ass*(bc3D(j,:)-2._RLEN*bccc3d(j,:)+D3STATE(j,:))
         END IF
      END DO

#if defined INCLUDE_SEAICE
      bbccc2D_ice=bccc2D_ice
      DO j=1,NO_D2_BOX_STATES_ICE
         IF (D2STATETYPE_ICE(j).ge.0) THEN
            bbccc2d_ice(j,:) = bccc2d_ice(j,:) + &
               ass*(bc2D_ice(j,:)-2._RLEN*bccc2d_ice(j,:)+D2STATE_ICE(j,:))
         END IF
      END DO
#endif

#if defined INCLUDE_BEN
      bbccc2D_ben=bccc2D_ben
      DO j=1,NO_D2_BOX_STATES_BEN
         IF (D2STATETYPE_BEN(j).ge.0) THEN
            bbccc2d_ben(j,:) = bccc2d_ben(j,:) + &
               ass*(bc2D_ben(j,:)-2._RLEN*bccc2d_ben(j,:)+D2STATE_BEN(j,:))
         END IF
      END DO
#endif

   ENDIF
   nstep=nmaxdelt
   nmin=0
   ntime=ntime+1
   dtm1=.5_RLEN*delt
   delt=2*nstep*mindelt
   timesec=delt/2._RLEN*ntime
#ifdef DEBUG
   LEVEL1 'ntime: ',ntime
   LEVEL1 'Integration time: ',timesec
#endif

 END subroutine integrationLf
!EOC
!-----------------------------------------------------------------------
