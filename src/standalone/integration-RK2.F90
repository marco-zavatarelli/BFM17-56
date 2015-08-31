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
#ifdef EXPLICIT_SINK
   use mem, ONLY: D3SINK
#endif

#if defined INCLUDE_SEAICE
   use mem, ONLY: NO_D2_BOX_STATES_ICE, D2SOURCE_ICE, &
        D2STATE_ICE, D2STATETYPE_ICE
#ifdef EXPLICIT_SINK
   use mem, ONLY: D2SINK_ICE
#endif
#endif

#if defined INCLUDE_BEN
   use mem, ONLY: NO_D2_BOX_STATES_BEN, D2SOURCE_BEN, &
        D2STATE_BEN, D2STATETYPE_BEN
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
   real(RLEN),parameter    :: eps=0.0_RLEN
   real(RLEN)              :: min3D,min2D
   logical                 :: cut,cut_pel
   integer                 :: i,j,ll
   integer,dimension(2,2)  :: blccc
#if defined INCLUDE_SEAICE
   real(RLEN)              :: min2D_ice
   integer,dimension(2,2)  :: blccc_ice
   logical                 :: cut_ice
#endif
#if defined INCLUDE_BEN
   real(RLEN)              :: min2D_ben
   integer,dimension(2,2)  :: blccc_ben
   logical                 :: cut_ben
#endif
!
! !EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'integration RK: starting delt = ',delt
#endif
   bbccc3D=D3STATE
#if defined INCLUDE_SEAICE
   min2D_ice =  ONE
   bbccc2D_ice=D2STATE_ICE
#endif
#if defined INCLUDE_BEN
   min2D_ben =  ONE
   bbccc2D_ben=D2STATE_BEN
#endif
   TLOOP : DO
   ! Integration step:
#ifndef EXPLICIT_SINK
      bccc3D=D3SOURCE
#else
      bccc3D=sum(D3SOURCE-D3SINK,2)
#endif
      ccc_tmp3D=D3STATE

#if defined INCLUDE_SEAICE
#ifndef EXPLICIT_SINK
      bccc2D_ice=D2SOURCE_ICE
#else
      bccc2D_ice=sum(D2SOURCE_ICE-D2SINK_ICE,2)
#endif
      ccc_tmp2D_ice=D2STATE_ICE
#endif

#if defined INCLUDE_BEN
#ifndef EXPLICIT_SINK
      bccc2D_ben=D2SOURCE_BEN
#else
      bccc2D_ben=sum(D2SOURCE_BEN-D2SINK_BEN,2)
#endif
      ccc_tmp2D_ben=D2STATE_BEN
#endif

      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
#ifndef EXPLICIT_SINK
            D3STATE(j,:) = ccc_tmp3D(j,:) + delt*D3SOURCE(j,:)
#else
            D3STATE(j,:) = ccc_tmp3D(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         END IF
      END DO

#if defined INCLUDE_SEAICE
      DO j=1,NO_D2_BOX_STATES_ICE
         IF (D2STATETYPE_ICE(j).ge.0) THEN
#ifndef EXPLICIT_SINK
            D2STATE_ICE(j,:) = ccc_tmp2D_ice(j,:) + delt*D2SOURCE_ICE(j,:)
#else
            D2STATE_ICE(j,:) = ccc_tmp2D_ice(j,:) + delt*sum(D2SOURCE_ICE(j,:,:)-D2SINK_ICE(j,:,:),1)
#endif
         END IF
      END DO
#endif

#if defined INCLUDE_BEN
      DO j=1,NO_D2_BOX_STATES_BEN
         IF (D2STATETYPE_BEN(j).ge.0) THEN
#ifndef EXPLICIT_SINK
            D2STATE_BEN(j,:) = ccc_tmp2D_ben(j,:) + delt*D2SOURCE_BEN(j,:)
#else
            D2STATE_BEN(j,:) = ccc_tmp2D_ben(j,:) + delt*sum(D2SOURCE_BEN(j,:,:)-D2SINK_BEN(j,:,:),1)
#endif
         END IF
      END DO
#endif

      nmin=nmin+nstep 
      ! Check for negative concentrations
      min3D=minval(D3STATE)
      cut_pel = ( min3D.lt.eps )
      cut = cut_pel
#if defined INCLUDE_SEAICE
      min2D_ice=minval(D2STATE_ICE)
      cut_ice = ( min2D_ice .lt. eps ) 
      cut = ( cut .OR. cut_ice )
#endif
#if defined INCLUDE_BEN
      min2D_ben=minval(D2STATE_BEN)
      cut_ben = ( min2D_ben .lt. eps )
      cut = ( cut .OR. cut_ben )
#endif
      IF( cut ) THEN ! cut timestep
         IF (nstep.eq.1) THEN
            LEVEL1 'Necessary Time Step too small! Exiting from the following systems:'
            if ( cut_pel) then
               blccc(:,1)=minloc(D3STATE)
               LEVEL1 'Pelagic variable: ',var_names(stPelStateS+blccc(1,1)-1)
               LEVEL1 'Value=',ccc_tmp3D(blccc(1,1),blccc(2,1)), &
                              'Rate=',bccc3D(blccc(1,1),blccc(2,1))
               D3STATE=bbccc3D
            end if

#if defined INCLUDE_SEAICE
            if ( cut_ice) then
               blccc_ice(:,2)=minloc(D2STATE_ICE)
               LEVEL1 'Sea ice variable: ', var_names(stIceStateS+blccc_ice(1,2)-1)
               LEVEL1 'Value=',ccc_tmp2D_ice(blccc_ice(1,2),blccc_ice(2,2)), &
                              'Rate=',bccc2D_ice(blccc_ice(1,2),blccc_ice(2,2))
               D2STATE_ICE=bbccc2D_ice
            end if
#endif

#if defined INCLUDE_BEN
            if ( cut_ben) then
               blccc_ben(:,2)=minloc(D2STATE_BEN)
               LEVEL1 'Benthic Variable: ', var_names(stBenStateS+blccc_ben(1,2)-1)
               LEVEL1 'Value=',ccc_tmp2D_ben(blccc_ben(1,2),blccc_ben(2,2)), &
                              'Rate=',bccc2D_ben(blccc_ben(1,2),blccc_ben(2,2))
               D2STATE_BEN=bbccc2D_ben
            end if
#endif

            LEVEL1 'EXIT at time: ',timesec
            STOP 'integration-RK'
         END IF
         nstep=nstep/2
         nmin=0
         D3STATE=bbccc3D

#if defined INCLUDE_SEAICE
         D2STATE_ICE=bbccc2D_ice
#endif

#if defined INCLUDE_BEN
         D2STATE_BEN=bbccc2D_ben
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
#ifndef EXPLICIT_SINK
               D3STATE(j,:) = ccc_tmp3D(j,:) +.5*delt*(D3SOURCE(j,:)+bccc3D(j,:)) 
#else
               D3STATE(j,:) = ccc_tmp3D(j,:) + &
                  .5*delt*(sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)+bccc3D(j,:))
#endif
            END IF
         END DO

#if defined INCLUDE_SEAICE
         DO j=1,NO_D2_BOX_STATES_ICE
            IF (D2STATETYPE_ICE(j).ge.0) THEN
#ifndef EXPLICIT_SINK
               D2STATE_ICE(j,:) = ccc_tmp2D_ice(j,:) + &
                  .5*delt*(D2SOURCE_ICE(j,:)+bccc2D_ice(j,:))
#else
               D2STATE_ICE(j,:) = ccc_tmp2D_ice(j,:) + &
                  .5*delt*(sum(D2SOURCE_ICE(j,:,:)-D2SINK_ICE(j,:,:),1)+bccc2D_ice(j,:))
#endif
            END IF
         END DO
#endif

#if defined INCLUDE_BEN
         DO j=1,NO_D2_BOX_STATES_BEN
            IF (D2STATETYPE_BEN(j).ge.0) THEN
#ifndef EXPLICIT_SINK
               D2STATE_BEN(j,:) = ccc_tmp2D_ben(j,:) + &
                  .5*delt*(D2SOURCE_BEN(j,:)+bccc2D_ben(j,:))
#else
               D2STATE_BEN(j,:) = ccc_tmp2D_ben(j,:) + &
                  .5*delt*(sum(D2SOURCE_BEN(j,:,:)-D2SINK_BEN(j,:,:),1)+bccc2D_ben(j,:))
#endif
            END IF
         END DO
#endif


         ! check for negative values
         min3D=minval(D3STATE)
         cut_pel = ( min3D.lt.eps )
         cut = cut_pel
#if defined INCLUDE_SEAICE
         min2D_ice=minval(D2STATE_ICE)
         cut_ice = ( min2D_ice .lt. eps ) 
         cut = ( cut .OR. cut_ice )
#endif
#if defined INCLUDE_BEN
         min2D_ben=minval(D2STATE_BEN)
         cut_ben = ( min2D_ben .lt. eps )
         cut = ( cut .OR. cut_ben )
#endif
         IF( cut ) THEN ! cut timestep
            IF (nstep.eq.1) THEN
               LEVEL1 'Necessary Time Step too small! Exiting from the following systems:'
               if ( cut_pel) then
                  blccc(:,1)=minloc(D3STATE)
                  LEVEL1 'Pelagic variable: ',var_names(stPelStateS+blccc(1,1)-1)
                  LEVEL1 'value=',ccc_tmp3D(blccc(1,1),blccc(2,1)), &
                                 'rate=',bccc3D(blccc(1,1),blccc(2,1))
                  D3STATE=bbccc3D
               end if

#if defined INCLUDE_SEAICE
               if ( cut_ice) then
                  blccc_ice(:,2)=minloc(D2STATE_ICE)
                  LEVEL1 'Sea ice variable: ', var_names(stIceStateS+blccc_ice(1,2)-1)
                  LEVEL1 'value=',ccc_tmp2D_ice(blccc_ice(1,2),blccc_ice(2,2)), &
                                 'rate=',bccc2D_ice(blccc_ice(1,2),blccc_ice(2,2))
                  D2STATE_ICE=bbccc2D_ice
               end if
#endif

#if defined INCLUDE_BEN
               if ( cut_ben) then
                  blccc_ben(:,2)=minloc(D2STATE_BEN)
                  LEVEL1 'Benthic variable: ', var_names(stBenStateS+blccc_ben(1,2)-1)
                  LEVEL1 'value=',ccc_tmp2D_ben(blccc_ben(1,2),blccc_ben(2,2)), &
                                 'rate=',bccc2D_ben(blccc_ben(1,2),blccc_ben(2,2))
                  D2STATE_BEN=bbccc2D_ben
               end if
#endif

               LEVEL1 'EXIT at time: ',timesec
               STOP 'integration-RK'
            END IF

            nstep=nstep/2
            nmin=0
            D3STATE=bbccc3D

#if defined INCLUDE_SEAICE
            D2STATE_ICE=bbccc2D_ice
#endif

#if defined INCLUDE_BEN
            D2STATE_BEN=bbccc2D_ben
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
