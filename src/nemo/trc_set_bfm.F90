!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_set_bfm.F90
!
! !INTERFACE:
   subroutine trc_set_bfm(kt,m)
!
! !DESCRIPTION:
!  Computes additional boundary conditions and transfer 
!  sinking velocity
!
! !USES:
   ! NEMO
   use oce_trc          ! ocean dynamics and active tracers variables
   use trc              ! ocean passive tracers variables
   ! BFM
   use global_mem, only:RLEN,ZERO,ONE
   use mem_param,  only: AssignAirPelFluxesInBFMFlag,        &
                         AssignPelBenFluxesInBFMFlag
   use mem_PelGlobal, only: p_rR6m, KSINK_rPPY,              &
                            AggregateSink, depth_factor
   use mem
   use constants,    only: SEC_PER_DAY
   use mem_settling, only: p_burvel_R6,p_burvel_R2,p_burvel_PI
   use api_bfm

   implicit none
! OPA domain substitutions
#include "domzgr_substitute.h90"
!
! !INPUT PARAMETERS:
   integer, intent(IN)     ::  kt  ! ocean time-step index
   integer, intent(IN)     ::  m   ! BFM variable index
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): Marcello Vichi (CMCC-INGV)
!  Sinking velocity formulation: O. Aumont (PISCES model)
!
! !LOCAL VARIABLES:
   ! 3D sinking velocity field
   integer               :: ji, jj, jk,n
   real(RLEN)            :: zfact,timestep,wsmax
   real(RLEN)            :: wbio(jpi,jpj,jpk)   
   logical               :: dosink
!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Biological timestep (in days)
   !---------------------------------------------
   timestep  = rdt*FLOAT(nn_dttrc)/SEC_PER_DAY

   !---------------------------------------------
   ! Transfer sinking velocities 
   ! (negative, z-axis is positive upwards)
   !---------------------------------------------
   dosink = .TRUE. ! sinking is initially set to true
   select case (m)
      !
      ! Phytoplankton group with silicate (namely, diatoms)
#ifdef INCLUDE_PELFE
      case (ppP1c,ppP1n,ppP1p,ppP1s,ppP1l,ppP1f)
#else
      case (ppP1c,ppP1n,ppP1p,ppP1s,ppP1l)
#endif
        !---------------------------------------------
        ! Prescribe sinking velocity below depth 
        ! threshold KSINK_rPPY (usually below 150m)
        ! This accelerate diatoms sinking to balance 
        ! the reduction occruing in deeper layers
        ! where nutrient concentration is high
        ! It is alternatevly done by AggregationSink
        !---------------------------------------------
         if ( KSINK_rPPY > 0 .AND. .NOT. AggregateSink ) &
            where( EPR > KSINK_rPPY ) sediPPY(iiP1,:) = p_rR6m
#ifndef NOPACK
         wbio = -unpack(sediPPY(iiP1,:),SEAmask,ZEROS)
#else
         do n = 1,NO_BOXES
            wbio(iwet(n),jwet(n),kwet(n)) = -sediPPY(iiP1,n)
         end do
#endif
      !
      ! Detritus
#ifdef INCLUDE_PELFE
      case (ppR6c,ppR6n,ppR6p,ppR6s,ppR6f)
#else
      case (ppR6c,ppR6n,ppR6p,ppR6s)
#endif
#ifndef NOPACK
         wbio = -unpack(sediR6(:),SEAmask,ZEROS)
#else
         DO n = 1,NO_BOXES
            wbio(iwet(n),jwet(n),kwet(n)) = -sediR6(n)
         END DO
#endif
      case default
         dosink = .FALSE. ! sinking is disabled 
         wbio = 0.0_RLEN
   end select
   !
   !---------------------------------------------
   ! Sinking speeds increase with depth below 
   ! the turbocline depth (aggregation).
   ! Velocity is limited according to the depth 
   ! of the layer (80% of it).
   ! This modulates the assigned sinking rate
   !---------------------------------------------
   !
   if ( AggregateSink .and. dosink) then
      do jk=1,jpk-1
         do jj=1,jpj
            do ji=1,jpi
               wsmax=0.8_RLEN*fse3t(ji,jj,jk)/timestep
               zfact = max(ZERO,exp((fsdepw(ji,jj,jk)-hmld(ji,jj)) &
                                        /depth_factor)-ONE)
               wbio(ji,jj,jk) = min(wsmax,(ONE+zfact)*wbio(ji,jj,jk))

             end do
         end do
      end do
   endif

   !---------------------------------------------
   ! Compute vertical sinking with upwind scheme
   !---------------------------------------------
   if (dosink)  CALL trc_sink_bfm(wbio)

   return
   end subroutine trc_set_bfm

!EOC

