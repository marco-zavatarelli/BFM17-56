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
   use global_mem, only:RLEN
   use mem_param,  only: AssignAirPelFluxesInBFMFlag,        &
                         AssignPelBenFluxesInBFMFlag
   use mem_PelGlobal, only: p_rR6m
   use mem
   use constants,    only: SEC_PER_DAY
   use mem_settling, only: p_burvel_R6,p_burvel_R2,p_burvel_PI
   use api_bfm

   implicit none
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
   integer,parameter     :: KSINK=20 ! set to jpk to exclude
   real(RLEN),parameter  :: depth_factor = 2000.0_RLEN
   real(RLEN)            :: zfact,timestep,wsmax
   real(RLEN)            ::  wbio(jpi,jpj,jpk)   
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
   select case (m)
#ifdef INCLUDE_PELFE
      case (ppP1c,ppP1n,ppP1p,ppP1s,ppP1l,ppP1f)
#else
      case (ppP1c,ppP1n,ppP1p,ppP1s,ppP1l)
#endif
#ifdef USEPACK
         wbio = -unpack(sediPPY(iiP1,:),SEAmask,ZEROS)
#else
         DO n = 1,NO_BOXES
            wbio(iwet(n),jwet(n),kwet(n)) = -sediPPY(iiP1,n)
         END DO
#endif
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         ! Prescribe sinking velocity below 
         ! level KSINK (usually > 150 m)
         !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
         wbio(:,:,KSINK:jpk) = -p_rR6m
         CALL trc_sink_bfm(wbio)       ! vertical sinking
#ifdef INCLUDE_PELFE
      case (ppR6c,ppR6n,ppR6p,ppR6s,ppR6f)
#else
      case (ppR6c,ppR6n,ppR6p,ppR6s)
#endif
#ifdef USEPACK
         wbio = -unpack(sediR6(:),SEAmask,ZEROS)
#else
         DO n = 1,NO_BOXES
            wbio(iwet(n),jwet(n),kwet(n)) = -sediR6(n)
         END DO
#endif
         CALL trc_sink_bfm(wbio)       ! vertical sinking
      case default
         wbio = 0.0_RLEN
   end select

#ifdef AGGREGATION
   !---------------------------------------------
   ! Sinking speeds increase with depth below 
   ! the turbocline depth (aggregation)
   ! Velocity is limited according to the depth 
   ! of the layer
   !---------------------------------------------
   do jk=1,jpk-1
      do jj=1,jpj
         do ji=1,jpi
            wsmax=0.8*fse3t(ji,jj,jk)/timestep
            zfact = max(0.0_RLEN,exp((fsdepw(ji,jj,jk)-hmld(ji,jj)) &
                                     /depth_factor)-1.0_RLEN);
            wbio(ji,jj,jk) = min(wsmax,(1.0_RLEN+zfact)*wbio(ji,jj,jk))

          end do
      end do
   end do
#endif

   return
   end subroutine trc_set_bfm

!EOC

