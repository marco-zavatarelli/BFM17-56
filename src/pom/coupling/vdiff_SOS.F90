!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: vdiff_SOS
!
!DESCRIPTION
!
!  This routine calculates the vertical diffusivity of BFM
!  biochemical components and integrates BFM state var's
!  with Source Splitting (SoS) method.
!
! !INTERFACE
!
  SUBROUTINE vdiff_SOS
!
! !USES:
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
 use global_mem, ONLY:RLEN,ZERO
!
 use Service
!
 use api_bfm, ONLY:D3STATEB
!
 use constants, ONLY: SEC_PER_DAY
!
 use mem_Param, ONLY: AssignAirPelFluxesInBFMFlag,p_small
!
 use mem_Settling
!
 use Mem, ONLY:D3STATE,D3SOURCE,NO_D3_BOX_STATES,                     &
               ppO2o,                                                 &
               ppO3c,                                                 &
               ppN1p,ppN3n,ppN4n,ppN5s,                               &
               ppR1c,ppR1n,ppR1p,                                     &
               ppR6c,ppR6n,ppR6p,ppR6s,                               &
               ppP1c,ppP1n,ppP1p,ppP1s,ppP1l,                         &
               ppP2c,ppP2n,ppP2p,ppP2l,                               &
               ppP3c,ppP3n,ppP3p,ppP3l,                               &
               ppP4c,ppP4n,ppP4p,ppP4l,                               &
               ppZ3c,ppZ3n,ppZ3p,                                     &
               ppZ4c,ppZ4n,ppZ4p,                                     &
               ppZ5c,ppZ5n,ppZ5p,                                     &
               ppZ6c,ppZ6n,ppZ6p,                                     &
               sediR6,sediPPY,                                        &
               iiP1,iiP2,iiP3,iiP4,                                   &
               n1p,n3n,n4n,n5s,                                       &
               jsurO2o,jsurO3c,o2o,                                   &
               Depth,                                                 &
               iiPhytoPlankton
!
 use POM, ONLY:SMOTH,KB,H,DTI,DZR,NRT_O2o,NRT_N1p,NRT_N3n,NRT_N4n,    &
      NBCBFM,UMOLBFM,NTP,TIME
!-------------------------------------------------------------------------!
!
!BOC
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
  IMPLICIT NONE
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Local Variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!-----COUNTER & FLAGS-----
!
 integer                  :: k,m,n,nbc
!
!-----BFM STATE VAR. @ time t-DTI, t, T+DTI RESPECTIVELY-----
!
 real(RLEN)               :: fbio(KB), ffbio(KB), fbbio(KB)
!
!-----SURFACE FLUX STORAGE FOR NUT'S O2 & CO2-----
!
 real(RLEN)               :: surflux
 real(RLEN)               :: botflux
!
!-----RELAXATION VELOCITY FOR NUT'S-----
!
 real(RLEN)               :: trelax_o2o
 real(RLEN)               :: trelax_n1p
 real(RLEN)               :: trelax_n3n
 real(RLEN)               :: trelax_n4n
!
!-----SEDIMENTATION VELOCITY-----
!
 real(RLEN)               :: sink(KB)
 real(RLEN)               :: POCsink
 real(RLEN)               :: W1R6
!
!-----TWICE THE TIME STEP-----
!
 real(RLEN)               :: DTI2
 ! The input general cir. vertical vel. is suppose to be in m/s
 real(RLEN)               :: W_ON = 1.0
! The input eddy vertical vel. is provided in m/d
 real(RLEN)               :: WEddy_ON = 0.1/86400.0 ! to m/s
!
!
    DTI2 = DTI*2.
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Relaxation of nutrients at surface
! -=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
    trelax_o2o=NRT_O2o/SEC_PER_DAY
    trelax_n1p=NRT_N1p/SEC_PER_DAY
    trelax_n3n=NRT_N3n/SEC_PER_DAY
    trelax_n4n=NRT_N4n
!
!-----LOOP OVER BFM STATE VAR'S-----
!
  do m = 1 , NO_D3_BOX_STATES
!
!    -----ZEROING-----
!
      surflux    = ZERO
      botflux    = ZERO
      fbio(:)    = ZERO
      fbbio(:)   = ZERO
      ffbio(:)   = ZERO
      sink(:)    = ZERO
      POCsink    = ZERO
!
!         -----LOAD BFM STATE VAR.-----
!
      do k = 1 , KB - 1

          fbio(k) = D3STATE(m,k)
          fbbio(k) = D3STATEB(m,k)

      end do

      fbio(kb)=fbio(kb-1)
      fbbio(kb)=fbbio(kb-1)

      !do k = 1 , KB - 1
      do k = 1 , KB

          sink(k) = W_ON*WGEN(k) + WEddy_ON*WEDDY(k)

      end do
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Nutrients surface and bottom fluxes:
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      select case ( m )

          case (ppO2o)

              surflux= -(jsurO2o(1)/SEC_PER_DAY)
              botflux = (o2o(kb-1)-O2BOTT)*trelax_o2o

          case (ppO3c)

              surflux= ZERO

          case (ppN1p)

              surflux = ZERO
              botflux = (n1p(kb-1)-PO4BOTT)*trelax_n1p

          case (ppN3n)

              surflux = ZERO
              botflux = (n3n(kb-1)-NO3BOTT)*trelax_n3n

          case (ppN4n)
          ! Ammonium bottom flux is reminerilization from BATS PON data

              surflux = ZERO
              botflux = PONBOTTgrad*trelax_n4n

          case (ppN5s)

              surflux = ZERO

          case default

              surflux = ZERO

      end select
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Bottom Flux
! R1: Dissolved Organic Matter
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
! The botflux for Dissolved Organic Matter is left equal to ZERO
!

! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Bottom Flux and Sedimentation
! R6: Particulate Organic Matter
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
! The botflux for Particulate Organic Matter is left equal to ZERO
!
      select case ( m )

          case (ppR6c:ppR6s)

              do k = 1 , KB - 1

                  sink(k) = sink(k) - sediR6(k)/SEC_PER_DAY

              end do

              ! Final sink value
              sink(kb) = sink(kb-1)

      end select
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Sedimentation Phytoplankton
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
! The botflux for Phytoplankton is left equal to ZERO
!

      if (m.GE.ppP1c .AND. m.LE.ppP4l) then

          select case (m)

              case (ppP1c:ppP1s)

                  n = iiP1

              case (ppP2c:ppP2l)

                  n = iiP2

              case (ppP3c:ppP3l)

                  n = iiP3

              case (ppP4c:ppP4l)

                  n = iiP4

          end select

          do k = 1 , KB - 1

              sink(k) = sink(k) - sediPPY(n,k)/SEC_PER_DAY

          enddo

          ! Final sink value
          sink(kb) = sink(kb-1)

      end if

! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Bottom Flux
! Z: Zooplankton
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
! The bot flux for Zooplankton is left equal to ZERO.
!

! ----- SINKING: UPSTREAM VERTICAL ADVECTION -----

      call adverte(fbbio,fbio,ffbio,sink)

! ----- SOURCE SPLITTING LEAPFROG INTEGRATION -----

      do K=1,KB-1

          ffbio(k)=fbbio(k)+DTI2*((ffbio(k)/H)+D3SOURCE(m,k))

      end do

! ----- COMPUTE VERTICAL DIFFUSION AND TERMINATE INTEGRATION -----
! ----- IMPLICIT LEAPFROGGING -----

      CALL PROFTS(ffbio,surflux,botflux,ZERO,ZERO,NBCBFM,DTI2,NTP,UMOLBFM)

! ----- CLIPPING......IF NEEDED-----

      do k=1, KB-1

          ffbio(k)=max(p_small,ffbio(k))

      end do
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Mix the time step and restore time sequence
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
     do n = 1, KB-1

        D3STATEB(m,n)=fbio(n)+0.5_RLEN*smoth*(ffbio(n)+ fbbio(n)-2.0_RLEN*fbio(n))
        D3STATE(m,n)=ffbio(n)

     enddo

   enddo ! loop over NO_D3_BOX_STATES

   if (.NOT.AssignAirPelFluxesInBFMFlag) then
        jsurO2o(:) = ZERO
        jsurO3c(:) = ZERO
   endif

  end subroutine vdiff_sos
!
!
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: adverte
!
!DESCRIPTION
!
!    SUBROUTINE TO HANDLE THE SINKING OF BFM STATE VAR'S
!    SINKING IS TREATED AS DOWNWARD VERTICAL ADVECTION
!    COMPUTED WITH UPSTREAM FINITE DIFFERENCES.
!
!                                          Maco.Zavatarelli@unibo.it
!
! !INTERFACE
!
  SUBROUTINE adverte(FB,F,FF,W)
!
! !USES:
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use POM, ONLY:KB, DZR, H
!
      use global_mem, ONLY:RLEN
!-------------------------------------------------------------------------!
!
!BOC
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
  IMPLICIT NONE
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Local Variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
 real(RLEN) :: FB(KB),F(KB),FF(KB)
 real(RLEN) :: W(KB)
 real(RLEN) :: DTI2
 integer :: k
!
       F(KB)=F(KB-1)
       FB(KB)=FB(KB-1)
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
! Calculate vertical advection. Mind downward velocities are negative!
! Upwind scheme:
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=
!
      FF(1)=DZR(1)*F(1)*W(2)

      do K=2,KB-1
!
         FF(K)=DZR(K)*(F(K)*W(K+1)-F(K-1)*W(K))
!
      end do
!
      return
!
    end subroutine adverte
!
!EOC
