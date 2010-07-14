#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_sink_muscl_bfm.F90
!
! !INTERFACE:
   SUBROUTINE trc_sink_muscl_bfm(ws_in)
!
! !DESCRIPTION:
!  Compute the sinking terms with an upwind advection scheme
!  (mass conserving) or with MUSCL. 
!
! !USES:
! OPA modules
! ==================
   use oce_trc          ! ocean dynamics and active tracers variables
   use trc              ! ocean passive tracers variables
   use constants, ONLY: DAY_PER_SEC

      IMPLICIT NONE
#include "domzgr_substitute.h90"
!
! !INPUT PARAMETERS:
      REAL(wp), INTENT(IN)  :: ws_in(jpi,jpj,jpk)
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
! Original authors (PISCES) :  A. Estublier, O. Aumont
! BFM adaptation: Marcello Vichi (CMCC-INGV)
!
! !LOCAL VARIABLES:
      INTEGER   :: ji,jj,jk
      REAL(wp)  :: ztraz(jpi,jpj,jpk),zakz(jpi,jpj,jpk)
      REAL(wp)  :: zkz(jpi,jpj,jpk)
      REAL(wp)  :: zigma,zew,zstep,zign,limit,dt
      REAL(wp)  :: ws(jpi,jpj,jpk),sinktemp(jpi,jpj,jpk)
      LOGICAL   :: ln_sink_muscl = .FALSE.
!
!EOP
!-----------------------------------------------------------------------
!BOC

!-------------------------------
! Initialization
!-------------------------------
      zstep  = rdt*FLOAT(ndttrc)
      ztraz  = 0.0_wp
      zkz    = 0.0_wp
      zakz   = 0.0_wp

!-------------------------------
! Set first level velocity to 0
! and shift array downward
! Conversion m/d -> m/s
!-------------------------------
      ws(:,:,1) = 0.0_wp
      DO jk=1,jpk-1
         ws(:,:,jk+1) = ws_in(:,:,jk)*DAY_PER_SEC*tmask(:,:,jk+1)
!#    if defined key_off_degrad
!         ws(:,:,jk+1)=  ws(:,:,jk+1)*facvol(:,:,jk)
!#    endif
      END DO
      
      IF ( ln_sink_muscl ) THEN
        limit = 1.0_wp
      !-------------------------------
      ! First guess of the slopes
      !-------------------------------
        DO jk=2,jpkm1
           ztraz(:,:,jk) = (trn(:,:,jk-1,1) - trn(:,:,jk,1))  &
                               *tmask(:,:,1)
        END DO

      !-------------------------------
      ! Compute slopes
      !-------------------------------
        DO jk=2,jpkm1
          DO jj = 1,jpj
            DO ji = 1, jpi
            zign = 0.5_wp*(sign(1.,ztraz(ji,jj,jk)*ztraz(ji,jj,jk+1))+1)
            zakz(ji,jj,jk) = 0.5_wp*(ztraz(ji,jj,jk)   &
                             +ztraz(ji,jj,jk+1))*zign
            END DO
          END DO
        END DO        
      !-------------------------------
      ! Slope limitation
      !-------------------------------
        DO jk=2,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
              zakz(ji,jj,jk) = sign(1.,zakz(ji,jj,jk)) * &
                              min(abs(zakz(ji,jj,jk)),   &
                              2.*abs(ztraz(ji,jj,jk+1)), &
                              2.*abs(ztraz(ji,jj,jk)))
            END DO
          END DO
        END DO        
      ELSE
        limit = 0.0_wp
      END IF

!-------------------------------
! vertical sinking flux
!-------------------------------
      DO jk=1,jpkm1
         DO jj = 1,jpj      
            DO ji = 1, jpi    
              zigma = ws(ji,jj,jk+1)*zstep/fse3w(ji,jj,jk+1)
              zew   = ws(ji,jj,jk+1)
              sinktemp(ji,jj,jk+1) = -zew*(trn(ji,jj,jk,1)  &
                 -limit*0.5_wp*(1+zigma)*zakz(ji,jj,jk))
            END DO
         END DO
      END DO  

!-------------------------------
! Boundary conditions
!-------------------------------
      sinktemp(:,:,1)   = 0.0_wp
      sinktemp(:,:,jpk) = 0.0_wp

!-------------------------------
! Add vertical sinking trends
!-------------------------------
      DO jk=1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               dt    = rdttra(jk) * FLOAT(ndttrc)
               zigma= &
                                  dt*(sinktemp(ji,jj,jk) - &
                                  sinktemp(ji,jj,jk+1))      &
                                  /fse3t(ji,jj,jk)
               trn(ji,jj,jk,1) = trn(ji,jj,jk,1)+zigma
            END DO
         END DO
      END DO

      trb(:,:,:,1)=trn(:,:,:,1)

   RETURN
   
END SUBROUTINE trc_sink_muscl_bfm
