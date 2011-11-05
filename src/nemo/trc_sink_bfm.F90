!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_sink_bfm.F90
!
! !INTERFACE:
   subroutine trc_sink_bfm(ws_in)

! !DESCRIPTION:
!
! This subroutine computes the RHS of a one-dimensional advection equation
! for the sinking of biogeochemical components. 
! Conservative advection has to be applied when settling of sediment or
! rising of phytoplankton is considered. 
!
! Fluxes are defined at the grid faces, the variable $Y_i$ is defined at the
!  grid centers. The fluxes are computed in an upstream-biased way,
!  \begin{equation}
!   \label{upstream}
!   F^n_{i} = \dfrac{1}{\Delta t}
!   \int_{z^\text{Face}_{i} - w \Delta t}^{z^\text{Face}_{i}} Y(z') dz'
!   \point
!  \end{equation}
!
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
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf (GOTM team)
!  BFM-NEMO adaptation: Marcello Vichi (CMCC-INGV)
!
!EOP
!
! !LOCAL VARIABLES:
      integer   :: ji,jj,jk
      real(wp)  :: Yu,Yd,Yc,Dc,Dd
      real(wp)  :: ws(jpi,jpj,jpk),cu(jpi,jpj,jpk)
!
!-----------------------------------------------------------------------
!BOC

!  initialize interface fluxes with zero
   cu   = 0.0_wp
!-------------------------------
! Set first level velocity to 0
! and shift array downward
! Conversion m/d -> m/s
!-------------------------------
   ws(:,:,:) = 0.0_wp
   do jk=1,jpk-1
      ws(:,:,jk+1) = ws_in(:,:,jk)*DAY_PER_SEC*tmask(:,:,jk+1)
   end do

   do jj = 1, jpj      
      do ji = 1, jpi    
         do jk = 2, jpkm1
            !-------------------------------
            ! Check the velocity direction
            !-------------------------------
            if (ws(ji,jj,jk) .gt. 0.0_wp) then
               ! positive speed
               Yu=trn(ji,jj,jk+1,1)          ! upstream value
               Yc=trn(ji,jj,jk,1)            ! central value
               if (jk .gt. 1) then
                  Yd=trn(ji,jj,jk-1,1)       ! downstream value
               else
                  Yd=trn(ji,jj,1,1) 
               end if
            else
               ! negative speed
               if (jk .gt. 2) then
                  Yu=trn(ji,jj,jk-2,1)       ! upstream value
               else
                  Yu=trn(ji,jj,1,1)          ! upstream value
               end if
               Yc=trn(ji,jj,jk-1,1)          ! central value
               Yd=trn(ji,jj,jk,1)            ! downstream value
            end if
            ! compute the mass flow across the central interface 
            cu(ji,jj,jk)=ws(ji,jj,jk)*Yc
         end do
      end do
   end do

   cu(:,:,1) = 0.0_wp
   cu(:,:,jpk) = 0.0_wp

! add the vertical advection trend to general tracer trend 
   Yc=0.0_wp
   Yd=0.0_wp
   do jj = 1, jpj      
      do ji = 1, jpi    
         do jk = 1, jpkm1
            tra(ji,jj,jk,1) = tra(ji,jj,jk,1) + &
                 (cu(ji,jj,jk+1)-cu(ji,jj,jk))/fse3t(ji,jj,jk)
         end do
      end do 
   end do 

   return
   end subroutine trc_sink_bfm
!EOC

