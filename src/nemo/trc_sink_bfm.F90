!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_sink_bfm.F90
!
! !INTERFACE:
   subroutine trc_sink_bfm(ws_in)

! !DESCRIPTION:
!
! This subroutine solves a one-dimensional advection equation for
! the sinking of biogeochemical components. 
! Conservative advection has to be applied when settling of sediment or
! rising of phytoplankton is considered. In this case the advection is of
! the form
!  \begin{equation}
!   \label{Yadvection_cons}
!    \partder{Y}{t} = - \partder{F}{z}
!    \comma
!  \end{equation}
! where $F=wY$ is the flux caused by the advective velocity, $w$.
!
! The discretized form of \eq{Yadvection_cons} is
!  \begin{equation}
!   \label{advDiscretized_cons}
!   Y_i^{n+1} = Y_i^n
!   - \dfrac{\Delta t}{h_i}
!    \left( F^n_{i} - F^n_{i-1} \right)
!   \comma
!  \end{equation}
! where the integers $n$ and $i$ correspond to the present time and space
! level, respectively. 
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
! Several kinds of boundary conditions are implemented for the upper
! and lower boundaries. They are set by the integer values {\tt Bcup}
! and {\tt Bcdw}, that have to correspond to the parameters defined
! in the module {\tt util}, see \sect{sec:utils}. The
! following choices exist at the moment:
!
! For the value {\tt flux}, the boundary values {\tt Yup} and {\tt Ydw} are
! interpreted as specified fluxes at the uppermost and lowest interface.
! Fluxes into the boundary cells are counted positive by convention.
! For the value {\tt value}, {\tt Yup} and {\tt Ydw} specify the value
! of $Y$ at the interfaces, and the flux is computed by multiplying with
! the (known) speed  at the interface. For the value {\tt oneSided},
! {\tt Yup} and {\tt Ydw} are ignored and the flux is computed
! from a one-sided first-order upstream discretisation using the speed
! at the interface and the value of $Y$ at the center of the boundary cell.
! For the value {\tt zeroDivergence}, the fluxes into and out of the
! respective boundary cell are set equal.
! This corresponds to a zero-gradient formulation, or to zero
! flux divergence in the boundary cells.
!
! Be careful that your boundary conditions are mathematically well defined.
! For example, specifying an inflow into the boundary cell with the
! speed at the boundary being directed outward does not make sense.
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
#ifdef BC
!     do the upper boundary conditions
      select case (Bcup)
      case (flux)
         cu(N) = - Yup              ! flux into the domain is positive
      case (value)
         cu(N) =  ww(N)*Yup
      case (oneSided)
         if (ww(N).ge._ZERO_) then
            cu(N) =  ww(N)*Y(N)
         else
            cu(N) = _ZERO_
         end if
      case (zeroDivergence)
         cu(N) = cu(N-1)
      case default
         FATAL 'unkown upper boundary condition type in adv_center()'
         stop
      end select


!     do the lower boundary conditions
      select case (Bcdw)
      case (flux)
         cu(0) =   Ydw               ! flux into the domain is positive
      case (value)
         cu(0) =  ww(0)*Ydw
      case (oneSided)
         if(ww(0).le._ZERO_) then
            cu(0) =  ww(0)*Y(1)
         else
            cu(0) = _ZERO_
         end if
      case (zeroDivergence)
         cu(0) = cu(1)
      case default
         FATAL 'unkown lower boundary condition type in adv_center()'
         stop
      end select
#endif

!     add the vertical advection trend to general tracer trend 
!write(*,*) 'jk, cu(jk+1), cu(jk)' 
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

