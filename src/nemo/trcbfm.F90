#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:
!
! !INTERFACE:
   subroutine trc_bfm( kt )
!
! !DESCRIPTION:
!  BFM stepping inside NEMO
!  This routine computes the BFM trends, integrates the benthic 
!  model with a simple Euler Forward and the pelagic model if 
!  a "no transport" simulation is prescribed.
!
! !USES:
   use global_mem, only: RLEN
   use constants,  only: SEC_PER_DAY
   use mem, only: D3STATE,D3SOURCE,NO_D3_BOX_STATES, &
                  D3STATETYPE,NO_BOXES
#ifndef ONESOURCE
   use mem, only: D3SINK
#endif 
#ifdef INCLUDE_BEN
   use mem, only: D2STATE,D2SOURCE,NO_D2_BOX_STATES, &
                  D2STATETYPE
#ifndef ONESOURCE
   use mem, only: D2SINK
#endif 
#endif
   use mem_param, only: CalcTransportFlag, CalcBenthicFlag, &
                        CalcPelagicFlag
   ! NEMO
   use oce_trc          ! ocean dynamics and active tracers variables
   use trc              ! ocean passive tracers variables
   USE prtctl          ! Print control                    (prt_ctl routine)
   
   implicit none
!
! !INPUT PARAMETERS:
!
   integer, intent( IN ) ::  kt  ! ocean time-step index
   
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): M. Vichi
!
! !LOCAL VARIABLES:
   integer               :: j
   integer,save          :: first=0
   real(RLEN)            :: delt
   real(RLEN), allocatable,save,dimension(:,:) :: tmp
   real(RLEN), allocatable,save, dimension(:)  :: tmp2
   integer               :: AllocStatus,DeAllocStatus



!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Accumulation of environmental fields
   !---------------------------------------------
   if ( kt == nittrc000 ) then
      IF( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_bfm : BFM tracer dynamics'
      ENDIF
	  tn_io(:,:,:)=tn(:,:,:)
      sn_io(:,:,:)= sn(:,:,:)
      wn_io(:,:,:)= wn(:,:,:)
      un_io(:,:,:)= un(:,:,:)
      vn_io(:,:,:)= vn(:,:,:)
      avt_io(:,:,:)=avt(:,:,:)
      qsr_io(:,:)= qsr(:,:)
      fr_i_io(:,:)= fr_i(:,:)
      wndm_io(:,:)= wndm(:,:)
      rhop_io(:,:,:)=rhop(:,:,:)
   else 
      vn_io(:,:,:)= vn_io(:,:,:)+vn(:,:,:)      
      un_io(:,:,:)= un_io(:,:,:)+un(:,:,:)      
      tn_io(:,:,:)= tn_io(:,:,:)+tn(:,:,:)      
      sn_io(:,:,:)= sn_io(:,:,:)+sn(:,:,:)      
      wn_io(:,:,:)= wn_io(:,:,:)+wn(:,:,:)      
      avt_io(:,:,:)= avt_io(:,:,:)+avt(:,:,:)      
      qsr_io(:,:)= qsr_io(:,:)+qsr(:,:)      
      fr_i_io(:,:)= fr_i_io(:,:)+fr_i(:,:)      
      wndm_io(:,:)= wndm_io(:,:)+wndm(:,:)      
      rhop_io(:,:,:)= rhop_io(:,:,:)+rhop(:,:,:)  
   end if

   !---------------------------------------------
   ! Proceed only every ndttrc
   !---------------------------------------------   
   if ( MOD( kt , ndttrc ) /= 0 ) return

   !---------------------------------------------
   ! Compute field averages
   !---------------------------------------------   
   tn_io(:,:,:) = tn_io(:,:,:)/ real( ndttrc,RLEN )
   sn_io(:,:,:) = sn_io(:,:,:)/ real( ndttrc,RLEN )
   vn_io(:,:,:) = vn_io(:,:,:)/ real( ndttrc,RLEN )
   un_io(:,:,:) = un_io(:,:,:) / real( ndttrc,RLEN )
   wn_io(:,:,:) = wn_io(:,:,:) / real( ndttrc,RLEN )
   qsr_io(:,:) = qsr_io(:,:)/ real( ndttrc,RLEN )
   avt_io(:,:,:) = avt_io(:,:,:)/ real( ndttrc,RLEN )
   rhop_io(:,:,:) = rhop_io(:,:,:)/ real( ndttrc,RLEN )
   fr_i_io(:,:) = fr_i_io(:,:)/ real( ndttrc,RLEN )
   wndm_io(:,:) = wndm_io(:,:)/ real( ndttrc,RLEN )

   !---------------------------------------------
   ! Biological timestep 
   !---------------------------------------------
   delt  = rdt*real(ndttrc,RLEN)

   !---------------------------------------------
   ! Compute external forcing functions
   !---------------------------------------------
   call envforcing_bfm

   !---------------------------------------------
   ! Compute Biogeochemical trends
   !---------------------------------------------
   call EcologyDynamics

   !---------------------------------------------
   ! Basic ODE solver for benthic variables or for 
   ! 0D test simulations when coupled with NEMO
   !---------------------------------------------
   if (CalcPelagicFlag .AND. .NOT.CalcTransportFlag) then
      do j=1,NO_D3_BOX_STATES
         if (D3STATETYPE(j).ge.0) then
#if defined ONESOURCE && defined D1SOURCE
            D3State(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:),1)
#ifndef D1SOURCE
            D3State(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:),1)
#endif            
#else
            D3State(j,:) = D3STATE(j,:) + delt*sum(D3SOURCE(j,:,:)-D3SINK(j,:,:),1)
#endif
         end if
      end do
   end if
#ifdef INCLUDE_BEN
   if (CalcBenthicFlag /= 0) then
      do j=1,NO_D2_BOX_STATES
         if (D2STATETYPE(j).ge.0) then
#ifdef ONESOURCE
            D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:),1)
#else
            D2STATE(j,:) = D2STATE(j,:) + delt*sum(D2SOURCE(j,:,:)-D2SINK(j,:,:),1)
#endif
         end if
      end do
   end if
#endif
 
   return
   end subroutine trc_bfm

!EOC

