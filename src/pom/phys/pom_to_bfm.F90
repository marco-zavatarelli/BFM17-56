#include "INCLUDE.h"
!
! !ROUTINE: Service
!
! DESCRIPTION
!    This subroutine passes the physical variables to the BFM
!
! !INTERFACE
   subroutine pom_to_bfm
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use global_mem,ONLY: RLEN
!
#ifdef NOPOINTERS
!
      use mem
!
#else
!
      use Mem, ONLY: ETW, ESW, EIR, ESS,ERHO,ewind,depth
!
#endif
!
      use POM, ONLY: KB,rcp,TB,SB,RHO,H,DZ,SWRAD,wusurf,wvsurf,KH,KM,U,V,Q2,Q2L,RHO,L
!
      use Service, ONLY: ISM,WGEN,WEDDY
!
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
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  Scalar Arguments
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
      real(RLEN) :: tauw
!
!     -----LOOP COUNTER-----
!
      integer :: k
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! 1D ARRAYS FOR BFM
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      do k = 1 , KB - 1
!
             ETW(k)   = tb(k)
             ESW(k)   = sb(k)
             ERHO(k)  = (rho(k)*1.E3_RLEN)+1.E3_RLEN
             ESS(k)   = ISM(k)
             Depth(k) = dz(k)*h
             
!
      end do
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Surface radiation in deg.C transformed in W.m-2
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      EIR(1) = (-1.0_RLEN)*SWRAD*rcp
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Wind approximate velocity calculation
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      tauw=sqrt(wusurf**2+wvsurf**2)*1.e3_RLEN
      ewind=sqrt(tauw/(1.25_RLEN*0.0014_RLEN))
!
      return
!
      end subroutine pom_to_bfm
!
! EOC
