!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: pom_bfm
!
! !INTERFACE:
       subroutine pom_bfm_1d
!
! !DESCRIPTION:
!  BFM coupling with POM
!  Time-marching routines: trend computations (bio+transport)
!  and integration
!
!
! !USES:
!
       use api_bfm,ONLY           :out_delta
       use Service,ONLY           :ilong,savef
       use constants,ONLY         :SEC_PER_DAY
       use POM,ONLY               :time,time0,dti,intt
       use global_mem, ONLY       :RLEN
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
       IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): M. Butenschoen,  M. Vichi, 
!                      M. Zavatarelli, L. Polimene 
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Local Variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!      -----MODEL TIME IN DAYS-----
!
       real(RLEN),save       :: TT
!
!      -----FLAG FOR TT INITIALISATION-----
!
       logical, save          ::first
       !real(RLEN) :: start1, finish1
       data first /.true./
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Pass physical variables into bfm
! compute extinction coefficients
! compute vertical light distribution
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
       call env_forcing_pom_bfm_1d
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! calculate biological processes 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
       call EcologyDynamics
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Vertical diffusion and Integration of BFM state Vars    
! with source splitting method
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
        !call cpu_time(start1)
        call vdiff_SOS 
        !call cpu_time(finish1)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  leap frog integration of 2d state var's
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
#ifdef INCLUDE_BEN
!
      call lf2d
!
#endif
!
!      -----DEFINE TIME FOR BFM-----
!
       if(first) then
!
           TT=time-time0-(dti/SEC_PER_DAY)
           out_delta=savef
           first=.false.
!
       endif
!
       TT = TT + dti/SEC_PER_DAY
!
!      -----MANAGE OUTPUT-----
!
       call pom_dia_bfm(intt,TT)
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Reset trend arrays
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
       call ResetFluxes
       !print '("Time1 = ",f6.3," seconds.")',finish1-start1
!
       return
!
       end subroutine pom_bfm_1d

!EOC


