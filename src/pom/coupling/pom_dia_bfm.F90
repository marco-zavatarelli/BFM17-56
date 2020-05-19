#include "INCLUDE.h"
#include "cppdefs.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: pom_dia_bfm
!
!DESCRIPTION    
!  This routine calculates means and writes the output in diagnostic mode
!  writes also the restart
!
! !INTERFACE
   subroutine pom_dia_bfm(kt,TT)
!
! !USES:
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   
   use global_mem, only:RLEN
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use netcdf_bfm, only: save_rst_bfm, ncid_rst
   use api_bfm, only: out_delta
   use constants,  only:SEC_PER_DAY
   use pom, ONLY: time, DTI,IEND
!
!-------------------------------------------------------------------------!
!BOC
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
  IMPLICIT NONE
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!-----TIME COUNTER IN DAYS-----
!
   real(RLEN),intent(INOUT)    ::  TT 
!
!-----TIME MARCHING LOOP COUNTER-----
!
   integer, intent(in)      :: kt
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Local Variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!  -----TIME IN SECONDS-----
!
   real(RLEN)               :: localtime 
!
!  -----SAVING FREQUENCY-----
!
   integer                  :: time_to_save !time in seconds
!
!-----------------------------------------------------------------------
!
!BOC
!
! -----TIME ELLAPSED (IN SECONDS) SINCE THE BEGINNING OF THE SIMULATION----
!
   localtime = time*SEC_PER_DAY   
!
! -----SAVING FREQUENCY IN TIME MARCHING LOOP ITERATIONS-----
!
   time_to_save=nint(out_delta*SEC_PER_DAY/DTI)
!
!  -----SUMMING UP THE FIELDS TO BE SAVED-----
!
   call calcmean_bfm(ACCUMULATE)
!
!-----WRITE OUTPUT-----
!
   if(TT+(DTI/SEC_PER_DAY).gt.out_delta) then
!
      call calcmean_bfm(MEAN)
      call save_bfm(localtime)
!
!     -----RESET TIME COUNTER-----
!
         TT = TT - nint(TT)
!
      end if
!
#ifndef BFM_POM
!
!-----WHEN RUNNING IN COUPLING WITH POM------
!-----RESTART IS WRITTEN IN MAIN AFTER THE END-----
!-----OF THE TIME MARCHING LOOP-----
!----- 
!-----WRITE RESTART-----
!
    if ( kt >= IEND ) then
!
!
         if(-TT.le.DTI/SEC_PER_DAY) then
!
             call save_rst_bfm(localtime)       
             call close_ncdf(ncid_rst)
             call close_ncdf(ncid_bfm)
             write (6,*) 'POM_DIA: NETCDF RESTART WRITTEN, TIME--> ', time,kt, iend,DTI,TT
!
        endif
!
   endif
#endif
!
   return
   end subroutine pom_dia_bfm

!EOC

