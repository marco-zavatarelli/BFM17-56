#include "INCLUDE.h"
#include "cppdefs.h"
!
! !ROUTINE:  restart_BFM_inPOM
!
!DESCRIPTION    
!
!  This routine writes the file needed for the BFM restart 
!
! !INTERFACE
!
   subroutine restart_BFM_inPOM
!
! !USES:
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use netcdf_bfm, only: save_rst_bfm, ncid_rst
   use constants,  only:SEC_PER_DAY
   use pom, ONLY: time 
   use global_mem, ONLY: RLEN
   use Mem
   use api_bfm
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
!  -----TIME IN SECONDS-----
!
   real(RLEN)               :: localtime 
!
!BOC
!
! -----TIME ELLAPSED (IN SECONDS) SINCE THE BEGINNING OF THE SIMULATION----
!
   localtime = time*SEC_PER_DAY   
!
!  -----WRITE RESTART-----
!
   call save_rst_bfm(localtime)       
!
! -----CLOSE OUTPUT AND RESTART FILES------
!
   call close_ncdf(ncid_rst)
   call close_ncdf(ncid_bfm)

 write (6,*) 'NETCDF RESTART WRITTEN, TIME--> ', time!,kt, nitend   (G)
!
   return
   end subroutine restart_BFM_inPOM

!EOC

