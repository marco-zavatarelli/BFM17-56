#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_dia_bfm.F90
!
! !INTERFACE:
   subroutine trc_dia_bfm(kt)
!
! !DESCRIPTION:
!
! !USES:
   use oce_trc
   use global_mem, only:RLEN
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use netcdf_bfm, only: save_rst_bfm, ncid_rst
   use mem
   use api_bfm, only: out_delta

   implicit none
!
! !INPUT PARAMETERS:
   integer,intent(IN)     :: kt
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): Marcello Vichi (CMCC-INGV)
!
! !LOCAL VARIABLES:
   real(RLEN)             :: localtime !time in seconds
!
!EOP
!-----------------------------------------------------------------------
!BOC

   localtime = real(kt,RLEN)*rdt

   !---------------------------------------------
   ! Update means
   !---------------------------------------------
   call calcmean_bfm(ACCUMULATE)

   !---------------------------------------------
   ! Zoom areas
   !---------------------------------------------
   ! Equatorial Pacific: instantaneous, every 5 days (75 ts)
    ! call ncout_zoom(kt,(/30,103,56,90/),"eqp",75,"ave(only(x))")

   !---------------------------------------------
   ! Write diagnostic output
   !---------------------------------------------
   if ( MOD( kt, out_delta ) == 0 ) then
      if ( lwp ) then
         write(numout,*) 'trc_dia_bfm : write NetCDF passive tracer concentrations at ', kt, 'time-step'
         write(numout,*) '~~~~~~ '
      end if
      call calcmean_bfm(MEAN)
      call save_bfm(localtime)
   end if

   !---------------------------------------------
   ! Close the file (MAV: are we saving the last step?)
   !---------------------------------------------
   if ( kt == nitend ) then
      ! save and close the restart file
      call save_rst_bfm
      call close_ncdf(ncid_rst)
      call close_ncdf(ncid_bfm)
   endif

   !---------------------------------------------
   ! Reset the arrays
   !---------------------------------------------
   call ResetFluxes

   return
   end subroutine trc_dia_bfm

!EOC

