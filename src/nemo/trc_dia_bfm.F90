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
   use global_mem, only: RLEN, LOGUNIT
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use netcdf_bfm, only: save_rst_bfm, ncid_rst
   use mem
   use api_bfm,    only: out_delta, ave_ctl
   use time,       only: bfmtime
#ifdef INCLUDE_PELCO2
   use constants, ONLY:MW_C
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: O3h, O3c, DIC, Ac, ERHO
#endif
#endif

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
!  2012     : Tomas Lovato (CMCC)
!
! !LOCAL VARIABLES:
   real(RLEN)             :: localtime  !time in seconds
!
!EOP
!-----------------------------------------------------------------------
!BOC

   localtime = real((bfmtime%stepnow - bfmtime%step0),RLEN) * rdt
   !---------------------------------------------
   ! Update means
   !---------------------------------------------
#ifdef INCLUDE_PELCO2
   ! Update DIC and alkalinity from model units to diagnostic output
   ! after transport of model state variables O3c and O3h
   ! mg C/m3 --> umol/kg
   ! mmol eq/m3 --> umol/kg
   DIC(:) = O3c(:)/MW_C/ERHO(:)*1000.0_RLEN
   Ac(:)  = O3h(:)/ERHO(:)*1000.0_RLEN
#endif
   call calcmean_bfm(ACCUMULATE)

   !---------------------------------------------
   ! Write diagnostic output
   !---------------------------------------------
   if ( MOD( (bfmtime%stepnow - bfmtime%step0), out_delta ) == 0 ) then
      if ( lwp ) then
         write(numout,*) 'trc_dia_bfm : write NetCDF passive tracer concentrations at ', kt, 'time-step'
         write(numout,*) '~~~~~~ '
      end if
      if (ave_ctl) localtime = localtime - ( real(out_delta,RLEN) * rdt / 2.0)
      call calcmean_bfm(MEAN)
      call save_bfm(localtime)
   end if
   !---------------------------------------------
   ! Save Restart and Close the files
   !---------------------------------------------
   if ( bfmtime%stepnow == bfmtime%stepEnd ) then
      call save_rst_bfm
      call close_ncdf(ncid_rst)
      call close_ncdf(ncid_bfm)
      ! clear main memory
      call ClearMem
   else
   !---------------------------------------------
   ! Reset the arrays for next step
   !---------------------------------------------
      call ResetFluxes
   endif
   
   return
   end subroutine trc_dia_bfm

!EOC

