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
   use global_mem, only: RLEN, LOGUNIT, bfm_lwp
   use netcdf_bfm, only: save_bfm, close_ncdf, ncid_bfm
   use api_bfm,    only: out_delta, save_delta, time_delta, &
                         update_save_delta, unpad_out
   use time,       only: bfmtime
#ifdef INCLUDE_PELCO2
   use constants, ONLY:MW_C
   use mem_CO2, ONLY: CloseCO2
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: O3h, O3c, DIC, ALK, ERHO
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
!  2014     : Tomas Lovato (CMCC)
!
! !LOCAL VARIABLES:
   real(RLEN)             :: localtime  !time in seconds
!
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Update means
   !---------------------------------------------
#ifdef INCLUDE_PELCO2
   ! Update DIC and alkalinity from model units to diagnostic output
   ! after transport of model state variables O3c and O3h
   ! mg C/m3 --> umol/kg
   ! mmol eq/m3 --> umol/kg
   DIC(:) = O3c(:)/MW_C/ERHO(:)*1000.0_RLEN
   ALK(:)  = O3h(:)/ERHO(:)*1000.0_RLEN
#endif
   call calcmean_bfm(ACCUMULATE)

   !---------------------------------------------
   ! Write diagnostic output
   !---------------------------------------------
   if ( bfmtime%stepnow .eq. save_delta ) then
      if ( lwp ) then
         write(numout,*) 'trc_dia_bfm : write NetCDF passive tracer concentrations at ', kt, 'time-step'
         write(numout,*) '~~~~~~ '
      end if
      localtime = (time_delta - real(bfmtime%step0,RLEN)) * bfmtime%timestep
      call calcmean_bfm(MEAN)
      call save_bfm(localtime)
      if ( unpad_out ) then
        localtime = real((bfmtime%stepnow - bfmtime%step0),RLEN) * bfmtime%timestep
        call write_rst_bfm(localtime)
      end if
      call update_save_delta(out_delta,save_delta,time_delta)
   end if
   !---------------------------------------------
   ! Save Restart and Close the files
   !---------------------------------------------
   if ( bfmtime%stepnow == bfmtime%stepEnd ) then
      if ( .NOT. unpad_out ) then
        localtime = real((bfmtime%stepEnd - bfmtime%step0),RLEN) * bfmtime%timestep
        call write_rst_bfm(localtime)
        call close_ncdf(ncid_bfm)
      endif
     !close systemforcings
#ifdef INCLUDE_PELCO2
      call CloseCO2()
#endif
      ! clear main memory
      call ClearMem

      LEVEL1 ' '
      LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
      LEVEL1 '             EXPERIMENT FINISHED               '
      LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
      LEVEL1 ' '

   else
   !---------------------------------------------
   ! Reset the arrays for next step
   !---------------------------------------------
      call ResetFluxes
   endif
   
   return
   end subroutine trc_dia_bfm

!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Create the restart NetCDF file
!
! !INTERFACE:
  subroutine write_rst_bfm(savetime)
!
! !DESCRIPTION:
! Wrapper subroutine to store restart file of BFM variables
!
! !REVISION HISTORY:
!  Original author(s): Tomas Lovato (2014)
!
! !USES:
   use global_mem, only: RLEN, LOGUNIT, bfm_lwp
   use time,       only: bfmtime
   use api_bfm,    only: out_rst_fname, parallel_rank,  &
                         ocepoint, surfpoint, botpoint
   use netcdf_bfm, only: init_netcdf_rst_bfm, save_rst_bfm, &
                         ncid_rst, close_ncdf
   use dom_oce
!
   IMPLICIT NONE
  ! * Substitutions
#include "domzgr_substitute.h90"
!
! !INPUT PARAMETERS:
  real(RLEN),intent(in)     :: savetime
! !LOCAL VARIABLES:
  character(len=PATH_MAX)   :: thistime, thisfname
  character(LEN=4)          :: str
!EOP
!-----------------------------------------------------------------------
!BOC
  LEVEL1 'write_rst_bfm: Create restart file at timestep: ',bfmtime%stepNow
  !
  ! variable parallel_rank must have been assigned previously
  ! in the coupling with the ocean model
  write(str,'(I4.4)') parallel_rank
  !
  ! Compose restart filename
  write(thistime,'(I8.8)') bfmtime%stepNow
  thisfname=TRIM(out_rst_fname)//'_'//TRIM(thistime)//'_restart_bfm_'//str

  ! Create the restart file
   call init_netcdf_rst_bfm(thisfname,TRIM(bfmtime%datestring),0,  &
             lat2d=gphit,lon2d=glamt,z=fsdept(1,1,:),  &
             oceanpoint=ocepoint,                      &
             surfacepoint=surfpoint,                   &
             bottompoint=botpoint,                     &
             mask3d=tmask)

  ! write restart data
  call save_rst_bfm(savetime)

  ! close the file
  call close_ncdf(ncid_rst)

  LEVEL1 'write_rst_bfm: restart file creation ... DONE!'

  end subroutine write_rst_bfm
!EOC

