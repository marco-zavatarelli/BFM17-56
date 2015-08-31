#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: standalone
!
! !INTERFACE:
   module standalone
!
! !DESCRIPTION:
! This module contains all the routines for the standard BFM
! simulation in standalone version, i.e. with a 0.5D setup
! (pelagic and benthic).
! It also contains all the ancillary functions for forcing functions.
!
! !USES:
!  default: all is private.
   use global_mem, only:RLEN,ONE
   use constants,  only:SEC_PER_DAY
   use init_var_bfm_local
#ifdef INCLUDE_PELCO2
   use mem_CO2, ONLY: CloseCO2
#endif   
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public timestepping,init_standalone,end_standalone

! !PUBLIC DATA MEMBERS:
   ! Note: all read from namelist
   !---------------------------------------------
   ! geographic and dimensional parameters
   !---------------------------------------------
   real(RLEN),public        :: latitude,longitude
   integer, public          :: nboxes
   real(RLEN),public        :: indepth
   !---------------------------------------------
   ! timestepping parameters
   ! real:
   ! mindelt: minimal time step allowed for computation
   ! endtim: endtime of integration
   ! time: actual time
   ! delt: actual time step of global integration
   ! maxdelt: maximal timestep
   !---------------------------------------------
   real(RLEN),public  :: maxdelt,mindelt,endtim, &
                         timesec,delt
   !---------------------------------------------
   ! integer:
   ! nmaxdelt: number of mindelts in maxdelts
   ! nendtim: number of total maxelts to endtim
   ! nmin: actual no. of mindelts in maxdelt intervall
   ! nstep: actual time step in mindelts
   ! ntime: actual time in maxdelts
   ! method: integration method
   !---------------------------------------------
   integer,public     :: nmaxdelt,nendtim,nmin,nstep,ntime, &
                         method
   !---------------------------------------------
   ! arrays for integration routines
   !---------------------------------------------
   real(RLEN),public,dimension(:,:),allocatable :: bbccc3D,bccc3D,ccc_tmp3D
#if defined INCLUDE_SEAICE
   real(RLEN),public,dimension(:,:),allocatable :: bbccc2D_ice,bccc2D_ice,ccc_tmp2D_ice
#endif
#if defined INCLUDE_BEN
   real(RLEN),public,dimension(:,:),allocatable :: bbccc2D_ben,bccc2D_ben,ccc_tmp2D_ben
#endif

   real(RLEN),public                            :: dtm1
   logical,public                               :: sspflag

!
! !PRIVATE DATA MEMBERS:
   real(RLEN),parameter :: PI=3.14159265,RFACTOR=PI/180.
   integer,parameter    :: namlst=10,unit=11

!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV), Momme Buthenschoen (UNIBO)
!
!
! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the standalone BFM
!
! !INTERFACE:
   subroutine init_standalone()
!
! !DESCRIPTION:
!  Main communication of array dimensions between
!  BFM and the standalone setup
!  Read also integration limits
!
!
! !USES:
   use constants, only: E2W
   use mem,   only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_BOXES_XY, NO_D2_BOX_DIAGNOSS, &
                  NO_D3_BOX_DIAGNOSS, NO_STATES, Depth, D3STATE
   use mem,  only: Volume, Area, Area2d
   use global_mem, only:RLEN,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
   use api_bfm
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm,&
                         init_netcdf_rst_bfm,read_rst_bfm
   use time
#if defined INCLUDE_SEAICE
   use mem, only: D2STATE_ICE, NO_D2_BOX_STATES_ICE, &
                  NO_D2_BOX_DIAGNOSS_ICE, NO_BOXES_ICE, NO_BOXES_Z_ICE, &
                  NO_STATES_ICE
#endif
#ifdef INCLUDE_BEN
   use mem, only: D2STATE_BEN, NO_D2_BOX_STATES_BEN, &
                  NO_D2_BOX_DIAGNOSS_BEN, NO_BOXES_BEN, NO_BOXES_Z_BEN, &
                  NO_STATES_BEN
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   namelist /standalone_nml/ nboxes,indepth,maxdelt,    &
            mindelt,endtim,method,latitude,longitude
   namelist /time_nml/ timefmt,MaxN,start,stop,simdays
!
! !LOCAL VARIABLES:
   real(RLEN) :: tt
   integer    :: dtm1,i
   character(LEN=8)  :: datestr
!EOP
!-----------------------------------------------------------------------
!BOC
   call Date_And_Time(datestr,timestr)
   STDERR ' '
   STDERR LINE
   STDERR '          BFM STANDALONE EXPERIMENT '
   STDERR LINE
   LEVEL1 'INITIALIZE STANDALONE'
   LEVEL1 ' '
   LEVEL1 'Simulation started on  ',datestr,' ',timestr

   !---------------------------------------------
   ! Give initial default values
   ! (overwritten with namelist)
   !---------------------------------------------
   nboxes      = 1
   indepth     = 10.0
   latitude    = 0.0
   longitude   = 0.0
   maxdelt     = 100.0
   mindelt     = 1.0
   endtim      = 360.0
   method      = 1

   open(namlst,file='Standalone.nml',status='old',action='read',err=100)
   read(namlst,nml=standalone_nml,err=101)
   close(namlst)
   open(namlst,file='Standalone.nml',status='old',action='read',err=100)
   read(namlst,nml=time_nml,err=103)
   close(namlst)

   !---------------------------------------------
   ! set the dimensions
   !---------------------------------------------
   NO_BOXES_X  = nboxes
   NO_BOXES_Y  = 1
   NO_BOXES_Z  = 1
   NO_BOXES    = NO_BOXES_X * NO_BOXES_Y * NO_BOXES_Z
   NO_BOXES_XY = NO_BOXES_X * NO_BOXES_Y
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES + NO_BOXES_XY
#ifdef INCLUDE_BEN
#  ifdef INCLUDE_BENPROFILES
   ! Dirty method to cheat the standalone model
   NO_BOXES_Z_BEN  = nboxes
   NO_BOXES_BEN = NO_BOXES_XY 
#  else
   NO_BOXES_Z_BEN  = 1
   NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
#  endif
   NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN
#endif
#ifdef INCLUDE_SEAICE
   NO_BOXES_Z_ICE  = 1
   NO_BOXES_ICE = NO_BOXES_XY * NO_BOXES_Z_ICE
   NO_STATES_ICE = NO_BOXES_ICE * NO_D2_BOX_STATES_ICE
#endif

   LEVEL2 'Number of Boxes:',nboxes
   LEVEL2 'Box Depth:',indepth
   ! set where surface and bottom boxes are 
   ! (actually all boxes in standalone mode)
   allocate(SRFindices(NO_BOXES_XY))
   allocate(BOTindices(NO_BOXES_XY))
   SRFindices = (/(i,i=1,NO_BOXES_XY)/)
   BOTindices = (/(i,i=1,NO_BOXES_XY)/)

   !---------------------------------------------
   ! initialise the timestepping parameters
   ! Use the GOTM time-manager
   ! Time is given in Julian day and seconds of day
   ! (both integer values)
   !---------------------------------------------
   timestep = maxdelt
   call init_time(MinN,MaxN)
   if (HasRealTime) then
      timesec=julianday*SEC_PER_DAY+secondsofday
      simdays=nint(simtime/SEC_PER_DAY)
   else
      timesec=ZERO
   end if
   nmaxdelt=1
   !LEVEL3 'nmaxdelt: ',nmaxdelt
   tt=maxdelt/2.
   do while (tt.ge.mindelt)
      tt=tt/2.
      nmaxdelt=nmaxdelt*2
      !LEVEL3 'nmaxdelt: ',nmaxdelt
   end do
   mindelt=maxdelt/nmaxdelt ! maxdelt = nmaxdelt*mindelt
   nendtim=MaxN
   nstep=nmaxdelt
   ntime=0
   nmin=0
   dtm1=maxdelt
   delt=maxdelt
   LEVEL2 'max delta T (sec): ',maxdelt
   LEVEL2 'min delta T (sec): ',mindelt
   LEVEL2 'nmaxdelt: ',nmaxdelt
   LEVEL2 'Simulation time (days): ',simdays
   LEVEL2 'End time: ',nendtim

   if (method.eq.3) delt=2.0_RLEN*delt
   LEVEL1 ' '
   select case (method)
      case (1)
        LEVEL1 'Integration method : Euler forward'
      case (2)
        LEVEL1 'Integration method : Runge-Kutta 2nd order'
      case (3)
        LEVEL1 'Integration method : Leap-frog'
   end select

   !---------------------------------------------
   ! Initialise the BFM with standalone settings
   !---------------------------------------------
   call init_bfm(namlst)
   !---------------------------------------------
   ! Initialise state variable names and diagnostics
   !---------------------------------------------
   call set_var_info_bfm
   !---------------------------------------------
   ! Allocate memory and give initial values
   ! to the pelagic system
   !---------------------------------------------
   call init_var_bfm(bio_setup)
   !---------------------------------------------
   ! Initialize internal constitutents of functional groups
   !---------------------------------------------
   call init_organic_constituents()
   !---------------------------------------------
   ! Set output stepping
   !---------------------------------------------
   save_delta = bfmtime%step0
   call update_save_delta(out_delta,save_delta,time_delta)

   if ( out_delta .lt. 0 .and. (timefmt .eq. 1 .or. timefmt .eq. 4) ) &
      write (LOGUNIT, *) 'WARNING: If timefmt = 1 or 4, Check the start time to save monthly values!!!'
   !---------------------------------------------
   ! Assign depth
   !---------------------------------------------
   Depth = indepth
   ! assume area is 1m^2 (make a parameter in the future for 
   ! mesocosm simulations)
   Area = ONE
   Area2d = ONE
   Volume = Depth*Area
   !---------------------------------------------
   ! Initialise external forcing functions
   !---------------------------------------------
   call init_envforcing_bfm

   !---------------------------------------------
   ! Read restart file (if flag)
   ! Overwrite previous initialization
   !---------------------------------------------
   if (bfm_init == 1) call read_rst_bfm(in_rst_fname)

   !---------------------------------------------
   ! Initialise the diagnostic variables
   !---------------------------------------------
   call CalcVerticalExtinction( )
   call CalcChlorophylla( )

   !---------------------------------------------
   ! Initialise netcdf output
   !---------------------------------------------
   call calcmean_bfm(INIT)
   call calcmean_bfm(ACCUMULATE)
   call init_netcdf_bfm(out_fname,start,out_title, &
             0,lat=latitude,lon=longitude,z=Depth, &
             oceanpoint=(/(i,i=1,NO_BOXES)/),      &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/))
   call init_save_bfm
   call calcmean_bfm(RESET)
   !---------------------------------------------
   ! Initialise netcdf restart file
   !---------------------------------------------
   call init_netcdf_rst_bfm(out_rst_fname,start,0,  &
             lat=latitude,lon=longitude,z=Depth,   &
             oceanpoint=(/(i,i=1,NO_BOXES)/),      &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/))

   !---------------------------------------------
   ! allocate and initialise integration arrays
   !---------------------------------------------
   allocate(bbccc3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(bccc3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(ccc_tmp3D(NO_D3_BOX_STATES,NO_BOXES))
#if defined INCLUDE_SEAICE
   allocate(bccc2D_ice(NO_D2_BOX_STATES_ICE,NO_BOXES_ICE))
   allocate(bbccc2D_ice(NO_D2_BOX_STATES_ICE,NO_BOXES_ICE))
   allocate(ccc_tmp2D_ice(NO_D2_BOX_STATES_ICE,NO_BOXES_ICE))
#endif
#if defined INCLUDE_BEN
   allocate(bccc2D_ben(NO_D2_BOX_STATES_BEN,NO_BOXES_BEN))
   allocate(bbccc2D_ben(NO_D2_BOX_STATES_BEN,NO_BOXES_BEN))
   allocate(ccc_tmp2D_ben(NO_D2_BOX_STATES_BEN,NO_BOXES_BEN))
#endif

   ! Initialize prior time step for leap-frog:
   if (method == 3) then
      bbccc3d = D3STATE
      ccc_tmp3D = D3STATE
#if defined INCLUDE_SEAICE
      bbccc2d_ice = D2STATE_ICE
      ccc_tmp2D_ice = D2STATE_ICE
#endif
#if defined INCLUDE_BEN
      bbccc2d_ben = D2STATE_BEN
      ccc_tmp2D_ben = D2STATE_BEN
#endif
   end if

#ifdef DEBUG
   ! print out initial values (first grid point only)
   do i=1,NO_D3_BOX_STATES
     LEVEL1 trim(var_names(stPelStateS+i-1)),D3STATE(i,1)
   end do
#if defined INCLUDE_SEAICE
   do i=1,NO_D2_BOX_STATES_ICE
     LEVEL1 trim(var_names(stIceStateS+i-1)),D2STATE_ICE(i,1)
   end do
#endif
#if defined INCLUDE_BEN
   do i=1,NO_D2_BOX_STATES_BEN
     LEVEL1 trim(var_names(stBenStateS+i-1)),D2STATE_BEN(i,1)
   end do
#endif
#endif

   return

100   call error_msg_prn(NML_OPEN,"standalone.f90","Standalone.nml")
101   call error_msg_prn(NML_READ,"standalone.f90","Standalone_nml")
103   call error_msg_prn(NML_READ,"standalone.f90","time_nml")

   end subroutine init_standalone
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   SUBROUTINE timestepping()
!
! !DESCRIPTION:
!
!
! !USES:
   use global_mem, only:RLEN
   use netcdf_bfm, only: save_bfm
   use mem
   use api_bfm, only: out_delta, save_delta , time_delta, update_save_delta
#ifdef DEBUG
   use api_bfm, only: printDEBUG
#endif
   use time
   IMPLICIT NONE
! !INPUT PARAMETERS:
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
integer    :: i
real(RLEN) :: localtime

   LEVEL1 'timestepping'
   do while (ntime.le.nendtim)
#ifdef DEBUG
      LEVEL2 'ntime=',ntime
#endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Compute environmental forcing functions
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call envforcing_bfm(ntime)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Compute extinction coefficient
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call CalcVerticalExtinction( )
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Compute reaction terms
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call EcologyDynamics
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Integrate forward in time
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifdef DEBUG
      call printDEBUG(ntime)
#endif
      select case (method)
         case (2)
            call integrationRK2
         case (3)
            call integrationLf
         case default
            call integrationEfw
      end select
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! update internal bfm time 
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      bfmtime%stepnow = ntime
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Compute means and store results
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call calcmean_bfm(ACCUMULATE)
      if ( ntime .eq. save_delta ) then
         localtime = (time_delta - real(bfmtime%step0,RLEN)) * bfmtime%timestep
         LEVEL1 'OUTPUT' , localtime / SEC_PER_DAY
         call calcmean_bfm(MEAN)
         call save_bfm(localtime)
         call update_save_delta(out_delta,save_delta,time_delta)
      end if
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Reset source-sink arrays, update time
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call ResetFluxes
      call update_time(ntime)
#ifdef DEBUG
      LEVEL2 'julian, seconds=',julianday,secondsofday
#endif
   end do

   END SUBROUTINE timestepping
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalise the standalone BFM
!
! !INTERFACE:
   subroutine end_standalone()
!
! !DESCRIPTION:
!  Terminates the standalone simulation.
!
!
! !USES:
   use time
   use netcdf_bfm, only: close_ncdf, ncid_bfm, save_rst_bfm, ncid_rst
   use api_bfm, only: time_delta
   IMPLICIT NONE
! !INPUT PARAMETERS:
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
! !REVISION HISTORY:
!  Author(s): Karsten Bolding and Hans Burchard
!
! !LOCAL VARIABLES:
   character(LEN=8)          :: datestr
   real(RLEN)       :: localtime
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! save and close the restart file
   localtime = (time_delta - real(bfmtime%step0,RLEN)) * bfmtime%timestep
   call save_rst_bfm(localtime)
   call close_ncdf(ncid_rst)
   call close_ncdf(ncid_bfm)
#ifdef INCLUDE_PELCO2
   !close systemforcings
   call CloseCO2()
#endif
   call end_envforcing_bfm
   call ClearMem
   call Date_And_Time(datestr,timestr)
   STDERR LINE
   STDERR 'BFM standalone finished on  ',datestr,' ',timestr
   end subroutine end_standalone
!EOC
!-----------------------------------------------------------------------

   END MODULE standalone

