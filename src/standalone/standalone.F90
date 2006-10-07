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
   use global_mem, only:RLEN
   use constants,  only:SEC_PER_DAY
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public envforcing_bfm,timestepping,init_standalone
   public temperature,salinity,light,lightAtTime,daylength,instLight
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
   ! forcing function parameters
   !---------------------------------------------
   real(RLEN), public :: tw,ts,tde,sw,ss,lw,ls
   real(RLEN), public :: botdep_c,botdep_n,botdep_p,botdep_si,botox_o
   !---------------------------------------------
   ! arrays for integration routines
   !---------------------------------------------
   real(RLEN),public,dimension(:,:),allocatable :: bbccc3D,bccc3D,ccc_tmp3D
   real(RLEN),public,dimension(:,:),allocatable :: bbccc2D,bccc2D,ccc_tmp2D
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
!  $Log: $
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
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_STATES,Depth,Depth_ben, D3STATE, D2STATE
   use global_mem, only:RLEN,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
   use api_bfm
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   use time

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
   namelist /anforcings_nml/ lw,ls,sw,ss,tw,ts,tde,     &
            botdep_c,botdep_n,botdep_p,botdep_si,botox_o
   namelist /time_nml/ timefmt,MaxN,start,stop,simdays
!
! !LOCAL VARIABLES:
   real(RLEN) :: tt
   integer    :: dtm1,i
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_standalone'
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
   lw          = 9.0
   ls          = 11.0
   sw          = 33.0
   ss          = 37.0
   tw          = 10.0
   ts          = 25.0
   tde         = 1.0
   botdep_c    = 0.0
   botdep_n    = 0.0
   botdep_p    = 0.0
   botdep_si   = 0.0
   botox_o     = 0.0

   open(namlst,file='standalone.nml',status='old',action='read',err=100)
   read(namlst,nml=standalone_nml,err=101)
   close(namlst)
   open(namlst,file='standalone.nml',status='old',action='read',err=100)
   read(namlst,nml=anforcings_nml,err=102)
   close(namlst)
   open(namlst,file='standalone.nml',status='old',action='read',err=100)
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
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES * NO_BOXES_XY
   LEVEL3 'Number of Boxes:',nboxes
   LEVEL3 'Box Depth:',indepth


   !---------------------------------------------
   ! initialise the timestepping parameters
   ! Use the GOTM time-manager
   ! Time is given in Julian day and seconds of day
   ! (both integer values)
   !---------------------------------------------
   timestep = maxdelt
   call init_time(MinN,MaxN)
   if (HasRealTime==.true.) then
      timesec=julianday*SEC_PER_DAY+secondsofday
      simdays=nint(simtime/SEC_PER_DAY)
   else
      timesec=0.0
   end if
   nmaxdelt=1
   LEVEL3 'nmaxdelt: ',nmaxdelt
   tt=maxdelt/2.
   do while (tt.ge.mindelt)
      tt=tt/2.
      nmaxdelt=nmaxdelt*2
      LEVEL3 'nmaxdelt: ',nmaxdelt
   end do
   mindelt=maxdelt/nmaxdelt ! maxdelt = nmaxdelt*mindelt
   nendtim=MaxN
   nstep=nmaxdelt
   ntime=0
   nmin=0
   dtm1=maxdelt
   delt=maxdelt
   if (method.eq.3) delt=2*delt
   LEVEL3 'Integration method: ',method
   LEVEL3 'maxdelt (sec): ',maxdelt
   LEVEL3 'mindelt (sec): ',mindelt
   LEVEL3 'nmaxdelt: ',nmaxdelt
   LEVEL3 'Simulation time (days): ',simdays
   LEVEL3 'nendtim: ',nendtim
   LEVEL3 'Initial time (sec): ',timesec

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
   !---------------------------------------------
   ! the argument list is kept for compatibility with GOTM
   call init_var_bfm(namlst,'bfm.nml',unit,bio_setup)
   !---------------------------------------------
   ! Assign depth
   !---------------------------------------------
   Depth = indepth
   Depth_ben = Depth
   !---------------------------------------------
   ! initialise netcdf output
   !---------------------------------------------
   call calcmean_bfm(INIT)
   call init_netcdf_bfm(out_title,'01-01-0000',0,  &
             lat=latitude,lon=longitude,z=Depth,   &
             oceanpoint=(/(i,i=1,NO_BOXES)/),      &
             surfacepoint=(/(i,i=1,NO_BOXES_XY)/), &
             bottompoint=(/(i,i=1,NO_BOXES_XY)/))
   call init_save_bfm
   !---------------------------------------------
   ! allocate and initialise integration arrays
   !---------------------------------------------
   allocate(bbccc3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(bccc3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(ccc_tmp3D(NO_D3_BOX_STATES,NO_BOXES))
   allocate(bbccc2D(NO_D2_BOX_STATES,NO_BOXES))
   allocate(bccc2D(NO_D2_BOX_STATES,NO_BOXES))
   allocate(ccc_tmp2D(NO_D2_BOX_STATES,NO_BOXES))
   ! Initialize prior time step for leap-frog:
   if (method == 3) then
      bbccc3d = D3STATE
      bbccc2d = D2STATE
      ccc_tmp3D = D3STATE
      ccc_tmp2D = D3STATE
   end if

   return

100   call error_msg_prn(NML_OPEN,"standalone.f90","standalone.nml")
101   call error_msg_prn(NML_READ,"standalone.f90","standalone_nml")
102   call error_msg_prn(NML_READ,"standalone.f90","anforcings_nml")
103   call error_msg_prn(NML_READ,"standalone.f90","time_nml")

   end subroutine init_standalone
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcings used in the BFM
!
! !INTERFACE
   subroutine envforcing_bfm()
!
! !DESCRIPTION
!
! !USES
   use api_bfm
   use global_mem, only: RLEN
   use mem,        only: ETW,ESW,EIR,ESS,SUNQ,ThereIsLight, &
                         rutQ6c,rutQ6n,rutQ6p,rutQ6s,R6c,R6n,R6p,R6s,O2o
   use mem_Param,  only: LightForcingFlag,p_PAR
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   real(RLEN) :: dfrac,wlight,dtime
   integer    :: dyear
   real(RLEN),external :: GetDelta
   real(RLEN) :: biodelta
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   LEVEL1 'envforcing_bfm'
   LEVEL2 'time=',timesec
#endif
   !---------------------------------------------
   ! Computes all the forcings
   !---------------------------------------------
   dtime = timesec/SEC_PER_DAY
   sunq=daylength(dtime,latitude)
   dfrac=(dtime-floor(dtime)) ! fraction of the day
   dyear=mod(dtime,360.) ! Day of the year
   wlight=light(dyear,dfrac)
   select case(LightForcingFlag)
    case (3) ! light on/off distribution for daylight average
      ThereIsLight=lightAtTime(dfrac,sunq)
      wlight=wlight*ThereIsLight
    case (1) ! instantaneous light distribution
      wlight=instLight(wlight,sunq,dfrac)
    case default ! light constant during the day
   end select
   ESS = 0.
   ETW = temperature(dyear,dfrac)
   ESW = salinity(dyear,dfrac)
   EIR = wlight*p_PAR
#ifdef DEBUG
   LEVEL2 'ETW=',ETW
   LEVEL2 'ESW=',ESW
   LEVEL2 'EIR=',EIR
#endif
   call CalcVerticalExtinction

   if (bio_setup==2) then
      ! Bottom deposition and ventilation fluxes
      ! (mg C m^-2 d^-1 or mmol NUT m^-2 d^-1)
      ! currently constant deposition rates read from namelist
      ! (se to zero for no deposition)
      biodelta=GetDelta()
      R6c(:) = R6c(:)+botdep_c*biodelta
      R6n(:) = R6n(:)+botdep_n*biodelta
      R6p(:) = R6p(:)+botdep_p*biodelta
      R6s(:) = R6s(:)+botdep_si*biodelta
      O2o(:) = O2o(:)+botox_o*biodelta
   end if

   end subroutine envforcing_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION daylength(time,latitude)
!
! !DESCRIPTION:
! This function computes the length of the daylight period in hours
! as a function of time of the year (days) and latitude
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in) :: time
   real(RLEN),intent(in) :: latitude
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: daylength
!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   real(RLEN)           :: declination
   real(RLEN),parameter :: cycle=360.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   declination = -0.406*cos(2.*PI*int(time)/cycle)
   daylength = acos(-tan(declination)*tan(latitude*RFACTOR))/PI*24.
   return

   END FUNCTION daylength
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION lightAtTime(df,dl)
!
! !DESCRIPTION:
!  This function determines whether there is light at a certain time
!  of the day. Returns an integer value 0 or 1
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
 real(RLEN),intent(in) :: df,dl

!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
 integer :: lightAtTime

!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
 real(RLEN) :: daytime,daylength
!
!EOP
!-----------------------------------------------------------------------
!BOC
    daytime=df*24. ! time of the day = fraction of the day * 24
    daytime=abs(daytime-12.) ! distance from noon
    daylength=dl/2.
    if(daytime.lt.daylength) then
      lightAtTime=1
    else
      lightAtTime=0
    endif
    return
   END FUNCTION lightAtTime
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION instLight(l,dl,df)
!
! !DESCRIPTION:
!  This function computes the instantaneous light at a certain time of
!  the day
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in) :: df,dl,l
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: instLight
!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
   real(RLEN) :: daylength,daytime
!
!EOP
!-----------------------------------------------------------------------
!BOC
     daytime=df*24. ! time of the day = fraction of the day * 24
     daytime=abs(daytime-12.) ! distance from noon
     daylength=dl/2.
     if(daytime.lt.daylength) then
       daytime=daytime/daylength*PI
       instLight=l*cos(daytime)+l
     else
       instLight=0.
     endif
     return
   END FUNCTION instLight
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION salinity(dy,df)
!
! !DESCRIPTION:
!  This function provides an articial salinity value given the
!  parameters in the standalone.nml namelist
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
!
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: salinity

!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:

!
!EOP
!-----------------------------------------------------------------------
!BOC
     salinity=(ss+sw)/2.-(ss-sw)/2.*cos((dy+(df-.5))*RFACTOR)
   END FUNCTION salinity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION temperature(dy,df)
!
! !DESCRIPTION:
!  This function provides an articial temperature value given the
!  parameters in the standalone.nml namelist
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
!
!
! !INPUT/OUTPUT PARAMETERS:

!
! !OUTPUT PARAMETERS:
   real(RLEN) :: temperature

!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:

!
!EOP
!-----------------------------------------------------------------------
!BOC
     temperature=(ts+tw)/2.-(ts-tw)/2.*cos((dy+(df-.5))*RFACTOR) &
                    -tde*.5*cos(2*Pi*df)
   END FUNCTION temperature
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   FUNCTION light(dy,df)
!
! !DESCRIPTION:
!  This function provides an articial light value given the
!  parameters in the standalone.nml namelist
!
! !USES:
   use global_mem, only:RLEN
   IMPLICIT NONE
! !INPUT PARAMETERS:
   integer,intent(in)    :: dy
   real(RLEN),intent(in) :: df
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   real(RLEN) :: light
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO)
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
     light=(ls+lw)/2.-(ls-lw)/2.*cos(dy*RFACTOR)
   END FUNCTION light
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
   use api_bfm, only: out_delta
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

   LEVEL1 'timestepping'

   do while (ntime.le.nendtim)
#ifdef DEBUG
      LEVEL2 'ntime=',ntime
#endif
      call envforcing_bfm
      call EcologyDynamics
      select case (method)
         case (2)
            call integrationRK2
         case (3)
            call integrationLf
         case default
            call integrationEfw
      end select
      call calcmean_bfm(ACCUMULATE)
      if (mod(ntime,out_delta).eq.0) then
         LEVEL1 'OUTPUT' , timesec/SEC_PER_DAY
         call calcmean_bfm(MEAN)
         call save_bfm(timesec)
      end if
      call ResetFluxes
   end do

   END SUBROUTINE timestepping
!EOC

!-----------------------------------------------------------------------

   END MODULE standalone

