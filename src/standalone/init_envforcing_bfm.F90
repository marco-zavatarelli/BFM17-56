#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the standalone BFM
!
! !INTERFACE:
   subroutine init_envforcing_bfm
!
! !DESCRIPTION:
!  Read the namelist parameters for the analytical forcings
!  Also initialize additional components for other forcing methods
!
!
! !USES:
   use constants, only: E2W
   use mem,   only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_STATES,Depth, D3STATE, D2STATE
   use mem,  only: Volume, Area, Area2d
   use global_mem, only:RLEN,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
   use api_bfm
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   use time
   use envforcing
#ifdef INCLUDE_BEN
   use mem, only: Depth_ben
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   namelist /forcings_nml/ forcing_method,lw,ls,sw,ss,tw,ts,tde, & 
            ww,ws,botdep_c,botdep_n,botdep_p,botdep_si,botox_o,  &
            forcing_file, seaice_file
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_envforcing_bfm'
       !---------------------------------------------
       ! Give initial default values
       ! (overwritten with namelist)
       !---------------------------------------------
       lw          = 9.0
       ls          = 11.0
       sw          = 33.0
       ss          = 37.0
       tw          = 10.0
       ts          = 25.0
       ww          = 2.0
       ws          = 0.5
       tde         = 1.0
       botdep_c    = 0.0
       botdep_n    = 0.0
       botdep_p    = 0.0
       botdep_si   = 0.0
       open(namlst,file='standalone.nml',status='old',action='read', &
            err=100)
       read(namlst,nml=forcings_nml,err=102)
       close(namlst)
    select case (forcing_method)
    case (1) ! analytical forcings
    case (2) ! input data
       LEVEL2 'Reading forcing data from:'
       LEVEL3 trim(forcing_file)
       open(unit_forcing,file=forcing_file,action='read',status='old',err=106)
    case (3) ! interactive air-sea fluxes
      !call init_air_sea(data_file,latitude, longitude)
    end select

   return

100   call error_msg_prn(NML_OPEN,"standalone.f90","standalone.nml")
102   call error_msg_prn(NML_READ,"standalone.f90","anforcings_nml")
106   call error_msg_prn(NML_OPEN,"standalone.f90",trim(forcing_file))

   end subroutine init_envforcing_bfm
!EOC
