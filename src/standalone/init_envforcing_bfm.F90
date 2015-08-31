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
   use global_mem, only:RLEN,LOGUNIT,NML_OPEN,NML_READ,error_msg_prn
   use api_bfm
   use time
   use envforcing

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
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
! !LOCAL VARIABLES:
   logical :: use_seaice_data

   namelist /forcings_nml/ forcing_method,ltype,lw,ls,sw,ss,tw,ts,tde, & 
            ww,ws,CO2inc,botdep_c,botdep_n,botdep_p,botdep_si,botox_o,  &
            forcing_file, &
            use_seaice_data, seaice_file, &
            use_external_data, data_file, &
            use_event_data, event_file
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
       CO2inc      = 0.0
       tde         = 1.0
       botdep_c    = 0.0
       botdep_n    = 0.0
       botdep_p    = 0.0
       botdep_si   = 0.0
       open(namlst,file='Standalone.nml',status='old',action='read', &
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
#ifdef INCLUDE_SEAICE
    if (use_seaice_data) then
       LEVEL2 'Reading sea-ice forcing data from:'
       LEVEL3 trim(seaice_file)
       open(unit_seaice,file=seaice_file,action='read',status='old',err=108)
    else
       LEVEL2 'Skipping sea-ice forcing data'
    end if
#endif
    ! Read external data (if activated)
    if (use_external_data) then
       LEVEL2 'Reading external data from:'
       LEVEL3 trim(data_file)
       open(unit_data,file=data_file,action='read',status='old',err=107)
    end if
    ! Read event data (if activated)
    if (use_event_data) then
       LEVEL2 'Reading event data from:'
       LEVEL3 trim(event_file)
        open(unit_event,file=event_file,action='read',status='old',err=109)
    end if

   return

100   call error_msg_prn(NML_OPEN,"init_envforcing_bfm.f90","standalone.nml")
102   call error_msg_prn(NML_READ,"init_envforcing_bfm.f90","anforcings_nml")
106   call error_msg_prn(NML_OPEN,"init_envforcing_bfm.f90",trim(forcing_file))
107   call error_msg_prn(NML_OPEN,"init_envforcing_bfm.f90",trim(data_file))
108   call error_msg_prn(NML_OPEN,"init_envforcing_bfm.f90",trim(seaice_file))
109   call error_msg_prn(NML_OPEN,"init_envforcing_bfm.f90",trim(event_file))

   end subroutine init_envforcing_bfm
!EOC
