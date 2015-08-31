!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ModuleGlobalMem
!
! DESCRIPTION
!    !   Definiton of the runtime error messages.
!   This module contains global settings:
!   -general constants for controlling prescision,
!   -parameters defining fle streams and error message numbers
!   -the subroutine for printing the message
!   and aborting the simulation
!
! !INTERFACE
  MODULE global_mem
!

!  
!
! !AUTHORS
!   P. Ruardij (NIOZ)/ M. Vichi (INGV) 
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi 
!   (rua@nioz.nl, vichi@bo.ingv.it)
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
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Global Constants
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! default REALS are double precision
  integer, parameter :: sp = selected_real_kind(6, 37)  ! real 4 (6 digits)
  integer, parameter :: dp = selected_real_kind(12,307) ! real 8 (12 digits)
  integer, parameter :: RLEN = dp                       ! default
  integer, parameter :: NMLUNIT=310
  ! the unit of the LOG file is not a parameter to allow parallel writing
  integer               ::LOGUNIT=0
  logical               ::bfm_lwp = .TRUE. ! logical writing for proc 0
  real(RLEN), parameter ::ZERO=0.0_RLEN
  real(RLEN), parameter ::ONE=1.0_RLEN
  real(RLEN), parameter ::PI=3.14159265359_RLEN
  real(RLEN), parameter ::BASETEMP= 10.0_RLEN
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Next parameters are defined to control types of State variables
  integer,    parameter ::OFF=-100
  integer,    parameter ::SINKSOURCE=-1
  integer,    parameter ::NOTRANSPORT=0
  integer,    parameter ::NOOBCSTATES=0
  integer,    parameter ::OBCSTATES=1
  integer,    parameter ::HORTRANSPORT=10
  integer,    parameter ::ALLTRANSPORT=20
  real(RLEN), parameter :: DONE=1._RLEN
  ! Error codes:
  integer,    parameter ::ALLOC=10
  integer,    parameter ::NML_OPEN=11
  integer,    parameter ::NML_READ=12
  integer,    parameter ::DIM_MISMATCH=13
  contains

  subroutine error_msg_prn(code,infile,what)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer :: code
  character(LEN=*) :: what,infile
   write(LOGUNIT,*) "*********** RUN TIME ERROR BEGIN ***********"
   select case (code)
        case (ALLOC)
          write(LOGUNIT,*) "Unable to allocate ",trim(what)," in ", &
                                trim(infile)
        case (NML_OPEN)
          write(LOGUNIT,*) "Unable to open ",trim(what)," in ",  &
                                trim(infile)
        case (NML_READ)
          write(LOGUNIT,*) "Namelist mismatch in ",trim(what),  &
                         " opened by ",trim(infile)
        case (DIM_MISMATCH)
          write(LOGUNIT,*) "Dimension mismatch while reading ",  &
                trim(what)," in ",trim(infile)
    end select
    write(LOGUNIT,*) "***********  RUN TIME ERROR END  ***********"
    stop "BFM error (see logfile)"
  end subroutine error_msg_prn
  
  end module global_mem
  
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
