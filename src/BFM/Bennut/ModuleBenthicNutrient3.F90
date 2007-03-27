!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicNutrient3
!
! DESCRIPTION
!   This file holds all alternative parameter values for parameters
!	defined in the process BenthicNutrient3.p.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_BenthicNutrient3
!
! !USES:

  use global_mem
!  
!
! !AUTHORS
!   P Ruardij
!
!
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenthicNutrient3 PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer   :: bennut_messages=0
  real(RLEN)  ::p_max_shift_change=0.25
  real(RLEN)  ::p_max_state_change=0.25
 
  !  limited to half of its value
  !  this parameter inhibits too
  !  thick oxic layers.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitBenthicNutrient3
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitBenthicNutrient3

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenthicNutrient3_parameters/ bennut_messages,p_max_shift_change,p_max_state_change
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading BenthicNutrient3 parameters.."
   open(NMLUNIT,file='BenthicNutrient3.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=BenthicNutrient3_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=BenthicNutrient3_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenthicNutrient3.f90","BenthicNutrient3.nml")
101 call error_msg_prn(NML_READ,"InitBenthicNutrient3.f90","BenthicNutrient3_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitBenthicNutrient3
  end module mem_BenthicNutrient3
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
