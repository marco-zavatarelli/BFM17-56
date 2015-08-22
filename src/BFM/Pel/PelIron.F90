#include "DEBUG.h"
#include "INCLUDE.h"
#ifdef INCLUDE_PELFE
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelIronDynamics
!
! DESCRIPTION
!       This process describes the additional dynamics of dissolved
!       iron in the watercolumn. Parameterized processes are:
!       - remineralization of bioavailable iron
!       - scavenging of dissolved iron
!
! !INTERFACE
  subroutine PelIronDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: iiPel,NO_BOXES,ETW,flux_vector
  use mem,  ONLY: N7f,R6f,R1f,ppN7f,ppR6f,ppR1f
#endif
  use mem_Param,  ONLY: p_small
  use mem_PelChem

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used: eTq_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector

!  
!
! !AUTHORS
!   Original version by M. Vichi
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2011 BFM System Team
!   (marcello.vichi@bo.ingv.it)
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, save :: first =0
  integer       :: AllocStatus
  real(RLEN),allocatable,save,dimension(:) :: fR1N7f, fR6N7f
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if (first==0) then
     first=1
     allocate(fR6N7f(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fR6N7f"
     allocate(fR1N7f(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fR1N7f"
  end if

  !-=-==-=-=-=-=-===--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
  ! Linear regeneration of bioavailable iron
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  fR1N7f(:)  =  p_sR1N7* eTq_vector(  ETW(:),  p_q10R6N7)* R1f(:)
  call flux_vector( iiPel, ppR1f, ppN7f, fR1N7f(:) )

  fR6N7f(:)  =  p_sR6N7* eTq_vector(  ETW(:),  p_q10R6N7)* R6f(:)
  call flux_vector( iiPel, ppR6f, ppN7f, fR6N7f(:) )

  !-=-==-=-=-=-=-===--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=
  ! Scavenging to particles
  ! Linear relaxation to the solubility value p_N7fsol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  fscavN7f(:) = max(ZERO,p_scavN7f*(N7f-p_N7fsol))
  call flux_vector( iiPel, ppN7f, ppN7f, -fscavN7f(:) )

  end subroutine PelIronDynamics
#endif
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
