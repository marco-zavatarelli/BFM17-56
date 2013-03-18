#include "DEBUG.h"


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcCO2SatInField.f90
!
! DESCRIPTION
!   !
!
! !INTERFACE
  SUBROUTINE CalcCO2SatInField(kmax,numc,ERHO,ETW,ESW,cc)

!
#ifdef BFM_NS
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
  use mem, ONLY: ppN1p,ppN5s,ppO3h,ppO3c
  use CO2System,ONLY: CalcCO2System,HplusBASIS 
  use mem_PelCO2    

!  
!
! !AUTHORS
!   H. Thomas   and  P. Ruardij
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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
  integer                    :: kmax
  integer                    :: numc
  real(RLEN),intent(IN)      :: ERHO(1:kmax)
  real(RLEN),intent(IN)      :: ETW(1:kmax)
  real(RLEN),intent(IN)      :: ESW(1:kmax)
  real(RLEN),intent(INOUT)   :: cc(1:numc,1:kmax)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,external    ::D3toD1
  integer             ::error
  integer             ::i
  integer             ::j
  real(RLEN)          ::Ac
  real(RLEN)          ::DIC
  real(RLEN)          ::dumCO2
  real(RLEN)          ::dumCO3
  real(RLEN)          ::dumHCO3
  real(RLEN)          ::dumpH

  if ( ppO3c > 0 ) then
    do i=1,kmax
          j=ppO3h
          Ac   =   cc(j,i)+ HplusBASIS
          error= CalcCO2System(MethodCalcCO2,ESW(i),ETW(i),ERHO(i),&
                   cc(ppN1p,i),cc(ppN5s,i),Ac,&
                   dumCO2,dumHCO3,dumCO3,dumpH,&
                   pCO2_in=pCO2_Air,DIC_out=DIC)
          j=ppO3c
          cc(j,i)=DIC* 12.0;
    enddo
  endif

#endif
  end subroutine CalcCO2SatInField
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
