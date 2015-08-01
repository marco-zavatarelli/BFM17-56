!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: InitTransportStateTypes
!
! DESCRIPTION
!   Defining way of transport/integration of for Statevariables

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine InitTransportStateTypes
!
! USES:
  use global_mem
  use mem

!  
!
! !AUTHORS
!   P. Ruardij and M. Vichi
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
  ! Setting of type for transport/integration  Pelagic state variables
  ! All variables are transported by default
  ! D3STATETYPE(:)=ALLTRANSPORT
  ! Put D3STATETYPE(..)=NOTRANSPORT to exclude
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  D3STATETYPE(:)=ALLTRANSPORT
  do i = 1,iiMicroZooPlankton
     j = ppMicroZooPlankton(i,iiN)
     if ( i == 0 ) D3STATETYPE(i)=NOTRANSPORT
     j = ppMicroZooPlankton(i,iiP)
     if ( i == 0 ) D3STATETYPE(i)=NOTRANSPORT
  end do
  do i = 1,iiMesoZooPlankton
     j = ppMesoZooPlankton(i,iiN)
     if ( i == 0 ) D3STATETYPE(i)=NOTRANSPORT
     j = ppMesoZooPlankton(i,iiP)
     if ( i == 0 ) D3STATETYPE(i)=NOTRANSPORT
  end do
#ifdef INCLUDE_BEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Setting of type for transport/integration  Benthic state variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  D2STATETYPE_BEN(:)=NOTRANSPORT
#endif
#ifdef INCLUDE_SEAICE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Setting of type for transport/integration  Seaice state variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  D2STATETYPE_ICE(:)=NOTRANSPORT
#endif

  end subroutine InitTransportStateTypes
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
