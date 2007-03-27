#INCLUDE "DEBUG.h"
#INCLUDE "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CheckMassConservation
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine CheckMassConservationDynamics
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): P1p, P2p, P3p, P4p, &
  ! B1p, Z3p, Z4p, Z5p, Z6p, R1p, R6p, N1p, P1n, P2n, P3n, P4n, B1n, Z3n, Z4n, &
  ! Z5n, Z6n, R1n, R6n, N3n, N4n, O4n, P1s, R6s, N5s
  ! The following Benthic-states are used (NOT in fluxes): Q1p, Q6p, Q1n, Q6n, &
  ! Q6s, Y1p, Y2p, Y3p, Y4p, Y5p, H1p, H2p, Q11p, K1p, K11p, Y1n, Y2n, Y3n, &
  ! Y4n, Y5n, H1n, H2n, Q11n, K4n, K14n, K21p, G4n, K3n, K24n, K5s
  ! The following Benthic 1-d global boxvars got a value: totbenp, totbenn, &
  ! totbens
  ! The following 0-d global parameters are used: CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#IFDEF NOPOINTERS
  use mem,  ONLY: D2STATE,D3STATE
#ELSE
  use mem, ONLY: P1p, P2p, P3p, P4p, B1p,  R1p, R6p, N1p, &
    P1n, P2n, P3n, P4n, B1n, R1n, R6n, N3n, N4n, O4n, P1s, &
    R6s, N5s, D2STATE
  use mem, ONLY: Q1p, Q6p, Q1n, Q6n, Q6s, Y1p, Y2p, Y3p, Y4p, Y5p, H1p, H2p, &
    Q11p, K1p, K11p, Y1n, Y2n, Y3n, Y4n, Y5n, H1n, H2n, Q11n, K4n, K14n, K21p, &
    G4n, K3n, K24n, K5s, D3STATE
#ENDIF
  use mem, ONLY: ppP1p, ppP2p, ppP3p, ppP4p, ppB1p, ppZ3p, ppZ4p, ppZ5p, ppZ6p, &
    ppR1p, ppR6p, ppN1p, ppP1n, ppP2n, ppP3n, ppP4n, ppB1n, ppZ3n, ppZ4n, ppZ5n, &
    ppZ6n, ppR1n, ppR6n, ppN3n, ppN4n, ppO4n, ppP1s, ppR6s, ppN5s, D2STATE, Depth
  use mem, ONLY: ppQ1p, ppQ6p, ppQ1n, ppQ6n, ppQ6s, ppY1p, ppY2p, &
    ppY3p, ppY4p, ppY5p, ppH1p, ppH2p, ppQ11p, ppK1p, ppK11p, ppY1n, ppY2n, &
    ppY3n, ppY4n, ppY5n, ppH1n, ppH2n, ppQ11n, ppK4n, ppK14n, ppK21p, &
    ppG4n, ppK3n, ppK24n, ppK5s, totbenp, totbenn, totbens, NO_BOXES_XY, iiBen, &
    totpelp, totpeln, totpels, totsysp, totsysn, totsyss, &
    iiPel, flux_vector,ppMicroZooplankton,ppMesoZooPlankton,MicroZooplankton,MesoZooPlankton, &
    iiMicroZooplankton,iiMesoZooPlankton,NO_BOXES,iiN,iiP
  use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
  use mem_Param,  ONLY: CalcBenthicFlag
  use mem_MesoZoo, ONLY: p_qnMc=>p_qnc,p_qpMc=>p_qpc
  use mem_MicroZoo, ONLY: p_qn_mz,p_qp_mz

!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Mon Nov 21 09:44:23 CET 2005
!
!
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)  :: s
  integer                         ::i,j
  
  totpelp=0.0D+00;
  totpeln=0.0D+00;
  do i=1, iiMicroZooplankton
     j=max(1,ppMicroZooPlankton(i,iiN))
     s=MicroZooplankton(i,j)
     if ( j==1) s=s*p_qn_mz(i)
     totpeln=totpeln +sum(s*Depth)
     j=max(1,ppMicroZooPlankton(i,iiP))
     s=MicroZooplankton(i,j)
     if ( j==1) s=s*p_qp_mz(i)
     totpelp=totpelp +sum(s*Depth)
  enddo
  do i=1, iiMesoZooplankton
     j=max(1,ppMesoZooPlankton(i,iiN))
     s=MesoZooplankton(i,j)
     if ( j==1) s=s*p_qnMc(i)
     totpeln=totpeln +sum(s*Depth)
     j=max(1,ppMesoZooPlankton(i,iiP))
     s=MesoZooplankton(i,j)
     if ( j==1) s=s*p_qpMc(i)
     totpelp=totpelp +sum(s*Depth)
   enddo

  totpelp(1)=totpelp(1)+ sum((P1p+ P2p+ P3p+ P4p+ B1p+ R1p+ R6p+ N1p)* Depth)
  totpeln(1)=totpeln(1)+ sum((P1n+ P2n+ P3n+ P4n+ B1n+ R1n+ R6n+ N3n+ N4n+ O4n)* Depth)
  totpels(1)  = sum( (P1s+ R6s+ N5s)* Depth)


  select case ( CalcBenthicFlag)

    case ( 0 )
      totbenp(:)  =   0.0D+00
      totbenn(:)  =   0.0D+00
      totbens(:)  =   0.0D+00

    case ( BENTHIC_RETURN )  ! Simple benthic return
      ! Mass conservation variables
      totbenp(:)  =  ( Q1p(:)+ Q6p(:))
      totbenn(:)  =  ( Q1n(:)+ Q6n(:))
      totbens(:)  =  ( Q6s(:))

    case ( BENTHIC_BIO )  ! Intermediate benthic return
      ! Mass conservation variables
      totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ Y4p(:)+ Y5p(:)+ H1p(:)+ &
        H2p(:)+ Q1p(:)+ Q6p(:)+ Q11p(:)+ K1p(:)+ K11p(:))
      totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ Y4n(:)+ Y5n(:)+ H1n(:)+ &
        H2n(:)+ Q1n(:)+ Q6n(:)+ Q11n(:)+ K4n(:)+ K14n(:))
      totbens(:)  =  ( Q6s(:))

    case ( BENTHIC_FULL )  ! Full benthic nutrients
      totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ Y4p(:)+ Y5p(:)+ H1p(:)+ &
        H2p(:)+ Q1p(:)+ Q6p(:)+ Q11p(:)+ K1p(:)+ K11p(:)+ K21p(:))
      totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ Y4n(:)+ Y5n(:)+ H1n(:)+ &
        H2n(:)+ Q1n(:)+ Q6n(:)+ Q11n(:)+ G4n(:)+ K3n(:)+ K4n(:)+ K14n(:)+ K24n(:))
      totbens(:)  =   K5s(:)+ Q6s(:)

  end select

  totsysn(:)=totpeln(:)+totbenn(:)
  totsysp(:)=totpelp(:)+totbenp(:)
  totsyss(:)=totpels(:)+totbens(:)

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
