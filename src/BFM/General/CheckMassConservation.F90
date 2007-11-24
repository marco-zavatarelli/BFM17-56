#include "DEBUG.h"
#include "INCLUDE.h"

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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO
  use constants, ONLY: MW_P, MW_N, MW_SI
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE,D3STATE
#else
  use mem, ONLY: B1c, B1p, B1n, N1p, N3n, N4n, O4n, &
    N5s, O3c, D2STATE
#ifdef BFM_BENTHIC
  use mem, ONLY: Q1c, Q6c, Q1p, Q6p, Q1n, Q6n, Q6s, Y1p, Y2p, Y3p, Y4p, Y5p, H1p, H2p, &
    Q11p, K1p, K11p, Y1n, Y2n, Y3n, Y4n, Y5n, H1n, H2n, Q11n, K4n, K14n, K21p, &
    G4n, K3n, K24n, K5s, D3STATE
#endif
#endif
  use mem, ONLY: ppB1c,ppB1p,ppB1n, &
     ppN1p, ppN3n, ppN4n, ppO4n, ppN5s, ppO3c, D2STATE, &
    Depth, Volume, Area, Area2d
  use mem, ONLY: &
    totpelc, totpelp, totpeln, totpels, totsysc, totsysp, totsysn, totsyss, &
    iiPel, flux_vector,ppMicroZooplankton,ppMesoZooPlankton,MicroZooplankton,MesoZooPlankton, &
    iiMicroZooplankton,iiMesoZooPlankton,NO_BOXES,iiC,iiN,iiP,iiS,&
    PhytoPlankton,iiPhytoPlankton,ppPhytoPlankton,PelDetritus,iiPelDetritus,ppPelDetritus
  use mem_MesoZoo, ONLY: p_qnMc=>p_qnc,p_qpMc=>p_qpc
  use mem_MicroZoo, ONLY: p_qn_mz,p_qp_mz
#ifdef BFM_BENTHIC
  use mem, ONLY: ppQ1c, ppQ6c, ppQ1p, ppQ6p, ppQ1n, ppQ6n, ppQ6s, ppY1p, ppY2p, &
    ppY3p, ppY4p, ppY5p, ppH1p, ppH2p, ppQ11p, ppK1p, ppK11p, ppY1n, ppY2n, &
    ppY3n, ppY4n, ppY5n, ppH1n, ppH2n, ppQ11n, ppK4n, ppK14n, ppK21p, &
    ppG4n, ppK3n, ppK24n, ppK5s, totbenc, totbenp, totbenn, totbens, NO_BOXES_XY, iiBen
  use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
#endif
  use mem_Param,  ONLY: CalcBenthicFlag,p_d_tot

!  
!
! !AUTHORS
!   Piet Ruardij and Marcello Vichi
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)  :: s
  integer                         ::i,j
  
  totpelc=ZERO
  totpelp=ZERO
  totpeln=ZERO
  totpels=ZERO
  do i=1, iiPhytoPlankton
     s=PhytoPlankton(i,iiC)
     totpelc=totpelc + s
     s=PhytoPlankton(i,iiN)
     totpeln=totpeln + s
     s=PhytoPlankton(i,iiP)
     totpelp=totpelp + s
     j=ppPhytoPlankton(i,iiS)
     if ( j/=0) then
        s=PhytoPlankton(i,iiS)
        totpels=totpels + s
     end if
  end do
  do i=1, iiMicroZooplankton
     j=max(1,ppMicroZooPlankton(i,iiN))
     s=MicroZooplankton(i,j)
     if ( j==1) s=s*p_qn_mz(i)
     totpeln=totpeln + s
     j=max(1,ppMicroZooPlankton(i,iiP))
     s=MicroZooplankton(i,j)
     if ( j==1) s=s*p_qp_mz(i)
     totpelp=totpelp + s
  end do
  do i=1, iiMesoZooplankton
     j=max(1,ppMesoZooPlankton(i,iiN))
     s=MesoZooplankton(i,j)
     if ( j==1) s=s*p_qnMc(i)
     totpeln=totpeln + s
     j=max(1,ppMesoZooPlankton(i,iiP))
     s=MesoZooplankton(i,j)
     if ( j==1) s=s*p_qpMc(i)
     totpelp=totpelp + s
   end do
  do i=1, iiPelDetritus
     if ( ppPelDetritus(i,iiC)/=0) then
        s=PelDetritus(i,iiC)
        totpelc=totpelc + s
     end if
     if ( ppPelDetritus(i,iiN)/=0) then
        s=PelDetritus(i,iiN)
        totpeln=totpeln + s
     end if
     if ( ppPelDetritus(i,iiP)/=0) then
        s=PelDetritus(i,iiP)
        totpelp=totpelp + s
     end if
     if ( ppPelDetritus(i,iiS)/=0) then
        s=PelDetritus(i,iiS)
        totpels=totpels + s
     end if
  end do

  ! Convert from default units to g and multiply for the water volume
  totpelc = (totpelc+ ( O3c + B1c ))*Volume/1000.0_RLEN
  totpeln = (totpeln+ ( B1n + N3n + N4n + O4n))*Volume*MW_N/1000.0_RLEN
  totpelp = (totpelp+ ( B1p + N1p))*Volume*MW_P/1000.0_RLEN
  totpels = (totpels+ ( N5s ))*Volume*MW_SI/1000.0_RLEN

  totsysc = sum(totpelc(:))
  totsysn = sum(totpeln(:))
  totsysp = sum(totpelp(:))
  totsyss = sum(totpels(:))

#ifdef BFM_BENTHIC
  select case ( CalcBenthicFlag)

    case ( 0 )
      totbenc(:)  =  ZERO
      totbenp(:)  =  ZERO
      totbenn(:)  =  ZERO
      totbens(:)  =  ZERO

    case ( BENTHIC_RETURN )  ! Simple benthic return
      ! Mass conservation variables
      totbenc(:)  =  ( Q1c(:)+ Q6c(:))
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

  ! Convert from default units to g and multiply for the sediment volume
  totbenc = totbenc(:)/1000.0_RLEN*Area2d*p_d_tot
  totbenn = totbenn(:)*MW_N/1000.0_RLEN*Area2d*p_d_tot
  totbenp = totbenp(:)*MW_P/1000.0_RLEN*Area2d*p_d_tot
  totbens = totbens(:)*MW_Si/1000.0_RLEN*Area2d*p_d_tot

  ! Add benthic mass to the total
  totsysc = totsysc+sum(totbenc(:))
  totsysn = totsysn+sum(totbenn(:))
  totsysp = totsysp+sum(totbenp(:))
  totsyss = totsyss+sum(totbens(:))
#endif

  end
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
