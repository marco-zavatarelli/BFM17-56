#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CheckMassConservation
!
! DESCRIPTION
!
! !INTERFACE
  subroutine CheckMassConservationDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,ONE,LOGUNIT
  use constants, ONLY: MW_P, MW_N, MW_SI
  use mem


  use mem_MesoZoo, ONLY: p_qncMEZ,p_qpcMEZ
  use mem_MicroZoo, ONLY: p_qncMIZ,p_qpcMIZ
  use mem_Param,  ONLY: CalcBenthicFlag,p_d_tot
#ifdef INCLUDE_BEN
  use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
#endif

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
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES)  :: s
  integer         :: i,j
  real(RLEN),save :: prevsysc,prevsysn,prevsysp,prevsyss
  real(RLEN),save :: initialc,initialn,initialp,initials
  integer         :: prec
  logical,save    :: first=.TRUE.
  logical,save    :: flag=.FALSE.
  real(RLEN),parameter :: p_prec=1.e-12_RLEN
!EOP
!-------------------------------------------------------------------------!
!BOC
  totpelc(:)=ZERO
  totpelp(:)=ZERO
  totpeln(:)=ZERO
  totpels(:)=ZERO
  do i=1, iiPhytoPlankton
     s=PhytoPlankton(i,iiC)
     totpelc(:)=totpelc(:) + s
     s=PhytoPlankton(i,iiN)
     totpeln(:)=totpeln(:) + s
     s=PhytoPlankton(i,iiP)
     totpelp(:)=totpelp(:) + s
     j=ppPhytoPlankton(i,iiS)
     if ( j/=0) then
        s=PhytoPlankton(i,iiS)
        totpels(:)=totpels(:) + s
     end if
  end do
  do i=1, iiMicroZooplankton
     s=MicroZooplankton(i,iiC)
     totpelc(:) = totpelc(:) + s
     j=ppMicroZooPlankton(i,iiN)
     if ( j /= 0) then
        s = MicroZooplankton(i,iiN)
     else
        s = MicroZooplankton(i,iiC)*p_qncMIZ(i)
     end if
     totpeln(:)=totpeln(:) + s
     j=ppMicroZooPlankton(i,iiP)
     if ( j /= 0) then
        s = MicroZooplankton(i,iiP)
     else
        s = MicroZooplankton(i,iiC)*p_qpcMIZ(i)
     end if
     totpelp(:)=totpelp(:) + s
  end do
  do i=1, iiMesoZooPlankton
     s=MesoZooPlankton(i,iiC)
     totpelc(:)=totpelc(:) + s
     j=ppMesoZooPlankton(i,iiN)
     if ( j /= 0) then
        s = MesoZooPlankton(i,iiN)
     else
        s = MesoZooPlankton(i,iiC)*p_qncMEZ(i)
     end if
     totpeln(:)=totpeln(:) + s
    j=ppMesoZooPlankton(i,iiP)
     if ( j /= 0) then
        s = MesoZooPlankton(i,iiP)
     else
        s = MesoZooPlankton(i,iiC)*p_qpcMEZ(i)
     end if
     totpelp(:)=totpelp(:) + s
   end do
  do i=1, iiPelDetritus
     if ( ppPelDetritus(i,iiC)/=0) then
        s=PelDetritus(i,iiC)
        totpelc(:)=totpelc(:) + s
     end if
     if ( ppPelDetritus(i,iiN)/=0) then
        s=PelDetritus(i,iiN)
        totpeln(:)=totpeln(:) + s
     end if
     if ( ppPelDetritus(i,iiP)/=0) then
        s=PelDetritus(i,iiP)
        totpelp(:)=totpelp(:) + s
     end if
     if ( ppPelDetritus(i,iiS)/=0) then
        s=PelDetritus(i,iiS)
        totpels(:)=totpels(:) + s
     end if
  end do
  do i=1, iiPelBacteria
     if ( ppPelBacteria(i,iiC)/=0) then
        s=PelBacteria(i,iiC)
        totpelc(:)=totpelc(:) + s
     end if
     if ( ppPelBacteria(i,iiN)/=0) then
        s=PelBacteria(i,iiN)
        totpeln(:)=totpeln(:) + s
     end if
     if ( ppPelBacteria(i,iiP)/=0) then
        s=PelBacteria(i,iiP)
        totpelp(:)=totpelp(:) + s
     end if
  end do

#ifdef INCLUDE_PELCO2
  totpelc(:) = totpelc(:)+ O3c(:)
#endif
  ! Convert from default units to g and multiply for the water volume
  totpelc(:) = totpelc(:)*Volume(:)/1000.0_RLEN
  ! Convert from default units to g and multiply for the water volume
  totpeln(:) = (totpeln(:)+ ( N3n(:) + N4n(:) + O4n(:))) &
               *Volume(:)*MW_N/1000.0_RLEN
  totpelp(:) = (totpelp(:)+ N1p(:)) &
               *Volume(:)*MW_P/1000.0_RLEN
  totpels(:) = (totpels(:)+ N5s(:)) &
               *Volume(:)*MW_SI/1000.0_RLEN

  totsysc(:) = sum(totpelc(:))
  totsysn(:) = sum(totpeln(:))
  totsysp(:) = sum(totpelp(:))
  totsyss(:) = sum(totpels(:))

  ! Store and check previous value
  if (first) then
     write(LOGUNIT,*) "Initializing Mass Conservation"
     first = .FALSE.
     flag  = .FALSE.
     initialc = totsysc(1)
     initialn = totsysn(1)
     initialp = totsysp(1)
     initials = totsyss(1)
  else
     prec = precision(prevsysc)
     write(LOGUNIT,*) "---> CheckMassConservation"
     write(LOGUNIT,*) "---> defined precision digits: sp=6 dp=12; operational: ",prec
     write(LOGUNIT,"(A,2D22.15)") "---> C:",totsysc(1),prevsysc
     write(LOGUNIT,"(A,2D22.15)") "---> N:",totsysn(1),prevsysn
     if (abs(totsysn(1)/initialn-ONE)>p_prec) then
        flag = .TRUE.
        write(LOGUNIT,*) "------> Change in N larger than specified precision:",p_prec,totsysn(1)/initialn-ONE
     end if
     write(LOGUNIT,"(A,2D22.15)") "---> P:",totsysp(1),prevsysp
     if (abs(totsysp(1)/initialp-ONE)>p_prec) then
        flag = .TRUE.
        write(LOGUNIT,*) "------> Change in P larger than specified precision:",p_prec,totsysp(1)/initialp-ONE
     end if
     write(LOGUNIT,"(A,2D22.15)") "---> Si:",totsyss(1),prevsyss
     if (abs(totsyss(1)/initials-ONE)>p_prec) then
        flag = .TRUE.
        write(LOGUNIT,*) "------> Change in Si larger than specified precision:",p_prec,totsyss(1)/initials-ONE
     end if
     if (flag)  then
        call flush(LOGUNIT)
        stop "Mass conservation violation in BFM! Check log file."
     end if
  end if
  prevsysc = totsysc(1) !first element is sufficient
  prevsysn = totsysn(1)
  prevsysp = totsysp(1)
  prevsyss = totsyss(1)

#ifdef INCLUDE_BEN
  ! Mass conservation variables
  totbenc(:)  =  ( Q1c(:)+ Q6c(:))
  totbenp(:)  =  ( Q1p(:)+ Q6p(:))
  totbenn(:)  =  ( Q1n(:)+ Q6n(:))
  totbens(:)  =  ( Q6s(:))
  select case ( CalcBenthicFlag)

    case ( 0 )
      totbenc(:)  =  ZERO
      totbenp(:)  =  ZERO
      totbenn(:)  =  ZERO
      totbens(:)  =  ZERO

    case ( BENTHIC_RETURN )  ! Simple benthic return
    continue

    case ( BENTHIC_BIO )  ! Intermediate benthic return
      ! Mass conservation variables
      totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ &
        Y4p(:)+ Y5p(:)+ H1p(:)+ &
        H2p(:)+ Q1p(:)+ Q6p(:)+ &
        Q11p(:)+ K1p(:)+ K11p(:))
      totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ &
        Y4n(:)+ Y5n(:)+ H1n(:)+ &
        H2n(:)+ Q1n(:)+ Q6n(:)+ &
        Q11n(:)+ K4n(:)+ K14n(:))
      totbens(:)  =  ( Q6s(:))

    case ( BENTHIC_FULL )  ! Full benthic nutrients
      totbenp(:) = ( Y1p(:)+ Y2p(:)+ Y3p(:)+ &
        Y4p(:)+ Y5p(:)+ H1p(:)+ &
        H2p(:)+ Q1p(:)+ Q6p(:)+ &
        Q11p(:)+ K1p(:)+ K11p(:)+ K21p(:))
      totbenn(:) = ( Y1n(:)+ Y2n(:)+ Y3n(:)+ &
        Y4n(:)+ Y5n(:)+ H1n(:)+ &
        H2n(:)+ Q1n(:)+ Q6n(:)+ &
        Q11n(:)+ G4n(:)+ K3n(:)+&
         K4n(:)+ K14n(:)+ K24n(:))
      totbens(:)  =   K5s(:)+ Q6s(:)

  end select

#ifdef INCLUDE_BENCO2
  totbenc(:) = totbenc(:)+G3c(:)
#endif
  ! Convert from default units to g and multiply for the sediment volume
  totbenc(:) = totbenc(:)/1000.0_RLEN*Area2d(:)*p_d_tot
  totbenn(:) = totbenn(:)*MW_N/1000.0_RLEN*Area2d(:)*p_d_tot
  totbenp(:) = totbenp(:)*MW_P/1000.0_RLEN*Area2d(:)*p_d_tot
  totbens(:) = totbens(:)*MW_Si/1000.0_RLEN*Area2d(:)*p_d_tot

  ! Add benthic mass to the total
  totsysc(:) = totsysc(:)+sum(totbenc(:))
  totsysn(:) = totsysn(:)+sum(totbenn(:))
  totsysp(:) = totsysp(:)+sum(totbenp(:))
  totsyss(:) = totsyss(:)+sum(totbens(:))
#endif

  end subroutine CheckMassConservationDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
