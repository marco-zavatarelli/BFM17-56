#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelForcingForBen
!
! DESCRIPTION
!
! !INTERFACE
  subroutine PelForcingForBenDynamics
!
#ifdef INCLUDE_BEN
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem, ONLY: R6c, R6n, R6p, R6s, N1p, N3n, N4n, N5s, N6r, O2o, &
    iiPhytoPlankton, ppPhytoPlankton, PhytoPlankton, D2STATE_BEN, &
    iiMicroZooPlankton, ppMicroZooPlankton, MicroZooPlankton
  use mem, ONLY: ppR6c, ppR6n, ppR6p, ppR6s, ppN1p, ppN3n, &
    ppN4n, ppN5s, ppN6r, ppO2o, &
    ETW, ESW, ERHO, ETW_Ben, ESW_Ben, ERHO_Ben, &
    PI_Benc, PI_Benn, PI_Benp, PI_Bens, sediPPY_Ben,sediR6_Ben, sediPPY, sediR6, &
    Depth, RI_Fc, ZI_Fc, ZI_Fn, ZI_Fp, RI_Fn, RI_Fp, &
    RI_Fs, N1p_Ben, N3n_Ben, N4n_Ben, N5s_Ben, N6r_Ben, O2o_Ben, ETW_Ben, &
    Depth_Ben, iiP1, iiC, iiN, iiP, iiS, iiBen, iiPel, flux
#ifdef INCLUDE_BENCO2
    use mem, ONLY: O3c_Ben,O3c,O3h_Ben,O3h
#endif
#endif
  use mem_MicroZoo, ONLY:p_qncMIZ,p_qpcMIZ
  use mem_Param,  ONLY: p_small
#ifdef BFM_GOTM
  use bio_var, ONLY: BOTindices
#else
  use api_bfm, ONLY: BOTindices
#endif

!  
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   Created at Wed Jun 16 02:04:44 PM CEST 2004
!
!
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, the BFM team
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
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
  real(RLEN), dimension(:), pointer  ::lcl_MicroZooPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i
  integer  :: j

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Compute total phytoplankton conc. used as food for filtereeders
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      do i = 1 , ( iiPhytoPlankton)
        lcl_PhytoPlankton => PhytoPlankton(i,iiC)
        PI_Benc(i,:)  =   lcl_PhytoPlankton(BOTindices)
        lcl_PhytoPlankton => PhytoPlankton(i,iiN)
        PI_Benn(i,:)  =   lcl_PhytoPlankton(BOTindices)
        lcl_PhytoPlankton => PhytoPlankton(i,iiP)
        PI_Benp(i,:)  =   lcl_PhytoPlankton(BOTindices)
        j=ppPhytoPlankton(i,iiS)
        if ( j > 0 ) then
          lcl_PhytoPlankton => PhytoPlankton(i,iiS)
          PI_Bens(i,:)  =   lcl_PhytoPlankton(BOTindices)
        else
          PI_Bens(i,:)  =   0.0
        end if
      end do
      sediPPY_Ben(:,:)  =  sediPPY(:,BOTindices)    

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Compute total microzooplankton conc. used as food for filtereeders
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ZI_Fc(:)  =   ZERO
      ZI_Fn(:)  =   ZERO
      ZI_Fp(:)  =   ZERO

      do i = 1 , ( iiMicroZooPlankton)
        lcl_MicroZooPlankton => MicroZooPlankton(i,iiC)
        ZI_Fc(:)  =   ZI_Fc(:)+ lcl_MicroZooPlankton(BOTindices)
        j = ppMicroZooPlankton(i,iiN)
        if ( j> 0) then
           lcl_MicroZooPlankton => MicroZooPlankton(i,iiN)
           ZI_Fn(:)  =   ZI_Fn(:)+ lcl_MicroZooPlankton(BOTindices)
        else
           ZI_Fn(:)  =   ZI_Fn(:)+ lcl_MicroZooPlankton(BOTindices)*p_qncMIZ(i)
        endif
        j = ppMicroZooPlankton(i,iiP)
        if ( j> 0) then
          lcl_MicroZooPlankton => MicroZooPlankton(i,iiP)
          ZI_Fp(:)  =   ZI_Fp(:)+ lcl_MicroZooPlankton(BOTindices)
        else
          ZI_Fp(:)  =   ZI_Fp(:)+ lcl_MicroZooPlankton(BOTindices)*p_qpcMIZ(i)
        endif
      enddo

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Compute total detritus conc. used as food for filtereeders
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      RI_Fc(:)  =   R6c(BOTindices)
      RI_Fn(:)  =   R6n(BOTindices)
      RI_Fp(:)  =   R6p(BOTindices)
      RI_Fs(:)  =   R6s(BOTindices)
      sediR6_Ben(:)  =   sediR6(BOTindices)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive Forcing for benthos
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !Nutrient Forcing:
      N1p_Ben(:)  =   max(p_small,N1p(BOTindices))
      N3n_Ben(:)  =   max(p_small,N3n(BOTindices))
      N4n_Ben(:)  =   max(p_small,N4n(BOTindices))
      N5s_Ben(:)  =   max(p_small,N5s(BOTindices))
      N6r_Ben(:)  =   max(p_small,N6r(BOTindices))

      !Oxygen Forcing:
      O2o_Ben(:)  =   max(p_small,O2o(BOTindices))

      ! Physical environment in the benthos is equal to the one in the
      ! adjacent level (layer) of the pelagic
      ETW_Ben(:)  =   ETW(BOTindices)
      ESW_Ben(:)  =   ESW(BOTindices)
      ERHO_Ben(:) =   ERHO(BOTindices)

#ifdef INCLUDE_BENCO2
      O3c_Ben(:)  =   O3c(BOTindices)
      O3h_Ben(:)  =   O3h(BOTindices)
#endif

      ! depth of the level above the sediment
      Depth_Ben(:)  =   Depth(BOTindices)

#endif
  end subroutine PelForcingForBenDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
