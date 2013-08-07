#include "DEBUG.h"
#include "cppdefs.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: SeaicetoPelCoup
!
! DESCRIPTION
!
! !INTERFACE
  subroutine SeaicetoPelCoup
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   use global_mem, only: RLEN,ZERO,ONE
   use mem_Param,  only: p_PAR,CalcSeaiceAlgae,CalcPhytoPlankton
   use mem_Param,  only: CalcSeaiceZoo,CalcSeaiceBacteria 
   use mem_Param,  only: p_eps0,p_epsR6,p_epsChla, p_small
   use constants,  only: E2W, SEC_PER_DAY
   ! seaice forcings
   use mem,        only: EICE,EVB,ETB,ESB,EIB,EHB,ESI,EDH,EDS,F3c, &
                         F2o,I1p,I3n,I4n,I5s,S1l,S2l,U6c
   use mem,        only: iiIce,ppI1p,N1p,N3n,N4n,N5s,O2o,O3c,P1l,P2l,NO_BOXES_XY, &
                         flux_vector, ppI3n,ppI4n,ppI5s,ppF2o,ppF3c,&
                         ppS1l,ppS2l
   use mem,        only: SeaiceAlgae,ppSeaiceAlgae,PhytoPlankton,ppPhytoPlankton,PELSURFACE, &
                         iiS1,iiS2,iiP1,iiP2,iiSeaiceAlgae,iiPhytoPlankton,iiPel
   use mem,        only: SeaiceDetritus,ppSeaiceDetritus,PelDetritus,ppPelDetritus, &
                         iiU1,iiU6,iiR1,iiR6
   use mem,        only: SeaiceBacteria,ppSeaiceBacteria,PelBacteria,ppPelBacteria, &
                         iiT1,iiB1
   use mem,        only: SeaiceZoo,ppSeaiceZoo,MicroZooPlankton,ppMicroZooPlankton, &
                         iiX1,iiZ5
   use mem,        only: iiC,iiN,iiP,iiS,iiL
   use mem,        only: jsurN1p,jsurN3n,jsurN4n,jsurN5s,jsurO3c, &
                         jsurO2o,jsurP1l,jsurP2l,Depth,NO_BOXES_ICE,ppN1p,ppN3n,ppN4n,ppN5s,ppO2o
   use mem,        only: ESW,EIR
   use api_bfm,    only: SRFindices

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   IMPLICIT NONE
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_N1,flux_pel_ice_N3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_N4,flux_pel_ice_N5
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_O2,flux_pel_ice_O3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice, flux_atm_N1, flux_atm_N3
   real(RLEN), dimension(NO_BOXES_XY) :: I1p_tilde,I3n_tilde,I4n_tilde,I5s_tilde
   real(RLEN), dimension(NO_BOXES_XY) :: F2o_tilde,F3c_tilde
   real(RLEN), dimension(NO_BOXES_XY) :: I1p_star,I3n_star,I4n_star,I5s_star
   real(RLEN), dimension(NO_BOXES_XY) :: F2o_star,F3c_star

   integer                            :: i,j,p
   real(RLEN), dimension(:), pointer  :: lcl_PelagicVar,lcl_SeaiceVar
   real(RLEN)                         :: tmpflux(NO_BOXES_ICE)
   real(RLEN)                         :: localdelta
   real(RLEN), external               :: GetDelta
!
! !AUTHORS
! Marcello Vichi and Letizia Tedesco
!
! !REVISION_HISTORY

! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2007 the BFM team
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
    localdelta = GetDelta()

    ! set to zero all the fluxes in case there is no ice
    flux_pel_ice_N1(:)= ZERO
    flux_pel_ice_N3(:)= ZERO
    flux_pel_ice_N4(:)= ZERO
    flux_pel_ice_N5(:)= ZERO
    flux_pel_ice_O2(:)= ZERO
    flux_pel_ice_O3(:)= ZERO
    tmpflux(:)= ZERO

    ! prescription of atmospheric nutrient fluxes
    flux_atm_N1(:)=ZERO
    flux_atm_N3(:)=ZERO


    where (EHB(:) > ZERO)
       I1p_tilde(:) = I1p(:)/EHB(:)
       I3n_tilde(:) = I3n(:)/EHB(:)
       I4n_tilde(:) = I4n(:)/EHB(:)
       I5s_tilde(:) = I5s(:)/EHB(:)
       F2o_tilde(:) = F2o(:)/EHB(:)
       F3c_tilde(:) = F3c(:)/EHB(:)
   
       I1p_star(:) = max(ZERO, N1p(SRFindices))*EHB(:)
       I3n_star(:) = max(ZERO, N3n(SRFindices))*EHB(:)
       I4n_star(:) = max(ZERO, N4n(SRFindices))*EHB(:)
       I5s_star(:) = max(ZERO, N5s(SRFindices))*EHB(:)
       ! the oxygen and DIC flux do not consider the repartition coefficient otherwise concentration is too low
       F2o_star(:) = max(ZERO, O2o(SRFindices))*EHB(:)
       F3c_star(:) = max(ZERO, O3c(SRFindices))*EHB(:)

       ! flux is divided in positive (to-ice) and negative (to-water)
       flux_pel_ice_N1(:)= max(ZERO,-EDH(:)*min(ZERO,I1p(:) - I1p_star(:))) &
                           + min(ZERO, (EDH(:)*I1p_tilde(:)))
       flux_pel_ice_N3(:)= max(ZERO,-EDH(:)*min(ZERO,I3n(:) - I3n_star(:))) &
                           + min(ZERO, (EDH(:)*I3n_tilde(:)))
        flux_pel_ice_N4(:)= max(ZERO,-EDH(:)*min(ZERO,I4n(:) - I4n_star(:))) &
                           + min(ZERO, (EDH(:)*I4n_tilde(:)))
       flux_pel_ice_N5(:)= max(ZERO,-EDH(:)*min(ZERO,I5s(:) - I5s_star(:))) &
                           + min(ZERO, (EDH(:)*I5s_tilde(:)))
       flux_pel_ice_O2(:)= max(ZERO,-EDH(:)*min(ZERO,F2o(:) - F2o_star(:))) &
                           + min(ZERO, (EDH(:)*F2o_tilde(:)))
       flux_pel_ice_O3(:)= max(ZERO,-EDH(:)*min(ZERO,F3c(:) - F3c_star(:))) &
                           + min(ZERO, (EDH(:)*F3c_tilde(:)))

       ! compute irradiance at the bottom of sea ice accounting for seaice algae
       ! it is assumed that irradiance in the BAL is located at the middle
       EIR(:) = EIR(:)*(ONE-EICE(:)) +  &
                (EIB(:)*EICE)*exp(-(1.5_RLEN+p_epsR6*U6c(:)+p_epsChla*(S1l(:)+S2l(:)))*0.5_RLEN*EHB(:))
    elsewhere
       ! additional flux to the water when there is no ice to ensure complete emptyness
       flux_pel_ice_N1(:) = min(ZERO,p_small-I1p(:))/localdelta
       flux_pel_ice_N3(:) = min(ZERO,p_small-I3n(:))/localdelta
       flux_pel_ice_N4(:) = min(ZERO,p_small-I4n(:))/localdelta
       flux_pel_ice_N5(:) = min(ZERO,p_small-I5s(:))/localdelta
       flux_pel_ice_O2(:) = min(ZERO,p_small-F2o(:))/localdelta
       flux_pel_ice_O3(:) = min(ZERO,p_small-F3c(:))/localdelta
       ! first release everything then force the value to be p_small
       
       I1p(:)=max(p_small,I1p)
       I3n(:)=max(p_small,I3n)
       I4n(:)=max(p_small,I4n)
       I5s(:)=max(p_small,I5s)
       F2o(:)=max(p_small,F2o)
       F3c(:)=max(p_small,F3c)

    end where

    call flux_vector( iiIce, ppI1p,ppI1p, flux_pel_ice_N1 + flux_atm_N1) 
    call flux_vector( iiIce, ppI3n,ppI3n, flux_pel_ice_N3 + flux_atm_N3) 
    call flux_vector( iiIce, ppI4n,ppI4n, flux_pel_ice_N4 ) 
    call flux_vector( iiIce, ppI5s,ppI5s, flux_pel_ice_N5 ) 
    call flux_vector( iiIce, ppF2o,ppF2o, flux_pel_ice_O2 ) 
    call flux_vector( iiIce, ppF3c,ppF3c, flux_pel_ice_O3 ) 

    ! assign fluxes to boundary variables 
    ! Compute rates for the pelagic variables and assign to D3SOURCES/SINKS
    jsurN1p(:)=jsurN1p(:) - flux_pel_ice_N1(:)
    tmpflux(SRFindices) = jsurN1p(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN1p,ppN1p, tmpflux(:) )

    jsurN3n(:)=jsurN3n(:) - flux_pel_ice_N3(:)
    tmpflux(SRFindices) = jsurN3n(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN3n,ppN3n, tmpflux(:) )

    jsurN4n(:)=jsurN4n(:) - flux_pel_ice_N4(:)
    tmpflux(SRFindices) = jsurN4n(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN4n,ppN4n, tmpflux(:) )

    jsurN5s(:)=jsurN5s(:) - flux_pel_ice_N5(:)
    tmpflux(SRFindices) = jsurN5s(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN5s,ppN5s, tmpflux(:) )

    ! The asignement of oxygen and CO2 fluxes is done in the pelagic routines
    jsurO2o(:)=jsurO2o(:) - flux_pel_ice_O2(:)
    jsurO3c(:)=jsurO3c(:) - flux_pel_ice_O3(:)

    ! Boundary conditions for the adapted SI algae
    if (CalcSeaiceAlgae(iiS1)) then
       ! loop over iiC, iiN, iiP, iiL,iiS
       do i = 1, iiS
           lcl_SeaiceVar => SeaiceAlgae(iiS1,i)
           lcl_PelagicVar => PhytoPlankton(iiP1,i)
           where (EHB(:)>ZERO)
              flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices))
           elsewhere
              flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/localdelta
           end where
           ! add the flux and assign it to the boundary variable of the pelagic system
            j = ppSeaiceAlgae(iiS1,i)
            call flux_vector( iiIce, j,j, flux_pel_ice(:) ) 
            j = ppPhytoPlankton(iiP1,i)
            PELSURFACE(j,:) =  PELSURFACE(j,:) - flux_pel_ice(:)  
            ! map the flux into a 3D temporary array
            tmpflux(SRFindices) = PELSURFACE(j,:) / Depth(SRFindices)
            call flux_vector(iiPel, j, j, tmpflux(:) )
       end do
    end if

    ! Boundary conditions for the surviving SI algae
    if (CalcSeaiceAlgae(iiS2)) then
       ! loop over iiC, iiN, iiP, iiL,iiS
       do i = 1, iiL
           lcl_SeaiceVar => SeaiceAlgae(iiS2,i)
           lcl_PelagicVar => PhytoPlankton(iiP2,i)
           where (EHB(:)>ZERO)
              flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices))  &
               + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
           elsewhere
              flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/localdelta
           end where
           ! add the flux and assign it to the boundary variable of the pelagic system
            j = ppSeaiceAlgae(iiS2,i)
            call flux_vector( iiIce, j,j, flux_pel_ice(:) )
            j = ppPhytoPlankton(iiP2,i)
            PELSURFACE(j,:) =  PELSURFACE(j,:) - flux_pel_ice(:)
            ! map the flux into a 3D temporary array
            tmpflux(SRFindices) = PELSURFACE(j,:) / Depth(SRFindices)
            call flux_vector(iiPel, j, j, tmpflux(:) )
       end do
    end if

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Detritus Fluxes to Pelagic from Seaice
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    do i = 1, iiS
       lcl_SeaiceVar => SeaiceDetritus(iiU6,i)
       lcl_PelagicVar => PelDetritus(iiR6,i)
       where (EHB(:)>ZERO)
          flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices))  &
               + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
       elsewhere
          flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/localdelta
       end where
       ! add the flux and assign it to the boundary variable of the pelagic system
       j = ppSeaiceDetritus(iiU6,i)
       call flux_vector( iiIce, j,j, flux_pel_ice(:) ) 
       j = ppPelDetritus(iiR6,i)
       PELSURFACE(j,:) =  PELSURFACE(j,:) - flux_pel_ice(:)  
       ! map the flux into a 3D temporary array
       tmpflux(SRFindices) = PELSURFACE(j,:) / Depth(SRFindices)
       call flux_vector(iiPel, j, j, tmpflux(:) )
    end do

    do i = 1, iiP
       lcl_SeaiceVar => SeaiceDetritus(iiU1,i)
       lcl_PelagicVar => PelDetritus(iiR1,i)
       where (EHB(:)>ZERO)
          flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices))  &
               + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
       elsewhere
          flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/localdelta
       end where
       ! add the flux and assign it to the boundary variable of the pelagic system
       j = ppSeaiceDetritus(iiU1,i)
       call flux_vector( iiIce, j,j, flux_pel_ice(:) )
       j = ppPelDetritus(iiR1,i)
       PELSURFACE(j,:) =  PELSURFACE(j,:) - flux_pel_ice(:)
       ! map the flux into a 3D temporary array
       tmpflux(SRFindices) = PELSURFACE(j,:) / Depth(SRFindices)
       call flux_vector(iiPel, j, j, tmpflux(:) )
    end do
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !Bacteria Flux to Pelagic from Seaice
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   if (CalcSeaiceBacteria(iiT1)) then
       ! loop over iiC, iiN, iiP
       do i = 1, iiP
           lcl_SeaiceVar => SeaiceBacteria(iiT1,i)
           lcl_PelagicVar => PelBacteria(iiB1,i)
           where (EHB(:)>ZERO)
              flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices))  &
               + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
           elsewhere
              flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/localdelta
           end where
           ! add the flux and assign it to the boundary variable of the pelagic system
            j = ppSeaiceBacteria(iiT1,i)
            call flux_vector( iiIce, j,j, flux_pel_ice(:) ) 
            j = ppPelBacteria(iiB1,i)
            PELSURFACE(j,:) =  PELSURFACE(j,:) - flux_pel_ice(:)  
            ! map the flux into a 3D temporary array
            tmpflux(SRFindices) = PELSURFACE(j,:) / Depth(SRFindices)
            call flux_vector(iiPel, j, j, tmpflux(:) )
       end do
    end if

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !Microzooplankton Flux to Pelagic from Seaice
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   if (CalcSeaiceZoo(iiX1)) then
       ! loop over iiC, iiN, iiP
       do i = 1, iiP
           lcl_SeaiceVar => SeaiceZoo(iiX1,i)
           lcl_PelagicVar => MicroZooPlankton(iiZ5,i)
           where (EHB(:)>ZERO)
              flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices))  &
               + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
           elsewhere
              flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/localdelta
           end where
           ! add the flux and assign it to the boundary variable of the pelagic system
            j = ppSeaiceZoo(iiX1,i)
            call flux_vector( iiIce, j,j, flux_pel_ice(:) ) 
            j = ppMicroZooPlankton(iiZ5,i)
            PELSURFACE(j,:) =  PELSURFACE(j,:) - flux_pel_ice(:)  
            ! map the flux into a 3D temporary array
            tmpflux(SRFindices) = PELSURFACE(j,:) / Depth(SRFindices)
            call flux_vector(iiPel, j, j, tmpflux(:) )
       end do
    end if



#ifdef DEBUG
   LEVEL2 'F3 flux ',flux_pel_ice_O3
#endif
 
!EOC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Other Seaice diagnostics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

end subroutine SeaicetoPelCoup
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
