#include "DEBUG.h"
#include "cppdefs.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Seaicealgae
!
! DESCRIPTION
!   Parameter values for the sea ice algae group
!
! !INTERFACE
  module mem_SeaicetoPel
!
! !USES:

  use global_mem
  use mem, only: iiSeaiceAlgae,iiSeaiceDetritus,iiSeaiceBacteria,iiSeaiceZoo
  use mem, only: iiP1,iiP2,iiZ5,iiB1,iiR1,iiR6
!
!
! !AUTHORS
!   Letizia Tedesco and Marcello Vichi
!
! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is not public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Sea ice - pelagic coupling PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! NAME         [UNIT]/KIND            DESCRIPTION
  !        :     --------- Physical parameters -----------------
  !  p_epsIce    [m-1]          Attenuation coefficient for pure ice
  !  CalcSeaiceAtmFlux logical  Turn on atmospheric nutrient flux
  !  p_atmN1     [?]            Atmospheric phosphate flux
  !  p_atmN3     [?]            Atmospheric nitrate flux
  !  p_atmN4     [?]            Atmospheric ammonium flux
  !  CalcSeaiceFloodFlux logical  Turn on nutrient flux due to flooding
  !  p_floN1     [?]            Flooding phosphate flux
  !  p_floN3     [?]            Flooding nitrate flux
  !  p_floN4     [?]            Flooding ammonium flux
  !  CalcSeaiceAddPelFlux logical Turn on additional nutrient flux from 
  !                               pelagic (note: to be used when there is no 
  !                               benthic model)
  !  p_pelN1     [mmol P/m3]      Reference PO4 concentration for restoration
  !  p_relaxN1   [d-1]            Relaxation time scale for restoration
  real(RLEN)  :: p_epsIce
  logical     :: CalcSeaiceAtmFlux, CalcSeaiceFloodFlux, CalcSeaiceAddPelFlux
  real(RLEN)  :: p_atmN1, p_atmN3, p_atmN4
  real(RLEN)  :: p_floN1, p_floN3, p_floN4, p_floN5
  real(RLEN)  :: p_pelN1, p_relaxN1

  integer,dimension(iiSeaiceAlgae)   :: PPY
  integer,dimension(iiSeaiceDetritus):: DET
  integer,dimension(iiSeaiceZoo)     :: ZOO
  integer,dimension(iiSeaiceBacteria):: BAC
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitSeaicetoPel, SeaicetoPelCoup
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitSeaicetoPel()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Seaicecoup_parameters/p_epsIce, &
    CalcSeaiceAtmFlux, p_atmN1, p_atmN3, p_atmN4, &
    CalcSeaiceFloodFlux, p_floN1, p_floN3, p_floN4, p_floN5, &
    CalcSeaiceAddPelFlux, p_pelN1, p_relaxN1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! set default values
    p_epsIce = 1.5_RLEN
    CalcSeaiceAtmFlux = .FALSE.
    p_atmN1 = ZERO
    p_atmN3 = ZERO
    p_atmN4 = ZERO
    CalcSeaiceFloodFlux = .FALSE.
    p_floN1 = ZERO
    p_floN3 = ZERO
    p_floN4 = ZERO
    p_floN5 = ZERO
    CalcSeaiceAddPelFlux = .FALSE.
    ! set the arrays of correspondance between sea ice and pelagic components
    !       iiS1  iiS2
    PPY = (/iiP1, iiP2/)
    !       iiU1 iiU6
    DET = (/iiR1,iiR6/)
    !       iiX1
    ZOO = (/iiZ5/)
    !       iiT1
    BAC = (/iiB1/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading Sea ice-pelagic parameters.."
    open(NMLUNIT,file='Seaice_Coupling.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=Seaicecoup_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Seaicecoup_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"SeaicetopelCoup.F90","Seaice_Coupling.nml")
101 call error_msg_prn(NML_READ,"SeaicetopelCoup.F90","Seaicecoup_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitSeaicetoPel

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
   use mem_PAR,    only: p_PAR, p_epsR6
   use mem_Param,  only: CalcSeaiceAlgae,CalcSeaiceZoo,CalcSeaiceBacteria 
   use mem_Param,  only: p_small
   use constants,  only: E2W, SEC_PER_DAY
   use mem
   use api_bfm,    only: SRFindices, BOTindices 
   use mem_SeaiceAlgae, only: p_epsSAL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   IMPLICIT NONE
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_N1,flux_pel_ice_N3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_N4,flux_pel_ice_N5
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice_O2,flux_pel_ice_O3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_pel_ice, flux_atm_N1, flux_atm_N3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_lat_O3, flux_atm_N4
   real(RLEN), dimension(NO_BOXES_XY) :: flux_flood_N1,flux_flood_N3
   real(RLEN), dimension(NO_BOXES_XY) :: flux_flood_N4,flux_flood_N5
   real(RLEN), dimension(NO_BOXES_XY) :: flux_bott_pel_N1
   real(RLEN), dimension(NO_BOXES_XY) :: I1p_tilde,I3n_tilde,I4n_tilde,I5s_tilde
   real(RLEN), dimension(NO_BOXES_XY) :: F2o_tilde,F3c_tilde
   real(RLEN), dimension(NO_BOXES_XY) :: I1p_star,I3n_star,I4n_star,I5s_star
   real(RLEN), dimension(NO_BOXES_XY) :: F2o_star,F3c_star
   real(RLEN), dimension(NO_BOXES_XY) :: EICE1D,epsIce,EIRsrf

   integer                            :: i,j,p
   real(RLEN), dimension(:), pointer  :: lcl_PelagicVar,lcl_SeaiceVar
   real(RLEN)                         :: tmpflux(NO_BOXES)
   real(RLEN)                         :: delta
   real(RLEN), external               :: GetDelta
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!

    delta = GetDelta()

    ! set to zero all the fluxes in case there is no ice
    flux_pel_ice_N1(:)= ZERO
    flux_pel_ice_N3(:)= ZERO
    flux_pel_ice_N4(:)= ZERO
    flux_pel_ice_N5(:)= ZERO
    flux_pel_ice_O2(:)= ZERO
    flux_pel_ice_O3(:)= ZERO

    flux_bott_pel_N1(:)=ZERO

    tmpflux(:)= ZERO

    ! prescription of atmospheric nutrient fluxes
    flux_atm_N1(:)=ZERO
    flux_atm_N3(:)=ZERO
    flux_atm_N4(:)=ZERO

    ! prescription of flooding fluxes
    flux_flood_N1(:)=ZERO
    flux_flood_N3(:)=ZERO
    flux_flood_N4(:)=ZERO
    flux_flood_N5(:)=ZERO
    flux_lat_O3(:)=ZERO

       if (CalcSeaiceAtmFlux) then
          ! Enhancement factor for nutrient flux from atmosphere (to be defined by the user, below an example)
          flux_atm_N1(:)= EDS(:)*p_atmN1*EHB(:)
          flux_atm_N3(:)= EDS(:)*p_atmN3*EHB(:)
          flux_atm_N4(:)= EDS(:)*p_atmN4*EHB(:)
       end if

       if (CalcSeaiceFloodFlux) then
          ! Enhancement factor for nutrient flux from flooding events (to be defined by the user, below an example)
          flux_flood_N1(:)=EDS(:)*p_floN1*EHB(:)                          
          flux_flood_N3(:)=EDS(:)*p_floN3*EHB(:)
          flux_flood_N4(:)=EDS(:)*p_floN4*EHB(:)                         
          flux_flood_N5(:)=EDS(:)*p_floN5*EHB(:)
          !Additional lateral flux of DIC (to be defined by the user, below an example)
          flux_lat_O3(:)=-0.5_RLEN*min(ZERO, (F3c(:)-(9.0_RLEN*EHB(:))))
       end if

       if (CalcSeaiceAddPelFlux) then
          ! additional nutrient flux to simulate bottom inputs when no benthic model is available
          flux_bott_pel_N1(:) = p_relaxN1*max(ZERO,p_pelN1 - N1p(BOTindices))*Depth(BOTindices)
       end if

    ! Assign surface EIR to a local variables
    EIRsrf(:) = EIR(SRFindices)
    ! compute attenuation in the BAL. Note that this is non dimensional
    ! as the sea ice components are integrated over depth
    epsIce(:) = p_epsIce*EHB(:) + p_epsR6*U6c(:)
    do i = 1 , ( iiSeaiceAlgae)
      lcl_SeaiceVar => SeaiceAlgae(i,iiL)
      epsIce(:) = epsIce(:) + p_epsSAL(i) * lcl_SeaiceVar
    end do

    ! fluxes to the BAL from atmosphere, flooding and pelagic
    where (EHB(:) > ZERO)
       ! volume concentration
       I1p_tilde(:) = I1p(:)/EHB(:)
       I3n_tilde(:) = I3n(:)/EHB(:)
       I4n_tilde(:) = I4n(:)/EHB(:)
       I5s_tilde(:) = I5s(:)/EHB(:)
       F2o_tilde(:) = F2o(:)/EHB(:)
       F3c_tilde(:) = F3c(:)/EHB(:)
   
       I1p_star(:) = max(ZERO, N1p(SRFindices) - I1p_tilde(:))
       I3n_star(:) = max(ZERO, N3n(SRFindices) - I3n_tilde(:))
       I4n_star(:) = max(ZERO, N4n(SRFindices) - I4n_tilde(:))
       I5s_star(:) = max(ZERO, N5s(SRFindices) - I5s_tilde(:))
       F2o_star(:) = max(ZERO, O2o(SRFindices) - F2o_tilde(:))
       F3c_star(:) = max(ZERO, O3c(SRFindices) - F3c_tilde(:))


       ! flux is divided in positive (to-ice) and negative (to-water)
       flux_pel_ice_N1(:)= max(ZERO,EDH(:)*I1p_star(:)) &
                           + min(ZERO, EDH(:)*I1p_tilde(:))
       flux_pel_ice_N3(:)= max(ZERO,EDH(:)*I3n_star(:)) &
                           + min(ZERO, EDH(:)*I3n_tilde(:))
        flux_pel_ice_N4(:)= max(ZERO,EDH(:)*I4n_star(:)) &
                           + min(ZERO, EDH(:)*I4n_tilde(:))
       flux_pel_ice_N5(:)= max(ZERO,EDH(:)*I5s_star(:)) &
                           + min(ZERO, EDH(:)*I5s_tilde(:))
       flux_pel_ice_O2(:)= max(ZERO,EDH(:)*F2o_star(:)) &
                           + min(ZERO, EDH(:)*F2o_tilde(:))
       flux_pel_ice_O3(:)= max(ZERO,EDH(:)*F3c_star(:)) &
                           + min(ZERO, EDH(:)*F3c_tilde(:))

       ! compute irradiance at the bottom of sea ice (top of pelagic)
       ! accounting for seaice algae
       ! it is assumed that irradiance in the BAL is located at the middle
       EICE1D(:) = ONE
       EIRsrf(:) = EIRsrf(:)*(ONE-EICE1D(:)) + (EIB(:)*EICE1D(:))*exp(-0.5_RLEN*epsIce(:))

    elsewhere

       ! additional flux to the water when there is no ice to ensure complete emptyness
       flux_pel_ice_N1(:) = min(ZERO,p_small-I1p(:))/delta
       flux_pel_ice_N3(:) = min(ZERO,p_small-I3n(:))/delta
       flux_pel_ice_N4(:) = min(ZERO,p_small-I4n(:))/delta
       flux_pel_ice_N5(:) = min(ZERO,p_small-I5s(:))/delta
       flux_pel_ice_O2(:) = min(ZERO,p_small-F2o(:))/delta
       flux_pel_ice_O3(:) = min(ZERO,p_small-F3c(:))/delta
       ! first release everything then force the value to be p_small
       I1p(:)=max(p_small,I1p)
       I3n(:)=max(p_small,I3n)
       I4n(:)=max(p_small,I4n)
       I5s(:)=max(p_small,I5s)
       F2o(:)=max(p_small,F2o)
       F3c(:)=max(p_small,F3c)

       ! compute irradiance at the bottom of sea ice (top of pelagic)
       ! accounting for seaice algae
       ! it is assumed that irradiance in the BAL is located at the middle
       EICE1D(:) = ZERO
       EIRsrf(:) = EIRsrf(:)*(ONE-EICE1D(:)) + (EIB(:)*EICE1D(:))*exp(-0.5_RLEN*epsIce(:))
    end where

    ! assign back the modified surface irradiance at the sea ice-pelagic interface
    !MAV: to check if this surface value is then propagated at depth 
    EIR(SRFindices) = EIRsrf

    ! add fluxes to D2SOURCE_ICE (on the diagonal)
    call flux_vector( iiIce, ppI1p,ppI1p, flux_pel_ice_N1 + flux_atm_N1 + flux_flood_N1) 
    call flux_vector( iiIce, ppI3n,ppI3n, flux_pel_ice_N3 + flux_atm_N3 + flux_flood_N3) 
    call flux_vector( iiIce, ppI4n,ppI4n, flux_pel_ice_N4 + flux_flood_N4) 
    call flux_vector( iiIce, ppI5s,ppI5s, flux_pel_ice_N5 + flux_flood_N5) 
    call flux_vector( iiIce, ppF2o,ppF2o, flux_pel_ice_O2 ) 
    call flux_vector( iiIce, ppF3c,ppF3c, flux_pel_ice_O3 + flux_lat_O3) 

    ! assign fluxes to boundary variables 
    ! Compute rates for the pelagic variables and assign to D3SOURCE
    jsurN1p(:)=jsurN1p(:) - flux_pel_ice_N1(:)
    tmpflux(SRFindices) = jsurN1p(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN1p,ppN1p, tmpflux(:) )

    jbotN1p(:)=jbotN1p(:) + flux_bott_pel_n1(:)
    tmpflux(BOTindices) = jbotN1p(:)/Depth(BOTindices)
    call flux_vector( iiPel, ppN1p, ppN1p, tmpflux )

    jsurN3n(:)=jsurN3n(:) - flux_pel_ice_N3(:)
    tmpflux(SRFindices) = jsurN3n(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN3n,ppN3n, tmpflux(:) )

    jsurN4n(:)=jsurN4n(:) - flux_pel_ice_N4(:)
    tmpflux(SRFindices) = jsurN4n(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN4n,ppN4n, tmpflux(:) )

    jsurN5s(:)=jsurN5s(:) - flux_pel_ice_N5(:)
    tmpflux(SRFindices) = jsurN5s(:) / Depth(SRFindices)
    call flux_vector(iiPel,ppN5s,ppN5s, tmpflux(:) )

    ! oxygen and CO2 fluxes are assigned in the pelagic routines
    jsurO2o(:)=jsurO2o(:) - flux_pel_ice_O2(:)
    jsurO3c(:)=jsurO3c(:) - flux_pel_ice_O3(:)

   do j = 1, iiSeaiceAlgae
     if (CalcSeaiceAlgae(j)) then
       do i = 1,iiLastElement
          lcl_SeaiceVar => SeaiceAlgae(j,i)
          lcl_PelagicVar => PhytoPlankton(PPY(j),i)
          if( associated(lcl_SeaiceVar) .AND. associated(lcl_PelagicVar) ) then
             where (EHB(:)>ZERO)
                flux_pel_ice(:) = max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices)*EVB(:))  &
                     + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
             elsewhere
                flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/delta
             end where
          end if
           ! add the flux and assign it to the boundary variable of the pelagic system
            p = ppSeaiceAlgae(j,i)
            if (p>0) call flux_vector( iiIce, p, p, flux_pel_ice(:) ) 
            p = ppPhytoPlankton(PPY(j),i)
            if (p>0) then
               PELSURFACE(p,:) =  PELSURFACE(p,:) - flux_pel_ice(:)  
               ! map the flux into a 3D temporary array and send to D3SOURCE
               tmpflux(SRFindices) = PELSURFACE(p,:) / Depth(SRFindices)
               call flux_vector(iiPel, p, p, tmpflux(:) )
            end if
       end do
     end if
   end do

   ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! Detritus Fluxes to Pelagic from Seaice
   ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   do j = 1, iiSeaiceDetritus
      do i = 1,iiLastElement
         lcl_SeaiceVar => SeaiceDetritus(j,i)
         lcl_PelagicVar => PelDetritus(DET(j),i)
         if( associated(lcl_SeaiceVar) .AND. associated(lcl_PelagicVar) ) then
            where (EHB(:)>ZERO)
               flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices)*EVB(:))  &
                    + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
            elsewhere
               flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/delta
            end where
         end if
         ! add the flux and assign it to the boundary variable of the pelagic system
         p = ppSeaiceDetritus(j,i)
         if (p>0) call flux_vector( iiIce, p, p, flux_pel_ice(:) ) 
         p = ppPelDetritus(DET(j),i)
         if (p>0) then
            PELSURFACE(p,:) =  PELSURFACE(p,:) - flux_pel_ice(:)  
            ! map the flux into a 3D temporary array
            tmpflux(SRFindices) = PELSURFACE(p,:) / Depth(SRFindices)
            call flux_vector(iiPel, p, p, tmpflux(:) )
         end if
      end do
   end do

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   !Bacteria Flux to Pelagic from Seaice
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   do j = 1, iiSeaiceBacteria
      if (CalcSeaiceBacteria(j)) then
        do i = 1,iiLastElement
           lcl_SeaiceVar => SeaiceBacteria(j,i)
           lcl_PelagicVar => PelBacteria(BAC(j),i)
           if( associated(lcl_SeaiceVar) .AND. associated(lcl_PelagicVar) ) then
              where (EHB(:)>ZERO)
                 flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices)*EVB(:))  &
                      + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
              elsewhere
                 flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/delta
              end where
           end if
           ! add the flux and assign it to the boundary variable of the pelagic system
           p = ppSeaiceBacteria(j,i)
           if (p>0) call flux_vector( iiIce, p,p, flux_pel_ice(:) ) 
           p = ppPelBacteria(BAC(j),i)
           if (p>0) then
              PELSURFACE(p,:) =  PELSURFACE(p,:) - flux_pel_ice(:)  
              ! map the flux into a 3D temporary array
              tmpflux(SRFindices) = PELSURFACE(p,:) / Depth(SRFindices)
              call flux_vector(iiPel, p, p, tmpflux(:) )
           end if
        end do
      end if
   end do

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !Microzooplankton Flux to Pelagic from Seaice
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   do j = 1, iiSeaiceZoo
      if (CalcSeaiceZoo(j)) then
        do i = 1,iiLastElement
           lcl_SeaiceVar => SeaiceZoo(j,i)
           lcl_PelagicVar => MicroZooPlankton(ZOO(j),i)
           if( associated(lcl_SeaiceVar) .AND. associated(lcl_PelagicVar) ) then
              where (EHB(:)>ZERO)
                 flux_pel_ice(:)= max(ZERO, EDH(:)*lcl_PelagicVar(SRFindices)*EVB(:))  &
                      + min(ZERO, EDH(:)*lcl_SeaiceVar(:)/EHB(:))
              elsewhere
                 flux_pel_ice(:) = min(ZERO,p_small-lcl_SeaiceVar(:))/delta
              end where
           end if
           ! add the flux and assign it to the boundary variable of the pelagic system
           p = ppSeaiceZoo(j,i)
           if (p>0) call flux_vector( iiIce, p,p, flux_pel_ice(:) ) 
           p = ppMicroZooPlankton(ZOO(j),i)
           if (p>0) then
              PELSURFACE(p,:) =  PELSURFACE(p,:) - flux_pel_ice(:)  
              ! map the flux into a 3D temporary array
              tmpflux(SRFindices) = PELSURFACE(p,:) / Depth(SRFindices)
              call flux_vector(iiPel, p, p, tmpflux(:) )
           end if
        end do
      end if
   end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Other Seaice diagnostics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#ifdef DEBUG
   LEVEL2 'F3 flux ',flux_pel_ice_O3
   LEVEL2 'flux_pel_ice_N1',flux_pel_ice_N1
   LEVEL2 'F2o',F2o
   LEVEL2 'flux_pel_ice_O2',flux_pel_ice_O2
   LEVEL2 'EICE=',EICE
   LEVEL2 'EHB=',EHB
#endif
 
!EOC

end subroutine SeaicetoPelCoup
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

end module mem_SeaicetoPel
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
