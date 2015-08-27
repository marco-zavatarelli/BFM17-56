!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Initialize
!
! DESCRIPTION
!   Initialization of model
!   Allocation of memory for variables, reading of data files 
!
! !INTERFACE
  SUBROUTINE Initialize
!
! USES:
  use mem, only: InitializeModel,ppMicroZooplankton,ppMesoZooPlankton, &
                 iiMicroZooplankton,iiMesoZooPlankton,NO_BOXES,        &
                 iiN,iiP,qpcMEZ,qncMEZ,qpcMIZ,qncMIZ
  use mem_Param
  use mem_PelGlobal
  use mem_PelChem
  use mem_PelBac
  use mem_MesoZoo
  use mem_MicroZoo
  use mem_Phyto
  use mem_PAR
  use mem_Settling
#ifdef INCLUDE_BEN
  use mem_BenOrganism
  use mem_FilterFeeder
  use mem_BenBac
  use mem_Bioturbation
  use mem_BenthicReturn1
  use mem_BenthicReturn2
  use mem_BenthicNutrient3
  use mem_BenAmmonium
  use mem_BenNitrate
  use mem_BenOxygen
  use mem_BenAnoxic
  use mem_BenDenitriDepth
  use mem_BenPhosphate
  use mem_BenSilica
  use mem_BenQ1Transport
  use mem_ControlBennutBuffers
#endif
#ifdef INCLUDE_PELCO2
  use mem_CO2
#endif
#ifdef INCLUDE_BENCO2
  use mem_BenCO2Transport
  use mem_BenAlkalinity
#endif
#ifdef INCLUDE_SILT
  use mem_Silt
#endif
#ifdef INCLUDE_SEAICE
  use mem_SeaiceAlgae
  use mem_SeaiceBac
  use mem_SeaiceZoo
  use mem_SeaicetoPel
#endif

!  
!
! !AUTHORS
!   mfstep ERSEM team
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
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
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      InitializeModel=0

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Allocate Memory for All global variables
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call AllocateMem

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Read all data files:(namelist files)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call InitParam
      call InitPelGlobal
      call InitPelChem
      call InitPelBac
      call InitMesoZoo
      call InitMicroZoo
      call InitPhyto
      call InitPAR
      call InitSettling
#ifdef INCLUDE_BEN
      ! Benthic initialization is done only if there is an active model
      ! When INCLUDE_BEN is defined, 
      ! CalcBenthicFlag=0 is used to test the benthic memory only
      if ( CalcBenthicFlag > 0 ) then
         call InitBenOrganism
         call InitFilterFeeder
         call InitBenBac
         call InitBioturbation
         call InitBenthicReturn1
         call InitBenthicReturn2
         call InitBenthicNutrient3
         call InitBenAmmonium
         call InitBenNitrate
         call InitBenOxygen
         call InitBenAnoxic
         call InitBenDenitriDepth
         call InitBenPhosphate
         call InitBenSilica
         call InitBenQ1Transport
         call InitControlBennutBuffers
#ifdef INCLUDE_BENCO2
         call InitBenCO2Transport
         call InitBenAlkalinity
#endif
      end if
#endif
#ifdef INCLUDE_PELCO2
      call InitCO2
#endif
#ifdef INCLUDE_SILT
      call InitSilt
#endif
#ifdef INCLUDE_SEAICE
      call InitSeaiceAlgae
      call InitSeaiceBac
      call InitSeaiceZoo
      call InitSeaicetoPel
#endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Read all other Init* files
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call InitTransportStateTypes
      call InitBoxParams

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Initialize nutrient quota in Microzooplankton and Mesozooplankton
      ! with the parameter values in the namelists.
      ! These are constant values when running with fixed quota in zoo
      ! In case of variable quota these values are recomputed every time-step
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       do i = 1 , ( iiMicroZooPlankton)
         if ( ppMicroZooPlankton(i,iiP) == 0 ) qpcMIZ(i,:)  =  p_qpcMIZ(i) 
         if ( ppMicroZooPlankton(i,iiN) == 0 ) qncMIZ(i,:)  =  p_qncMIZ(i)
       end do
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Nutrient quota in omnivorous and herbivorous mesozooplankton
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       do i = 1 , ( iiMesoZooPlankton)
         if ( ppMesoZooPlankton(i,iiP) == 0 ) qpcMEZ(i,:)  =   p_qpcMEZ(i)
         if ( ppMesoZooPlankton(i,iiN) == 0 ) qncMEZ(i,:)  =   p_qncMEZ(i)
       end do


    END SUBROUTINE Initialize
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
