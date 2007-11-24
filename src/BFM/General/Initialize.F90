!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Initialize
!
! DESCRIPTION
!   Initialization of model
!   Allocation of memory for variables, reading of data files 

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE Initialize
!
! USES:
  use mem, only: InitializeModel,ppMicroZooplankton,ppMesoZooPlankton, &
                 MicroZooplankton,MesoZooPlankton,iiMicroZooplankton,  &
                 iiMesoZooPlankton,NO_BOXES,iiN,iiP,qpZc,qnZc,qp_mz,qn_mz
  use mem_Param
  use mem_WindOxReaeration_3
  use mem_PelGlobal
  use mem_PelChem
  use mem_CO2
  use mem_PelBac
  use mem_MesoZoo,p_qnMc=>p_qnc,p_qpMc=>p_qpc
  use mem_MicroZoo
  use mem_Phyto
  use mem_PhotoAvailableRadiation
  use mem_LightAdaptation
  use mem_Settling
#ifdef BFM_BENTHIC
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
#ifdef BFM_SI
  use mem_SeaiceAlgae
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
      call InitWindOxReaeration_3
      call InitPelGlobal
      call InitPelChem
      call InitCO2
      call InitPelBac
      call InitMesoZoo
      call InitMicroZoo
      call InitPhyto
      call InitPhotoAvailableRadiation
      call InitLightAdaptation
      call InitSettling
#ifdef BFM_BENTHIC
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
#endif
#ifdef BFM_SI
      call InitSeaiceAlgae
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
         if ( ppMicroZooPlankton(i,iiP) == 0 ) qp_mz(i,:)  =  p_qp_mz(i) 
         if ( ppMicroZooPlankton(i,iiN) == 0 ) qn_mz(i,:)  =  p_qn_mz(i)
       end do

       do i = 1 , ( iiMesoZooPlankton)
         if ( ppMesoZooPlankton(i,iiP) == 0 ) qpZc(i,:)  =   p_qpMc(i)
         if ( ppMesoZooPlankton(i,iiN) == 0 ) qnZc(i,:)  =   p_qnMc(i)
       end do


    END SUBROUTINE Initialize
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
