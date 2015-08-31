!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PAR
!
! DESCRIPTION
!
! PROCESS VALUES
!	PAR -- Irradiance in Water
!
! FILE
!
!
! !INTERFACE
  module mem_PAR
!
! !USES:

  use global_mem
  use mem,  ONLY: iiPhytoPlankton, NO_BOXES

!  
!
! !AUTHORS
!
!
! !REVISION_HISTORY
!   !
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
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! PAR PARAMETERS (read from nml)
  ! LightPeriodFlag               numeric Choose the light averaging period
  !                                       1 = Instantanous irradiance
  !                                       2 = Daily average
  !                                       3 = Daylight average with explicit
  !                                           photoperiod
  ! LightLocationFlag             numeric Choose the parameterization of light
  !                                       location in the discrete grid
  !                                       1 = Light at the top of the cell
  !                                       2 = Light in the middle of the cell
  !                                       3 = Average Light in the cell
  ! Light attenuation parameters
  !  ChlAttenFlag                 numeric Choose the type of Chl attenuation
  !                                       1 = Broadband linear 
  !                                       2 = 3-band tabulated 
  !  p_PAR        [-]           Fraction of Photosynthetically Available Radiation
  !  p_eps0       [1/m]         Background extinction coefficient
  !  p_epsESS     [m2/g]        Specific attenuation coefficient of
  !                             suspended sediments
  !  p_epsR6      [m2/mgC]      Specific attenuation coefficient of particulate
  !                             detritus
  !                       
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: LightPeriodFlag=1
  integer  :: LightLocationFlag=3
  real(RLEN) :: &
      p_PAR=0.40_RLEN,       &
      p_eps0=0.0435_RLEN,    &
      p_epsIR=0.35_RLEN,     &
      p_epsESS=0.04e-3_RLEN, &
      p_epsR6=0.1e-3_RLEN
  integer  :: ChlAttenFlag=1      
  ! shared variables
  real(RLEN), dimension(3,61) :: xepsRGB  ! Tabulated attenuation coefficients 
                                          ! for RGB parameterization(Lengaigne)
  real(RLEN)                  :: p_PARRGB ! portion of PAR for RGB
  ! arrays for the 3-band parameterization
  real(RLEN),allocatable,dimension(:) :: B_eps, G_eps, R_eps
  real(RLEN),allocatable,dimension(:) :: EIRB, EIRG, EIRR

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")
  public InitPAR

  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPAR()
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer :: AllocStatus

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /PAR_parameters/ LightPeriodFlag, LightLocationFlag,         &
            p_PAR, p_eps0, p_epsIR, p_epsESS, p_epsR6, ChlAttenFlag
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading PAR parameters.."
    open(NMLUNIT,file='Pelagic_Environment.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=PAR_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=PAR_parameters)
    write(LOGUNIT,*) "#  PAR specifications:"
    write(LOGUNIT,*) "#  Chl Attenuation Flag p_ChlAttenFlag =",ChlAttenFlag
    select case (ChlAttenFlag) 
      case (2) 
         write(LOGUNIT,*) "#   Use 3 bands tabulated RGB (Lengaigne et al, 2007)"
         ! Initialize the tabulated values and 3-band arrays
         call ChlAttenuation(xepsRGB)
         allocate(B_eps(NO_BOXES),stat=AllocStatus)
         if (AllocStatus  /= 0) stop "error allocating B_eps"
         allocate(EIRB(NO_BOXES),stat=AllocStatus)
         if (AllocStatus  /= 0) stop "error allocating EIRB"
         allocate(G_eps(NO_BOXES),stat=AllocStatus)
         if (AllocStatus  /= 0) stop "error allocating G_eps"
         allocate(EIRG(NO_BOXES),stat=AllocStatus)
         if (AllocStatus  /= 0) stop "error allocating EIGB"
         allocate(R_eps(NO_BOXES),stat=AllocStatus)
         if (AllocStatus  /= 0) stop "error allocating R_eps"
         allocate(EIRR(NO_BOXES),stat=AllocStatus)
         if (AllocStatus  /= 0) stop "error allocating EIRR"
         p_PARRGB = p_PAR/3._RLEN ! the visible part is divided equally in 3 bands
      case default 
         write(LOGUNIT,*) "#   Linear Chl-Specific attenuation coefficient"
    end select
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPAR","Pelagic_Ecology.nml")
101 call error_msg_prn(NML_READ,"InitPAR","PAR_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPAR
  end module mem_PAR
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
