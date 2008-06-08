!$Id: $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_bfm --- BFM bio model \label{sec:bio_bfm}
!
! !INTERFACE:
   module bio_bfm
!
! !DESCRIPTION:
!
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_bfm, pointers_gotm_bfm,            &
          var_info_bfm, envforcing_bfm, do_bio_bfm,   &
          allocate_memory_bfm,reset_diagonal,         &
          test_on_negative_states, end_bio_bfm,       &
          do_bfm_river_loads, settling_vel_bfm,       &
          CalcVertFluxAtLev
!
!
! !PRIVATE DATA MEMBERS:
   REALTYPE,public,dimension(:),allocatable :: cdepth,wx
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
!  $Log: $
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the template bio module
!
! !INTERFACE:
   subroutine init_bio_bfm(nlev,out_unit)
!
! !DESCRIPTION:
!  Here, the main communication of array dimensions between GOTM
!  and BFM is done.
!
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_D2_BOX_STATES, NO_BOXES_XY,       &
                  NO_D2_BOX_DIAGNOSS, NO_D3_BOX_DIAGNOSS,&
                  NO_D2_BOX_FLUX, NO_D3_BOX_FLUX,&
                  NO_STATES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: nlev
   integer,          intent(in)   :: out_unit
!
   integer :: i,rc,n
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_bfm'


   ! BFM  --> GOTM
   numc  = NO_D3_BOX_STATES
   numbc = NO_D2_BOX_STATES
   numc_diag  = NO_D3_BOX_DIAGNOSS
   numbc_diag = NO_D2_BOX_DIAGNOSS
   numc_flux  = NO_D3_BOX_FLUX
   numbc_flux = NO_D2_BOX_FLUX
   ! numcc is the number of transported variables
   numcc = numc

   ! GOTM --> BFM
   NO_BOXES_X  = 1
   NO_BOXES_Y  = 1
   NO_BOXES_Z  = nlev
   NO_BOXES    = NO_BOXES_X * NO_BOXES_Y * NO_BOXES_Z
   NO_BOXES_XY = NO_BOXES_X * NO_BOXES_Y
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES +   &
                 NO_D2_BOX_STATES * NO_BOXES_XY
   !LOGUNIT = out_unit

   ! assign the indices of the surface and bottom grid-points
   ! GOTM vertical index is from the bottom to the surface
   allocate(SRFindices(NO_BOXES_XY))
   allocate(BOTindices(NO_BOXES_XY))
   SRFindices(:) = nlev
   BOTindices(:) = 1
   
   LEVEL3 'pelagic variables =',numc
   LEVEL3 'pelagic transported variables =',numcc
   LEVEL3 'benthic variables =',numbc
   LEVEL3 'pelagic variables prepared for output',numc_diag
   LEVEL3 'benthic variables prepared for output',numbc_diag
   LEVEL3 'NO_BOXES_X=',NO_BOXES_X
   LEVEL3 'NO_BOXES_Y=',NO_BOXES_Y
   LEVEL3 'NO_BOXES_Z=',NO_BOXES_Z
   LEVEL3 'NO_BOXES=',NO_BOXES
   LEVEL3 'NO_BOXES_XY=',NO_BOXES_XY
   LEVEL3 'NO_STATES=',NO_STATES
   LEVEL3 'Step 1 of GOTM <-> BFM initialisation done ...'

!  sfl=_ZERO_
!  sfl_read=_ZERO_
   allocate(wx(1:NO_BOXES),stat=rc)
   allocate(cdepth(1:NO_BOXES),stat=rc)
   return

   end subroutine init_bio_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialize BFM and GETM shared memory
!
! !INTERFACE:
   subroutine pointers_gotm_bfm()
!
! !DESCRIPTION:
! Allocate pointers to GOTM memory
!
! !USES:
   use mem, only: D3STATE,D3SOURCE,D3SINK,D3STATETYPE, &
                  D3DIAGNOS,D2STATE,D2SOURCE,D2SINK,   &
                  D2STATETYPE,NO_BOXES,NO_BOXES_XY,    &
                  D2DIAGNOS,NO_D2_BOX_STATES,          &
                  NO_D2_BOX_DIAGNOSS

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:

   !---------------------------------------------
   ! Pelagic pointers
   !---------------------------------------------
   D3STATE  => cc(:,1:NO_BOXES)
   D3SOURCE => pp(:,:,1:NO_BOXES)
   D3SINK   => dd(:,:,1:NO_BOXES)
   D3STATETYPE => pelvar_type
   !---------------------------------------------
   ! 3D diagnostics
   !---------------------------------------------
   if (numc_diag > 0) D3DIAGNOS => diag(:,1:NO_BOXES)

   !---------------------------------------------
   ! Benthic pointers
   !---------------------------------------------
   if (bio_setup >=2 ) then
      D2STATE  => ccb(:,1:NO_BOXES_XY)
      D2SOURCE => ppb(:,:,1:NO_BOXES_XY)
      D2SINK   => ddb(:,:,1:NO_BOXES_XY)
      D2STATETYPE => benvar_type
   else
      ! allocate memory anyhow to avoid problems with BFM allocation
      allocate(D2STATE(1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
      allocate(D2SOURCE(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
      allocate(D2SINK(1:NO_D2_BOX_STATES,1:NO_D2_BOX_STATES,1:NO_BOXES_XY))
      allocate(D2STATETYPE(1:NO_D2_BOX_STATES ))
   end if

   !---------------------------------------------
   ! allocate 2D diagnostics (not only benthic)
   !---------------------------------------------
   if (numbc_diag>0) D2DIAGNOS => diagb(:,1:NO_BOXES_XY)

   end subroutine pointers_gotm_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Providing info on variables
!
! !INTERFACE:
   subroutine var_info_bfm()
!
! !DESCRIPTION:
!  This subroutine provides information on the variables. To be used
!  when storing data in NetCDF files.
!
! !USES:
   use mem
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call set_var_info_bfm
   return
   end subroutine var_info_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcing used in the BFM
!
! !INTERFACE
   subroutine envforcing_bfm(nlev, bioshade_feedback, h)
!
! !DESCRIPTION
!
! !USES
! BFM modules
use constants, ONLY: E2W
use mem_Param, ONLY: p_eps0, p_epsESS, p_PAR,p_small
use mem,       ONLY: NO_BOXES, R6c, PhytoPlankton, xEPS, ESS, ERHO, &
                     iiPhytoPlankton, iiL, Chla, ETW, ESW, &
                     Depth, EIR, ABIO_eps, EWIND, ETAUB
use mem,       ONLY: Volume, Area, Area2D
#ifdef INCLUDE_SILT
use mem,       ONLY: R9x
#endif
use mem_Param,  ONLY: p_eps0, p_epsESS,p_poro
use global_interface,   ONLY: eTq
use bio_var,    ONLY: wind_gotm => wind, u_taub
use global_mem, only: ONE
! GOTM module for light extinction
use observations, only: A,g1,g2

IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,intent(in)                   :: nlev
   logical, intent(in)                  :: bioshade_feedback
   REALTYPE,intent(in)                  :: h(0:nlev)
!
! !OUTPUT PARAMETERS:

! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from a template by Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer             :: i,n
   REALTYPE            :: psilt
!EOP
!-----------------------------------------------------------------------
!BOC

#ifdef DEBUG
   LEVEL2 'calculating environmental forcings for the BFM'
#endif
   !---------------------------------------------
   ! Update the depths and volume of layers
   !---------------------------------------------
   ! Assign a virtual value for the grid-point area
   ! In the future, the actual value from GETM might be used
   ! (put here in case of adaptive grids)
   Area(:)   = ONE
   Area2D(:) = ONE
   Depth(:) = h(1:nlev)
   Volume(:) = Depth(:)*Area(:)
   ! cdepth is cumulative depth
   cdepth(NO_BOXES) = Depth(NO_BOXES)
    do n=NO_BOXES-1,1,-1
       cdepth(n)=depth(n)+ cdepth(n+1)
    enddo
   !---------------------------------------------
   ! Assign physical variables to the BFM
   !---------------------------------------------
   EWIND = wind_gotm
   ETAUB = u_taub
   ETW(:) = t(1:nlev)
   ESW(:) = s(1:nlev)
   ERHO(:) = rho(1:nlev)
   psilt=(p_poro(1) - 0.38662 )/ 0.00415
#ifdef INCLUDE_SILT
   ESS(:) = R9x(:)
#else
   ESS(:) = ZERO
#endif
            
   !---------------------------------------------
   ! Compute biological extinction coefficient
   ! The same abiotic extinction coefficients of 
   ! GOTM are used for the visible part
   !---------------------------------------------
   p_PAR = (_ONE_ - A)
   p_eps0 = _ONE_/g2
!MAV: this cannot be generic for all GOTM applications
! If we want to use this kind of parameterization, it has
! to be the same also in GOTM
!   p_eps0=1.17692307692-0.0307692307692*ESW(1)

   ! This part assumes that ISM extinction
   ! is an external forcing
   if (abioshade_(nlev) /= _ZERO_ ) then
     ABIO_eps(:) = abioshade_(1:nlev)
     p_eps0 = _ZERO_
   end if

   call  CalcVerticalExtinction( )

   !---------------------------------------------
   ! Note that irradiance in the BFM is in
   ! uE/m2/s and is defined at the top of each
   ! layer (the derivation of the average
   ! EIR for production is done in the
   ! Phytoplankton routines)
   !---------------------------------------------
   EIR(nlev) = max(p_small,p_PAR*I_0/E2W)
   do i=nlev,2,-1
     EIR(i-1) = EIR(i)*exp(-xEPS(i)*Depth(i))
   end do

   !---------------------------------------------
   ! bioshade is instead derived in the
   ! middle of the layer and it's non-dimensional
   ! (p_eps0 must be removed again because GOTM
   ! already considers it)
   !---------------------------------------------
   if (bioshade_feedback) &
     bioshade_(1:nlev) =  EIR(:)*exp(-(xEPS(:)-p_eps0)* &
                          Depth(:)*0.5)/ EIR(nlev)

#ifdef DEBUG
   LEVEL3 'ETW',ETW(nlev)
   LEVEL3 'ESW',ESW(nlev)
   LEVEL3 'EIR',EIR(nlev)
   LEVEL3 'EWIND',EWIND
   LEVEL3 'ETAUB',ETAUB
#endif
   end subroutine envforcing_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of the BFM model
!
! !INTERFACE
   subroutine do_bio_bfm(first,nlev)
!
! !DESCRIPTION
!  This subroutine is a wrapper for the computing core of the BFM
!  All the input arguments are dummies, for compatibility with gotm.
!
! !USES
   use mem_param, only: AssignAirPelFluxesInBFMFlag,AssignPelBenFluxesInBFMFlag
   use mem, only: sediPI, sediR6, sediR2,iiC,iiN,iiP,iiS,iiL, &
                  ppR2c, ppR6c, ppR6n, ppR6p, ppR6s, NO_BOXES_Z,   &
                  ppR1c, ppR1n, ppR1p,   &
                  ppO2o,ppN1p,ppN3n,ppN4n,ppN5s,ppN6r,  &
                  NO_D3_BOX_STATES, Depth,              &
                  ppPhytoPlankton,iiPhytoPlankton, &
                  PELBOTTOM, PELSURFACE, D3STATE, &
                  jK3G4n,jK13K3n
   use mem, only: N1p, N3n, N4n, N5s
   use mem_PelGlobal, only: p_rR6m
   use constants,  only: RLEN,SEC_PER_DAY
   use gotm_error_msg, only:gotm_error
   use util,          only  : Dirichlet, Neumann

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)         :: first
   integer, intent(in)         :: nlev
!
! !OUTPUT PARAMETERS:

!See "USE association" above
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!  from template by Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   logical                     :: ll_larger
   integer                     :: n,k,i,j,l
   REALTYPE                    :: nudg_vel ! nudging velocity
   REALTYPE                    :: RelaxTime = 30.0 ! days
   REALTYPE                    :: topm3psec,corr,Nloss
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Reset source term arrays 
   !---------------------------------------------
   call ResetFluxes

   !---------------------------------------------
   ! Compute BFM terms
   !---------------------------------------------
#ifdef INCLUDE_SILT
   call SiltDynamics
#endif
   call EcologyDynamics

!RUA specific to North Sea
   Nloss=jK3G4n(1)-jK13K3n(1);
   !---------------------------------------------
   ! Surface fluxes (mmol/m2/day)
   !---------------------------------------------
   nudg_vel=Depth(NO_BOXES_Z)/RelaxTime
   if ( .NOT. AssignAirPelFluxesInBFMFlag ) then
     select case (surface_flux_method)
        case (-1)! absolutely nothing
        case (1) ! constant
           PELSURFACE(ppN1p,1) = 0.01
           PELSURFACE(ppN3n,1) = 0.01
           PELSURFACE(ppN4n,1) = 0.01
           PELSURFACE(ppN5s,1) = 0.01
!RUA
!           sfl(ppN3n) =   0.12  *topm3psec
!           sfl(ppN4n) =   0.09  *topm3psec
!           sfl(ppN1p) =   0.0  !0.0
        case (0) ! nudging
           PELSURFACE(ppN1p,1) = max(_ZERO_,1.0-N1p(NO_BOXES_Z))*nudg_vel
           PELSURFACE(ppN3n,1) = max(_ZERO_,7.0-N3n(NO_BOXES_Z))*nudg_vel
           PELSURFACE(ppN4n,1) = max(_ZERO_,4.0-N4n(NO_BOXES_Z))*nudg_vel
           PELSURFACE(ppN5s,1) = max(_ZERO_,2.0-N5s(NO_BOXES_Z))*nudg_vel
        case (2) ! from file via sfl_read
           ! fluxes are in mmol m-2 d-1
           sfl(ppN3n) =   1.0*sfl_read(1)/SEC_PER_DAY
           sfl(ppN4n) =   1.0*sfl_read(2)/SEC_PER_DAY
           sfl(ppN1p) =0.0
!          sfl(ppN1p) =   1.0*sfl_read(3)/SEC_PER_DAY
        case (3) ! sfl array filled externally - for 3D models
!RUA
           sfl(ppN3n)= Nloss *0.12/0.21 * topm3psec
           sfl(ppN3n)= Nloss *0.09/0.21 * topm3psec
        case default
     end select
     ! assign surface fluxes to gotm surface flux array (sfl)
     ! convert from days to seconds
     ! oxygen flux is computed in WindOxReaeration
     sfl(ppO2o) =  PELSURFACE(ppO2o,1)/SEC_PER_DAY
     sfl(ppN3n) =  PELSURFACE(ppN3n,1)/SEC_PER_DAY
     sfl(ppN4n) =  PELSURFACE(ppN4n,1)/SEC_PER_DAY
     sfl(ppN1p) =  PELSURFACE(ppN1p,1)/SEC_PER_DAY
     sfl(ppN5s) =  PELSURFACE(ppN5s,1)/SEC_PER_DAY
   endif

   !---------------------------------------------
   ! Collect Bottom fluxes (mmol/m2/day)
   !---------------------------------------------
   if ( .NOT.AssignPelBenFluxesInBFMFlag ) then
     select case (bottom_flux_method)
        case (-1)! absolutely nothing
        case (0) ! default BFM with benthic model
           bfl(ppR6c) = PELBOTTOM(ppR6c,1)/SEC_PER_DAY
           bfl(ppR6n) = PELBOTTOM(ppR6n,1)/SEC_PER_DAY
           bfl(ppR6p) = PELBOTTOM(ppR6p,1)/SEC_PER_DAY
           bfl(ppR6s) = PELBOTTOM(ppR6s,1)/SEC_PER_DAY
     
           bfl(ppR1c) =  PELBOTTOM(ppR1c,1)/SEC_PER_DAY
           bfl(ppR1n) =  PELBOTTOM(ppR1n,1)/SEC_PER_DAY
           bfl(ppR1p) =  PELBOTTOM(ppR1p,1)/SEC_PER_DAY
     
           bfl(ppO2o) = PELBOTTOM(ppO2o,1)/SEC_PER_DAY
           bfl(ppN1p) = PELBOTTOM(ppN1p,1)/SEC_PER_DAY
           bfl(ppN3n) = PELBOTTOM(ppN3n,1)/SEC_PER_DAY
           bfl(ppN4n) = PELBOTTOM(ppN4n,1)/SEC_PER_DAY
           bfl(ppN5s) = PELBOTTOM(ppN5s,1)/SEC_PER_DAY
           bfl(ppN6r) = PELBOTTOM(ppN6r,1)/SEC_PER_DAY
     
           do i=1,iiPhytoPlankton
             k=ppPhytoPlankton(i,iiC) 
             bfl(k) = PELBOTTOM(k,1)/SEC_PER_DAY
             k=ppPhytoPlankton(i,iiN) 
             bfl(k) = PELBOTTOM(k,1)/SEC_PER_DAY
             k=ppPhytoPlankton(i,iiP) 
             bfl(k) = PELBOTTOM(k,1)/SEC_PER_DAY
             k=ppPhytoPlankton(i,iiL) 
             bfl(k) = PELBOTTOM(k,1)/SEC_PER_DAY
             k=ppPhytoPlankton(i,iiS)
             if ( k > 0 ) bfl(k) = PELBOTTOM(k,1)/SEC_PER_DAY
           enddo
        case (1) ! prescribed boundary fluxes (user)
           pelvar_bbc(ppN1p) = Dirichlet
           bfl(ppN1p) = 0.2_RLEN
           pelvar_bbc(ppN3n) = Dirichlet
           bfl(ppN3n) = 4.0_RLEN
           pelvar_bbc(ppN5s) = Dirichlet
           bfl(ppN5s) = 5.0_RLEN
           ! detritus flux at the bottom
           PELBOTTOM(ppR6c,1) = -p_rR6m*D3STATE(ppR6c,1)
           PELBOTTOM(ppR6n,1) = -p_rR6m*D3STATE(ppR6n,1)
           PELBOTTOM(ppR6p,1) = -p_rR6m*D3STATE(ppR6p,1)
           PELBOTTOM(ppR6s,1) = -p_rR6m*D3STATE(ppR6s,1)
           bfl(ppR6c) = PELBOTTOM(ppR6c,1)/SEC_PER_DAY
           bfl(ppR6n) = PELBOTTOM(ppR6n,1)/SEC_PER_DAY
           bfl(ppR6p) = PELBOTTOM(ppR6p,1)/SEC_PER_DAY
           bfl(ppR6s) = PELBOTTOM(ppR6s,1)/SEC_PER_DAY
      end select
   endif
   end subroutine do_bio_bfm
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  settling_vel_bfm
!
! !INTERFACE
   subroutine settling_vel_bfm
!
! !DESCRIPTION
!  This subroutine assign the sinking velocities of BFM variables
!
! !USES
   use mem, only: sediPI, sediR6, sediR2,iiC,iiN,iiP,iiS,iiL, &
                  sediMeZ,sediMiZ, &
                  ppR2c, ppR6c, ppR6n, ppR6p, ppR6s, NO_BOXES_Z,   &
                  ppR1c, ppR1n, ppR1p,   &
                  NO_D3_BOX_STATES, Depth,              &
                  ppPhytoPlankton,iiPhytoPlankton, &
                  ppMesoZooPlankton,iiMesoZooPlankton, &
                  ppMicroZooPlankton,iiMicroZooPlankton, &
                  PELBOTTOM, PELSURFACE
   use constants,  only: SEC_PER_DAY
   use gotm_error_msg, only:gotm_error

   IMPLICIT NONE
!
!  Local variables
   logical                     :: ll_larger
   integer                     :: n,k,i,j,l
   REALTYPE                    :: topm3psec,corr2
   REALTYPE                    :: corr(NO_BOXES_Z)        
   REALTYPE                    :: c1dimz(NO_BOXES_Z)        

   !---------------------------------------------
   ! Transfer sinking velocities (m/d -> m/s)
   !---------------------------------------------
   if ( bio_setup ==2 ) return
   do i=1,iiPhytoPlankton
     ll_larger=(maxval(sediPI(i,1:NO_BOXES_Z))> _ZERO_)
     l=ppPhytoPlankton(i,iiC)
     c1dimz(NO_BOXES_Z)= _ZERO_
     c1dimz(1:NO_BOXES_Z-1)=(Depth(2:NO_BOXES_Z)*sediPI(i,1:NO_BOXES_Z-1) &
                            +Depth(1:NO_BOXES_Z-1)*sediPI(i,2:NO_BOXES_Z))/ &
                                (Depth(1:NO_BOXES_Z)+Depth(2:NO_BOXES_Z-1))
     corr=min(cdepth(1)*rel_max_sedi_rate,abs(c1dimz))/(1.0D-80+abs(c1dimz))


     ws(l,1:NO_BOXES_Z) = -c1dimz(1:NO_BOXES_Z)/SEC_PER_DAY*corr
     llws(l)=ll_larger
     ws(l,0)= ws(l,1)

     k=ppPhytoPlankton(i,iiN)
     ws(k,0:NO_BOXES_Z) =ws(l,0:NO_BOXES_Z)
     llws(k)=ll_larger
     k=ppPhytoPlankton(i,iiP)
     ws(k,0:NO_BOXES_Z) = ws(l,0:NO_BOXES_Z)
     llws(k)=ll_larger
     k=ppPhytoPlankton(i,iiL)
     ws(k,0:NO_BOXES_Z) = ws(l,0:NO_BOXES_Z)
     llws(k)=ll_larger
     k=ppPhytoPlankton(i,iiS)
     if ( i==1  ) then 
          ws(k,0:NO_BOXES_Z) = ws(l,0:NO_BOXES_Z)
          llws(k)=ll_larger
     endif
   enddo
   do i=1,iiMesoZooPlankton
     ll_larger=(maxval(abs(sediMeZ(i,1:NO_BOXES_Z)))> 0.001)
     l=ppMesoZooPlankton(i,iiC)
     c1dimz(NO_BOXES_Z)=0.0;
     c1dimz(1:NO_BOXES_Z-1)=(Depth(2:NO_BOXES_Z)*sediMeZ(i,1:NO_BOXES_Z-1) &
                            +Depth(1:NO_BOXES_Z-1)*sediMeZ(i,2:NO_BOXES_Z))/ &
                                (Depth(1:NO_BOXES_Z)+Depth(2:NO_BOXES_Z-1))

     corr=min(cdepth(1)*rel_max_sedi_rate,abs(c1dimz))/(1.0D-80+abs(c1dimz))

     ws(l,1:NO_BOXES_Z) = -c1dimz(1:NO_BOXES_Z)/SEC_PER_DAY*corr
     llws(l)=ll_larger
     ws(l,0)= ws(l,1)

     k=ppMesoZooPlankton(i,iiN)
     if ( k>0 ) then
       ws(k,0:NO_BOXES_Z) =ws(l,0:NO_BOXES_Z)
       llws(k)=ll_larger
       k=ppMesoZooPlankton(i,iiP)
       ws(k,0:NO_BOXES_Z) = ws(l,0:NO_BOXES_Z)
       llws(k)=ll_larger
     endif
   enddo
   do i=1,iiMicroZooPlankton
     ll_larger=(maxval(abs(sediMiZ(i,1:NO_BOXES_Z)))> 0.001)
     l=ppMesoZooPlankton(i,iiC)
     c1dimz(NO_BOXES_Z)=0.0;
     c1dimz(1:NO_BOXES_Z-1)=(Depth(2:NO_BOXES_Z)*sediMiZ(i,1:NO_BOXES_Z-1) &
                            +Depth(1:NO_BOXES_Z-1)*sediMiZ(i,2:NO_BOXES_Z))/ &
                                (Depth(1:NO_BOXES_Z)+Depth(2:NO_BOXES_Z-1))
     corr=min(cdepth(1)*rel_max_sedi_rate,abs(c1dimz))/(1.0D-80+abs(c1dimz))


     ws(l,1:NO_BOXES_Z) = -c1dimz(1:NO_BOXES_Z)/SEC_PER_DAY*corr
     llws(l)=ll_larger
     ws(l,0)= ws(l,1)

     k=ppMicroZooPlankton(i,iiN)
     if ( k>0 ) then
       ws(k,0:NO_BOXES_Z) =ws(l,0:NO_BOXES_Z)
       llws(k)=ll_larger
       k=ppMesoZooPlankton(i,iiP)
       ws(k,0:NO_BOXES_Z) = ws(l,0:NO_BOXES_Z)
       llws(k)=ll_larger
     endif
   enddo

   c1dimz(NO_BOXES_Z)=0.0;
   c1dimz(1:NO_BOXES_Z-1)=(Depth(2:NO_BOXES_Z)*sediR2(1:NO_BOXES_Z-1) &
                            +Depth(1:NO_BOXES_Z-1)*sediR2(2:NO_BOXES_Z))/ &
                                (Depth(1:NO_BOXES_Z)+Depth(2:NO_BOXES_Z-1))
   corr=min(cdepth(1)*rel_max_sedi_rate,abs(c1dimz))/(1.0D-80+abs(c1dimz))
   ws(ppR2c,1:NO_BOXES_Z) = -c1dimz/SEC_PER_DAY*corr
   ws(ppR2c,0) = ws(ppR2c,1)
   llws(ppR2c)=.TRUE.

   ws(ppR6c,NO_BOXES_Z)=0.0
   corr=min(cdepth(1)*rel_max_sedi_rate,abs(sediR6))/(1.0D-80+abs(sediR6))
   c1dimz(1:NO_BOXES_Z-1)=(Depth(2:NO_BOXES_Z)*sediR6(1:NO_BOXES_Z-1) &
                            +Depth(1:NO_BOXES_Z-1)*sediR6(2:NO_BOXES_Z))/ &
                                (Depth(1:NO_BOXES_Z)+Depth(2:NO_BOXES_Z-1))
   ws(ppR6c,1:NO_BOXES_Z) = -c1dimz/SEC_PER_DAY*corr
   ws(ppR6c,0) = ws(ppR6c,1)

   ws(ppR6n,0:NO_BOXES_Z) =ws(ppR6c,0:NO_BOXES_Z)
   ws(ppR6p,0:NO_BOXES_Z) =ws(ppR6c,0:NO_BOXES_Z)
   ws(ppR6s,0:NO_BOXES_Z) =ws(ppR6c,0:NO_BOXES_Z)

   llws(ppR6c)=.TRUE.
   llws(ppR6n)=.TRUE.
   llws(ppR6p)=.TRUE.
   llws(ppR6s)=.TRUE.

   return

   end subroutine settling_vel_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Compute the vertical flux at a certain level
!
! !INTERFACE
   subroutine CalcVertFluxAtLev(statenr,lev,nlev,Depth,dt,out)
!
! !DESCRIPTION
! !USES
   IMPLICIT NONE
!
   integer,intent(IN)          :: statenr
   integer,intent(IN)          :: lev
   integer,intent(IN)          :: nlev
   REALTYPE,intent(IN)         :: Depth(1:nlev)
   REALTYPE,intent(IN)         :: dt
   REALTYPE,intent(OUT)        :: out

   out= sum((cc_before_transport(statenr,lev+1:nlev) &
               -cc(statenr,lev+1:nlev))*Depth(lev+1:nlev))/dt
   return
   end subroutine CalcVertFluxAtLev
!EOC

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE
   subroutine do_bfm_river_loads(action ,n,var)
!
! !DESCRIPTION
!  Get info on riverloads and move values to BFM array 
!  in order to make them available for output

! !USES
   use mem, only: PELRIVER
   use constants,  only: SEC_PER_DAY

! !INPUT PARAMETERS:
   IMPLICIT NONE
   integer,intent(IN)        ::action
   integer,intent(IN)        ::n
   REALTYPE,intent(IN)       ::var(1:n)
   
!EOP
!-----------------------------------------------------------------------
!BOC
     if ( action >=1 ) then 
       PELRIVER(1:n,1)=var(1:n) *SEC_PER_DAY
     else
       PELRIVER(1:n,1)=0.0D+00;
     endif
   end subroutine do_bfm_river_loads
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Reset diagonal id a 3d array
!
! !INTERFACE:
   subroutine reset_diagonal(n,pp)
!
! !DESCRIPTION:
!    Reset of the diagonal
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                    :: n
   REALTYPE,dimension(:,:,:),intent(inout) :: pp
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!
! !LOCAL VARIABLES:
   integer                   :: i
!EOP
!-----------------------------------------------------------------------
!BOC
     do i=1,n
       pp(i,i,:) = _ZERO_
     end do

   return
   end subroutine reset_diagonal
!EOC
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocate_bfm
!
! !INTERFACE:
        subroutine allocate_memory_bfm(nlev)
!
! !INPUT PARAMETERS:
        implicit none
        integer,intent(IN)            ::nlev
!
! !LOCAL VARAIBELS:
   integer                   :: rc
!
! !DESCRIPTION:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!  28-04-2006  Piet Ruardij Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

   if ( numc_diag > 0 ) then
     allocate(diag(1:numc_diag,0:nlev),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (cc)'
     diag=_ZERO_
   endif
   ! temporary array variable for the computation of vertical fluxes
   allocate(cc_before_transport(1:numc,0:nlev),stat=rc)
   if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (cc_before_transport)'
   cc_before_transport=_ZERO_

   if (bio_setup >= 2) then
     ! allocate benthic state variables
     allocate(ccb(1:numbc,1),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ccb)'
     allocate(ppb(1:numbc,1:numbc,1),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ppb)'
     allocate(ddb(1:numbc,1:numbc,1),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (ppb)'
     ccb=_ZERO_
     ppb=_ZERO_
     ddb=_ZERO_
     ! allocate variable holding type and save attributes
     allocate(benvar_type(1:numbc),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (benvar_type)'
     benvar_type = 0
   end if

   if ( numbc_diag > 0 ) then
     allocate(diagb(1:numbc_diag,1),stat=rc)
     if (rc /= 0) STOP 'allocate_memory_bfm: Error allocating (cc)'
     diagb=_ZERO_
   endif

 end subroutine allocate_memory_bfm
!EOC
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Test negative concentrations
!
! !INTERFACE:
       subroutine test_on_negative_states ( statenr,lldeep, c1dimz, &
                                            h, nlev, after, error )
!
! !DESCRIPTION:
!   Routine to check for negative values.
!   Negative values are corrected with the average of neighbour
!   grid points. A warning is given.
!
! !USES:
       use mem_Param, ONLY: p_small
       use gotm_error_msg, only:set_warning_for_getm
       IMPLICIT NONE
!
! !INPUT PARAMETERS:
       integer,intent(IN)                      :: statenr
       integer,intent(IN)                      :: nlev
       logical,intent(IN)                      :: lldeep
       REALTYPE,intent(IN)                     :: h(0:nlev)
       character(len=*),intent(IN)             :: after
!
! !OUTPUT PARAMETERS:
       integer,intent(OUT)                     :: error
!      Array cldim is modified if necessary
       REALTYPE,intent(INOUT)                  :: c1dimz(0:nlev)
!
! !LOCAL VARAIBELS:
        integer              ::k
        integer              ::i,n
        REALTYPE             ::r
        REALTYPE             ::sumbefore,sumafter
        character(len=160)   ::msg
        character(len=20)    ::onem
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!      Created by P. Ruardij 21-06-2006
!      MAV: Restored the argument c1dim for compatibility 
!           with the rest of the code
!
!EOP
!-------------------------------------------------------------------------
!BOC
       error=0
       ! in this way NaN are directly found!
       if (minval(c1dimz(1:nlev)) .ge. 0.00D+00) then           !BFM
          continue
       else
          if ( .not.lldeep ) then
              r=max(0.0D+00,sum(h(1:nlev)* c1dimz(1:nlev))/sum(h(1:nlev)))
              if ( r > 1.0D-10) then
                error=1
                write(msg,'(''statenr:'',I4,'' Negative value after call to '',A)') & 
                      statenr, after
                i=len_trim(msg)
                STDERR msg(1:i)
                STDERR "Averaging over the vertical value:",r
                call set_warning_for_getm()
             endif
             c1dimz(1:nlev)=r;
          else
            k=0                                                !BFM
            n=0
            sumbefore=sum(c1dimz(1:nlev)*h(1:nlev))
            do i = 1,nlev                                      !BFM
              if ( c1dimz(i).lt.0.0D+00) then                   !BFM
                  write(onem,'(I2,'':'',F10.3,''/'')') i,c1dimz(i) 
                  if (index(onem,'-')>0 )n=n+1
                  k=-i                                          !BFM
                  if ( i == 1 ) then
                     if ( c1dimz(i+1) > 0.0 )  then
                       c1dimz(i)=0.1* c1dimz(i+1)
                       k=i
                     else
                       c1dimz(i)=0.0
                       k=i
                     endif
                  elseif ( i == nlev ) then
                     if ( c1dimz(i-1) > 0.0 )  then
                       c1dimz(i)=0.1* c1dimz(i-1)
                       k=i
                     else
                       c1dimz(i)=0.0
                       k=i
                     endif
                  else if ( (c1dimz(i-1) > 0.0) .and. ( c1dimz(i+1)>0.0 ) ) then
                     k=i
                     c1dimz(i)=(c1dimz(i-1)+c1dimz(i+1)) * 0.1
                  else if ( c1dimz(i-1) >= 0.0 ) then
                       c1dimz(i)=0.1* c1dimz(i-1)
                       k=i
                  else if ( c1dimz(i+1) >= 0.0 ) then
                       c1dimz(i)=0.1* c1dimz(i+1)
                       k=i
                  endif
               endif
               if ( error.ge.0) error=k
            end do                                    !BFM
            sumafter=p_small+sum(c1dimz(1:nlev)*h(1:nlev))
            r=sumbefore/sumafter
            c1dimz=c1dimz*min(1.0,max(0.0,r))
            if ( n > 0 ) then
              call set_warning_for_getm()
              write(msg,'(''statenr='',I3,'' Negative values after '' &
                    ,A,'' n='',I2,'' shift='',F10.3,''%'')') &
                    statenr,after,n,100.0D+00*(1.0-r) 
              i=len_trim(msg)
              STDERR msg(1:i)
            endif
         endif
       endif

     end subroutine test_on_negative_states


!EOC
!-------------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_bfm
!
! !DESCRIPTION:
!  Nothing done here --- supplied for completeness
!  with GOTM bio structure.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_bfm
!EOC

!-----------------------------------------------------------------------


   end module bio_bfm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License
! www.gnu.org
!-----------------------------------------------------------------------
