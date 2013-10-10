#include "cppdefs.h"
#include "INCLUDE.h"
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !ROUTINE: ResetFluxes.F90
  !
  ! !INTERFACE
  subroutine ResetFluxes
    !
    ! !DESCRIPTION
    !  Reset the arrays for the next integration
    !
    ! !USES
    use global_mem, only:ZERO,LOGUNIT
#ifdef NOPOINTERS
    use mem
#else

    use mem, ONLY: NO_D3_BOX_STATES,D3SOURCE, &
         PELBOTTOM, PELSURFACE, D3FLUX_FUNC
#ifdef EXPLICIT_SINK
    use mem, ONLY: D3SINK
#endif

#if defined INCLUDE_SEAICE
    use mem, ONLY: NO_D2_BOX_STATES_ICE, D2SOURCE_ICE, &
         D2FLUX_FUNC_ICE
#ifdef EXPLICIT_SINK
    use mem, ONLY: D2SINK_ICE
#endif
#endif


#if defined INCLUDE_BEN
    use mem, ONLY: NO_D2_BOX_STATES_BEN, D2SOURCE_BEN, &
         D2FLUX_FUNC_BEN
#ifdef EXPLICIT_SINK
    use mem, ONLY: D2SINK_BEN
#endif
#endif

#endif
    implicit none
    !
    ! !INPUT PARAMETERS:
    !
    ! !OUTPUT PARAMETERS:
    !
    ! !REVISION HISTORY:
    !  Original author(s): Momme Butenschoen (UNIBO)
    !
    ! !LOCAL VARIABLES:
    integer :: i
    !EOP
    !-----------------------------------------------------------------------
    !BOC
    !

    ! Reset the source and sink arrays
#ifndef EXPLICIT_SINK

    D3SOURCE(:,:) = ZERO
#if defined INCLUDE_SEAICE
    D2SOURCE_ICE(:,:) = ZERO
#endif
#if defined INCLUDE_BEN
    D2SOURCE_BEN(:,:) = ZERO
#endif

#else

    D3SOURCE(:,:,:) = ZERO
    ! Reset sink term arrays 
    D3SINK(:,:,:) = ZERO
#if defined INCLUDE_SEAICE
    D2SOURCE_ICE(:,:,:) = ZERO
    D2SINK_ICE(:,:,:) = ZERO
#endif
#if defined INCLUDE_BEN
    D2SOURCE_BEN(:,:,:) = ZERO
    D2SINK_BEN(:,:,:) = ZERO
#endif

#endif

    if (allocated(D3FLUX_FUNC)) D3FLUX_FUNC(:,:) = ZERO
#if defined INCLUDE_SEAICE
    if (allocated(D2FLUX_FUNC_ICE)) D2FLUX_FUNC_ICE(:) = ZERO
#endif
#if defined INCLUDE_BEN
    if (allocated(D2FLUX_FUNC_BEN)) D2FLUX_FUNC_BEN(:) = ZERO
#endif

    ! reset surface and bottom fluxes
    do i=1,NO_D3_BOX_STATES
       PELSURFACE(i,:) = ZERO 
       PELBOTTOM(i,:) = ZERO
    end do

  end subroutine ResetFluxes
  !EOC
  !-----------------------------------------------------------------------
