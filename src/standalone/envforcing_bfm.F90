#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light and other environmental forcings used in the BFM
!
! !INTERFACE
   subroutine envforcing_bfm(step)
!
! !DESCRIPTION
! This routine sets the environmental forcings according to user
! choice. The currently implemented methods are:
! 1) analytical sinusoidal forcings
! 2) external values from files
! 3) interactive slab ocean from meteorological data
! Any combination of the method above must be provided by the user
! in an additional subroutine.
!
! !USES
   use global_mem, only: RLEN
   use envforcing
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV), Jordi Sole (IMEDEA)
!
! !LOCAL VARIABLES:
   integer        :: step

!EOP
!-----------------------------------------------------------------------
!BOC
    select case (forcing_method)
    case (1) ! analytical forcings
      call analytical_forcing
    case (2) ! input data
      call analytical_forcing
      call external_forcing
    case (3) ! interactive air-sea fluxes
!      call do_air_sea(timesec,startime)
    end select
    ! Assign external data
    call external_data
    ! Assign external event data
    call event_data
#ifdef INCLUDE_SEAICE
    call external_seaice
#endif
    if (init_forcing_vars) init_forcing_vars=.false.
  end subroutine envforcing_bfm
!EOC
!-----------------------------------------------------------------------

