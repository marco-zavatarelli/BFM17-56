#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Close forcing files in the standalone BFM
!
! !INTERFACE:
   subroutine end_envforcing_bfm
!
! !DESCRIPTION:
!  Read the namelist parameters for the analytical forcings
!  Also initialize additional components for other forcing methods
!
!
! !USES:
   use envforcing

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'end_envforcing_bfm'
    select case (forcing_method)
    case (1) ! analytical forcings
    case (2) ! input data
       LEVEL2 'Closing forcing data file:'
       LEVEL2 trim(forcing_file)
       close(unit_forcing)
       if (use_external_data) then
          LEVEL2 'Closing external data file:'
          LEVEL2 trim(data_file)
          close(unit_data)
       end if
    case (3) ! interactive air-sea fluxes
    end select
#ifdef INCLUDE_SEAICE
    LEVEL2 'Closing sea-ice forcing data file:'
    LEVEL2 trim(seaice_file)
    close(unit_seaice)
#endif

   return

   end subroutine end_envforcing_bfm
!EOC
