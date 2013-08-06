#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ClearMem.F90
!
! !INTERFACE
subroutine ClearMem
!
! !DESCRIPTION
!  Deallocate all the arrays
!
! !USES
   use mem
   use api_bfm
   implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV)
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
!

   integer :: i,j,origin,destination

    ! free willy, free the fluxes
    origin=0
    do i=stPelStateS,stPelStateE
       origin=origin+1
       destination=0
       do j=stPelStateS,stPelStateE
          destination=destination+1
          if( allocated(D3FLUX_MATRIX(origin,destination)%p) ) deallocate(D3FLUX_MATRIX(origin,destination)%p)
       end do
    end do
    deallocate(D3FLUX_MATRIX)
    deallocate(D3FLUX_FUNC)
#if defined INCLUDE_SEAICE
    origin=0
    do i=stIceState2dS,stIceState2dE
       origin=origin+1
       destination=0
       do j=stIceState2dS,stIceState2dE
          destination=destination+1
          if( allocated(D2FLUX_MATRIX_ICE(origin,destination)%p) ) deallocate(D2FLUX_MATRIX_ICE(origin,destination)%p)
       end do
    end do
    deallocate(D2FLUX_MATRIX_ICE)
    deallocate(D2FLUX_FUNC_ICE)
#endif
#if defined INCLUDE_BEN
    origin=0
    do i=stBenState2dS,stBenState2dE
       origin=origin+1
       destination=0
       do j=stBenState2dS,stBenState2dE
          destination=destination+1
          if( allocated(D2FLUX_MATRIX_BEN(origin,destination)%p) ) deallocate(D2FLUX_MATRIX_BEN(origin,destination)%p)
       end do
    end do
    deallocate(D2FLUX_MATRIX_BEN)
    deallocate(D2FLUX_FUNC_BEN)
#endif

    ! from api_bfm 
    deallocate(var_ids)
    deallocate(var_ave)
    deallocate(var_names)
    deallocate(var_units)
    deallocate(var_long)
    deallocate(c1dim)
    if (allocated(D3ave)) deallocate(D3ave)
    if (allocated(D2ave)) deallocate(D2ave)

#ifndef NOT_STANDALONE
     deallocate(D3STATE)
     deallocate(D3SOURCE)
#ifdef EXPLICIT_SINK
     deallocate(D3SINK)
#if defined INCLUDE_SEAICE
     deallocate(D2SINK_ICE)
#endif
#if defined INCLUDE_BEN
     deallocate(D2SINK_BEN)
#endif
#endif

     deallocate(D3STATETYPE)
     deallocate(D3DIAGNOS)
     deallocate(D2DIAGNOS)
#ifdef BFM_NEMO
     deallocate(D3STATEOBC)
#endif

#if defined INCLUDE_SEAICE
     if ( allocated(D2ave_ice) ) deallocate(D2ave_ice)
     deallocate(D2DIAGNOS_ICE)
     deallocate(D2STATE_ICE)
     deallocate(D2SOURCE_ICE)
     deallocate(D2STATETYPE_ICE)
#ifdef BFM_NEMO
     deallocate(D2STATEOBC_ICE)
#endif
#endif

#if defined INCLUDE_BEN
     if ( allocated(D2ave_ben) ) deallocate(D2ave_ben)
     deallocate(D2DIAGNOS_BEN)
     deallocate(D2STATE_BEN)
     deallocate(D2SOURCE_BEN)
     deallocate(D2STATETYPE_BEN)
#ifdef BFM_NEMO
     deallocate(D2STATEOBC_BEN)
#endif
#endif

#endif

end subroutine ClearMem
!EOC
!-----------------------------------------------------------------------
