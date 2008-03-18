#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !SUBROUTINE: Initialise variable components
!
! !INTERFACE:
   subroutine init_cnps(c,n,p,s,l,nc,pc,sc,lc)
!
! !DESCRIPTION:
!  This subroutine initialises the other internal components
!  of biogeochemical variables
!
! !USES:
   use global_mem, only: RLEN,ZERO
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
    real(RLEN),dimension(:),intent(in)        :: c
    real(RLEN),intent(in),optional            :: nc,pc,sc,lc
!
! !OUTPUT PARAMETERS:
    real(RLEN),dimension(:),intent(inout),optional :: n,p,s,l
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
!LOCAL VARIABLES:
    real(RLEN)                   :: nc_ratio,pc_ratio,sc_ratio, &
                                    lc_ratio
!
!EOP
!-----------------------------------------------------------------------
!BOC

    nc_ratio = 0.0126_RLEN ! Redfield
    if ((present(nc)).AND.(nc>ZERO)) nc_ratio = nc

    pc_ratio = 0.7862e-3_RLEN ! Redfield
    if ((present(pc)).AND.(pc>ZERO)) pc_ratio = pc

    sc_ratio = 0.0145_RLEN ! Redfield
    if ((present(sc)).AND.(sc>ZERO)) sc_ratio = sc

    lc_ratio = 0.03_RLEN ! standard diatom value
    if ((present(lc)).AND.(lc>ZERO)) lc_ratio = lc

    if (present(n)) then
       where (n==ZERO) 
          n = nc_ratio*c
       end where
    end if
    if (present(p)) then
       where (p==ZERO) 
          p = pc_ratio*c
       end where
    end if
    if (present(s)) then
       where (s==ZERO) 
          s = sc_ratio*c
       end where
    end if
    if (present(l)) then
       where (l==ZERO) 
          l = sc_ratio*c
       end where
    end if

  end subroutine init_cnps
!EOC



!-----------------------------------------------------------------------

