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
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
    REALTYPE,dimension(:),intent(in)        :: c
    REALTYPE,intent(in),optional            :: nc,pc,sc,lc
!
! !OUTPUT PARAMETERS:
    REALTYPE,dimension(:),intent(inout),optional :: n,p,s,l
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
!LOCAL VARIABLES:
    REALTYPE                     :: nc_ratio,pc_ratio,sc_ratio, &
                                    lc_ratio
!
!EOP
!-----------------------------------------------------------------------
!BOC

    nc_ratio = 0.0126 ! Redfield
    if ((present(nc)).AND.(nc>_ZERO_)) nc_ratio = nc

    pc_ratio = 0.7862e-3 ! Redfield
    if ((present(pc)).AND.(pc>_ZERO_)) pc_ratio = pc

    sc_ratio = 0.0145 ! Redfield
    if ((present(sc)).AND.(sc>_ZERO_)) sc_ratio = sc

    lc_ratio = 0.03 ! standard diatom value
    if ((present(lc)).AND.(lc>_ZERO_)) lc_ratio = lc

    if (present(n)) then
       where (n==_ZERO_) 
          n = nc_ratio*c
       end where
    end if
    if (present(p)) then
       where (p==_ZERO_) 
          p = pc_ratio*c
       end where
    end if
    if (present(s)) then
       where (s==_ZERO_) 
          s = sc_ratio*c
       end where
    end if
    if (present(l)) then
       where (l==_ZERO_) 
          l = sc_ratio*c
       end where
    end if

  end subroutine init_cnps
!EOC



!-----------------------------------------------------------------------

