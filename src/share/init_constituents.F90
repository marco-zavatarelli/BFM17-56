#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !SUBROUTINE: Initialise internal constituents given carbon
!
! !INTERFACE:
   subroutine init_constituents(c,n,p,s,l,f,nc,pc,sc,lc,fc)
!
! !DESCRIPTION:
!  This subroutine initialises the internal constituents
!  of biogeochemical variables given carbon content
!
! !USES:
   use global_mem, only: RLEN,ZERO
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
    real(RLEN),dimension(:),intent(in)        :: c
    real(RLEN),intent(in),optional            :: nc,pc,sc,lc,fc
!
! !OUTPUT PARAMETERS:
    real(RLEN),dimension(:),intent(inout),optional :: n,p,s,l,f
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
!LOCAL VARIABLES:
    real(RLEN)                   :: nc_ratio,pc_ratio,sc_ratio, &
                                    lc_ratio,fc_ratio
!
! COPYING
!
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!EOP
!-----------------------------------------------------------------------
!BOC

    nc_ratio = 0.0126_RLEN ! Redfield
    if (present(nc)) then
       if (nc>ZERO) nc_ratio = nc
    end if

    pc_ratio = 0.7862e-3_RLEN ! Redfield
    if (present(pc)) then
       if (pc>ZERO) pc_ratio = pc
    end if

    sc_ratio = 0.0145_RLEN ! Redfield
    if (present(sc)) then 
       if (sc>ZERO) sc_ratio = sc
    end if

    lc_ratio = 0.03_RLEN ! standard diatom value
    if (present(lc)) then
       if (lc>ZERO) lc_ratio = lc
    end if

    fc_ratio = 3.e-04_RLEN ! standard diatom value
    if (present(fc)) then
       if (fc>ZERO) fc_ratio = fc
    end if

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
          l = lc_ratio*c
       end where
    end if
    if (present(f)) then
       where (f==ZERO) 
          f = fc_ratio*c
       end where
    end if

  end subroutine init_constituents
!EOC
!-----------------------------------------------------------------------

