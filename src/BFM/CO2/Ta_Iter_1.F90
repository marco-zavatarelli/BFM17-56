#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE:  ta_iter_1
!
! !DESCRIPTION
! This routine expresses TA as a function of DIC, hSWS (H+ on
! sea water scale) and constants.
! It also calculates the derivative of this function with respect to
! hSWS. It is used in the iterative solution for hSWS. In the call
! "x" is the input value for hSWS, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhSWS
!
!   INTENT(IN) x = H+ total on Sea Water Scale, 
!   INTENT(OUT) fn = calculated value for TA
!   INTENT(OUT) df = calculated value for dTA/dHtotal
!
!       fn = hco3(x) + co3(x) + borate(x) + oh(x) + hpo4(x) +
!            2*po4(x) + silicate(x) + hfree(x) + hso4(x) +
!            hf(x) + h3po4(x) - ta
!
!       df = dfn/dx
!
! !INTERFACE
  subroutine ta_iter_1(x,fn,df)
!
! !USES
  use mem_CO2
!
! AUTHORS
!   the OCMIP team
!   Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
!   - by A. Mouchet, 2004:
!   Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
!    1) fixed errors in signs, parenthesis and coefficient c in derivative
!    2) changed from Total to Seawater Scale
!       * c defined for seawater H scale;
!       * fn and df adapted to KF on free H scale
!       * comments have been adapted
! 
!EOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
IMPLICIT NONE

    real(RLEN),intent(IN)  :: x
    real(RLEN),intent(OUT) :: fn,df
    real(RLEN)             :: x2,x3,k12,k12p,k123p,c,a,a2,da,b,b2,db  
    real(RLEN),parameter   :: T1=1.0_RLEN,T2=2.0_RLEN,T3=3.0_RLEN

        x2 = x*x
        x3 = x2*x
        k12 = s_k1*s_k2
        k12p = s_k1p*s_k2p
        k123p = k12p*s_k3p
        c = T1 + s_st/s_ks + s_ft/s_kf
        a = x3 + s_k1p*x2 + k12p*x + k123p
        a2 = a*a
        da = T3*x2 + T2*s_k1p*x + k12p
        b = x2 + s_k1*x + k12
        b2=b*b
        db = T2*x + s_k1

        fn = s_k1*x*s_dic/b +        &
             T2*s_dic*k12/b +        &
             s_bt/(T1 + x/s_kb) +    &
             s_kw/x +                &
             s_pt*k12p*x/a +         &
             T2*s_pt*k123p/a +       &
             s_sit/(T1 + x/s_ksi) -  &
             x/c -                   &
             s_st/(T1 + s_ks/(x/c))- &
             s_ft/(T1 + s_kf/(x/c))- &
             s_pt*x3/a -             &
             s_ta

        df = ((s_k1*s_dic*b) - s_k1*x*s_dic*db)/b2 -     &
             T2*s_dic*k12*db/b2 -                        &
             s_bt/s_kb/(T1+x/s_kb)**T2 -                 &
             s_kw/x2 +                                   &
             (s_pt*k12p*(a - x*da))/a2 -                 &
             T2*s_pt*k123p*da/a2 -                       &
             s_sit/s_ksi/(T1+x/s_ksi)**T2 -              &
             T1/c -                                      &
             s_st*(T1 + s_ks/(x/c))**(-T2)*(s_ks*c/x2) - &
             s_ft/(T1 + s_kf/(x/c))**(-T2)*(s_kf/x2) -   &
             s_pt*x2*(T3*a-x*da)/a2

  return

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  end subroutine ta_iter_1
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!EOC
