!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: Calcdepth
!
!DESCRIPTION    
!
!  This subroutine establishes the vertical resolution log distributions
!  at the top and bottom, and a linear distribution between KL1 and KL2.
!  Default values: KL1 = .3*KB AND KL2 = KB-2.
!  Yields a log distribution at the top and none at the bottom.
!
! !INTERFACE
      Subroutine CALCDEPTH(Z,ZZ,DZ,DZZ,KB,KL1,KL2)
!
! USES:

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      use global_mem,ONLY: RLEN

!-------------------------------------------------------------------------!
!BOC
!
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Implicit typing is never allowed
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      IMPLICIT NONE
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Scalar Arguments
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      INTEGER :: KB
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Array Arguments 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      REAL(RLEN) :: DZ(KB),DZZ(KB),Z(KB),ZZ(KB)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Scalars
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      REAL(RLEN) :: BB,CC,DEL1,DEL2
      INTEGER :: K,KL1,KL2
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Intrinsic Functions
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      INTRINSIC EXP,FLOAT

      BB = FLOAT(KL2-KL1) + 4.
      CC = FLOAT(KL1) - 2.
      DEL1 = 2./BB/EXP(.693147*FLOAT(KL1-2))
      DEL2 = 2./BB/EXP(.693147*FLOAT(KB-KL2-1))
      Z(1) = 0.
      ZZ(1) = -DEL1/2.

      DO 2 K = 2,KL1 - 2
          Z(K) = -DEL1*EXP(.693147*FLOAT(K-2))
    2 ZZ(K) = -DEL1*EXP(.693147* (FLOAT(K)-1.5))

      DO 4 K = KL1 - 1,KL2 + 1
          Z(K) = - (FLOAT(K)-CC)/BB
    4 ZZ(K) = - (FLOAT(K)-CC+0.5)/BB

!      DO 5 K = KL2 + 1,KB - 1
!          Z(K) = (1.0-DEL2*EXP(.693147*FLOAT(KB-K-1)))* (-1.)
!    5 ZZ(K) = (1.0-DEL2*EXP(.693147* (FLOAT(KB-K)-1.5)))* (-1.)
!      Z(KB) = -1.0
!      ZZ(KB-1) = -1.* (1.-DEL2/2.)
!      ZZ(KB) = -1.* (1.+DEL2/2.)

      DO K = 1,KB - 1
          DZ(K) = Z(K) - Z(K+1)
          DZZ(K) = ZZ(K) - ZZ(K+1)

      END DO
      RETURN
      
      end subroutine CALCDEPTH
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

