!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: Profq
!
!DESCRIPTION
!
!   Momme Butenschon, February 2005
!   Dipartimento di Fisica
!   Universita di Bologna
!
! This subroutine solves for the turbulent closure.
! Turbulent kinetic energy (Q2/2)
! Turbulent length scale (Q2l)
!
! !INTERFACE
      SUBROUTINE PROFQ(DT2)
!
! USES:

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Modules (use of ONLY is strongly encouraged!)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      use global_mem,ONLY: RLEN
      use POM,ONLY: h,a,c,kb,kq,dzz,dz,vh,vhp,wusurf,wvsurf,wubot,wvbot &
           ,q2f,s,t,q2lb,rho,dtef,sprod,km,u,v,bprod,prod,q2lf,z &
           ,l,sh,sm,kn,kh,gm,gh,zz,q2b,q2,UMOL

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
      REAL(RLEN) :: DT2
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Scalars
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      REAL(RLEN) :: A1,A2,B1,B2,C1,CIWC,COEF1,COEF2,COEF3,COEF4,COEF5,CONST1,DH, &
          E1,E2,E3,GEE,KAPPA,P,SMALL,SQ,zco,coef6
      INTEGER :: K,KI
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Arrays
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      REAL(RLEN) :: BOYGR(KB),CC(KB),TEMP1(KB),TEMP2(KB),TEMP3(KB)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Intrinsic Functions
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      INTRINSIC ABS,MIN,SQRT
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Data Statements
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      DATA A1,B1,A2,B2,C1/0.92,16.6,0.74,10.1,0.08/
      DATA E1/1.8/,E2/1.33/,E3/1.0/
      DATA KAPPA/0.40/,SQ/0.2/,CIWC/1.0/
      DATA GEE/9.806/
      DATA SMALL/1.E-8/
!      SM=KB*0.39
!      SH=KB*0.49
!      GM=KB*0.154
!      GH=KB*0.154
!     ..

      DH = H
      DO K = 2,KB - 1
          A(K) = -DT2* (KQ(K+1)+KQ(K)+2.*UMOL)*.5/ (DZZ(K-1)*DZ(K)*DH*DH)
          C(K) = -DT2* (KQ(K-1)+KQ(K)+2.*UMOL)*.5/ (DZZ(K-1)*DZ(K-1)*DH*DH)
      END DO
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  THE FOLLOWING SECTION SOLVES THE EQUATION
     !  DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      CONST1 = 16.6**.6666667*CIWC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Boundary Conditions
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      VH(1) = 0.0
      VHP(1) = SQRT(WUSURF**2+WVSURF**2)*CONST1
      Q2F(KB) = .5*SQRT((WUBOT+WUBOT)**2+ (WVBOT+WVBOT)**2)*CONST1
      DO K = 1,KB - 1

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Calculate pressure in units of decibars
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!          CC(K) = 1449.2 + 1.34* (S(K)-35.) + 4.55*T(K) - &
!                 0.045*T(K)**2 + 0.00821*P + (15.0**1.e-9*P**2)
!          TEMP1(K) = 2./CC(K)
!          TEMP2(K) = (0.00821*P)
!          TEMP3(K) = (1.-0.40* (P/CC(K)**2))

           CC(K) = CC(K)* (1.-TEMP1(K)* (TEMP2(K)+15.*1.e-9*P**2)* &
                  TEMP3(K))** (-0.5)
               p = -gee*1.025*ZZ(k)*dh*.1
               cc(k) = 1449.1 + .00821*p + 4.55*t(k) - .045*t(k)**2 + &
                          1.34*(S(k)-35.)
               cc(k) = cc(k) &
                          /sqrt((1.-.01642*p/cc(k))*(1.-0.40*p/cc &
                          (k)**2))
      END DO

      DO K = 2,KB - 1
          Q2B(K) = ABS(Q2B(K))
          Q2LB(K) = ABS(Q2LB(K))
          BOYGR(K) = GEE* (RHO(K-1)-RHO(K))/ (DZZ(K-1)*DH)! & (G)
!                     +GEE**2*2.*1.025/ (CC(K-1)**2+CC(K)**2) (G)
          DTEF(K) = Q2B(K)*SQRT(Q2B(K))/ (B1*Q2LB(K)+SMALL)
          SPROD(K) = .25*KM(K)* ((U(K)+U(K)-U(K-1)-U(K-1))**2+ &
                    (V(K)+V(K)-V(K-1)-V(K-1))**2)/ (DZZ(K-1)*DH)**2* &
                    CIWC**2
          BPROD(K) = KH(K)*BOYGR(K)
          PROD(K) = SPROD(K) + BPROD(K)
      END DO

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Sweep downward
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      DO K = 2,KB - 1
          VHP(K) = 1./ (A(K)+C(K)* (1.-VH(K-1))- (2.*DT2*DTEF(K)+1.))
          VH(K) = A(K)*VHP(K)
          VHP(K) = (-2.*DT2*PROD(K)+C(K)*VHP(K-1)-Q2B(K))*VHP(K)
      END DO

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Sweep upward
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      DO 104 K = 1,KB - 1
          KI = KB - K
          Q2F(KI) = VH(KI)*Q2F(KI+1) + VHP(KI)

  104 CONTINUE

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  The following section solves for the equation
     !  DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Boundary Conditions
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      VH(1) = 0.
      VHP(1) = 0.
      Q2LF(KB) = 0.

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Sweep downward
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      DO K = 2,KB - 1
          DTEF(K) = DTEF(K)* (1.+E2* ((1./ABS(Z(K)-Z(1))+1./ABS(Z(K)- &
                   Z(KB)))*L(K)/ (DH*KAPPA))**2)
          VHP(K) = 1./ (A(K)+C(K)* (1.-VH(K-1))- (DT2*DTEF(K)+1.))
          VH(K) = A(K)*VHP(K)
          VHP(K) = (DT2* (- (SPROD(K)+E3*BPROD(K))*L(K)*E1)+ &
                  C(K)*VHP(K-1)-Q2LB(K))*VHP(K)
      END DO

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Sweep upward
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      DO K = 1,KB - 1
          KI = KB - K
          Q2LF(KI) = VH(KI)*Q2LF(KI+1) + VHP(KI)
      END DO

      DO K = 2,KB - 1
          IF (Q2F(K).GT.SMALL .OR. Q2LF(K).GT.SMALL) GO TO 112
          Q2F(K) = SMALL
          Q2LF(K) = SMALL
  112 END DO

      DO k =1,kb-1
      q2f(k)=abs(q2f(k))
      q2lf(k)=abs(q2lf(k))
      END DO

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  The following section solves for KM and KH
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      COEF1 = A2* (1.-6.*A1/B1)
      COEF2 = 3.*A2*B2 + 18.*A1*A2
      COEF3 = A1* (1.-3.*C1-6.*A1/B1)
      COEF4 = 18.*A1*A1 + 9.*A1*A2
      COEF5 = 9.*A1*A2

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Note that SM and SH limit to infinity when GH approaches 0.0288
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      L(1) = 0.
      L(KB) = 0.
      GH(1) = 0.
      GH(KB) = 0.
      DO K = 2,KB - 1
          L(K) = Q2LF(K)/Q2F(K)
          GH(K) = L(K)**2/Q2F(K)*BOYGR(K)
      END DO
      DO K = 1,KB
          GH(K) = MIN(GH(K),.028)
          SH(K) = COEF1/ (1.-COEF2*GH(K))
          SM(K) = COEF3 + SH(K)*COEF4*GH(K)
          SM(K) = SM(K)/ (1.-COEF5*GH(K))
      END DO

      DO K = 1,KB
          KN(K) = L(K)*SQRT(ABS(Q2(K)))
          KQ(K) = (KN(K)*.41*SM(K)+KQ(K))*.5
!          KQ(K) = (KN(K)*.41*SH(K)+KQ(K))*.5
          KM(K) = (KN(K)*SM(K)+KM(K))*.5
          KH(K) = (KN(K)*SH(K)+KH(K))*.5
      END DO

      RETURN

      end subroutine PROFQ

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
