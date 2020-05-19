!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: Profe
!
!DESCRIPTION    
!   This subroutine solves for vertical diffusivity
!
!   Momme Butenschon, May 2004
!   Dipartimento di Fisica
!   Universita di Bologna
!
! !INTERFACE
      SUBROUTINE PROFE(FF,WFSURF,FSURF,NBC,DT2)
!
! USES:

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Modules (use of ONLY is strongly encouraged!)
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      use global_mem,ONLY: RLEN
      use POM,ONLY: H,KB,A,C,L,DZ,DZZ,VH,VHP,Z,UMOL,KH
!
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
      real(RLEN) :: DT2,FSURF,WFSURF
      integer :: NBC

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Array Arguments 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      real(RLEN) :: FF(KB)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Scalars
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      real(RLEN) :: DH,UMOLPR
      integer :: K,KI
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Arrays
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      real(RLEN) :: EXTC(5),RAD(KB),TR(5)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Intrinsic Functions
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      intrinsic EXP
!     
      UMOLPR = 1.e-5
!
      DH = H
      DO K = 2,KB - 1
          A(K-1) = -DT2* (KH(K)+UMOLPR)/ (DZ(K-1)*DZZ(K-1)*DH*DH)
          C(K) = -DT2* (KH(K)+UMOLPR)/ (DZ(K)*DZZ(K-1)*DH*DH)
      END DO

      GO TO (51,52) NBC

   51 CONTINUE
      VH(1) = A(1)/ (A(1)-1.)
      VHP(1) = -DT2* WFSURF/ (-DZ(1)*DH) - FF(1)
      VHP(1) = VHP(1)/ (A(1)-1.)
      GO TO 53
   52 CONTINUE
      VH(1) = 0.
      VHP(1) = FSURF
   53 CONTINUE

     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! THE FOLLOWING SECTION SOLVES THE EQUATION                     
     ! DT2*(KH*FF')' -FF = FB                                        
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      DO K = 2,KB - 2
          VHP(K) = 1./ (A(K)+C(K)* (1.-VH(K-1))-1.)
          VH(K) = A(K)*VHP(K)
          VHP(K) = (C(K)*VHP(K-1)-FF(K))*VHP(K)
      END DO
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! INSTEAD OF MATCHING SOLUTION TO A LOWER LAYER VALUE OF F(KB-1),  
     ! ONE MAY IMPOSE AN ADIABATIC BOTTOM B.C. BY REMOVING C'S FROM     
     ! COL. 1 OF THE NEXT TWO LINES.                                    
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      FF(KB-1) = (C(KB-1)*VHP(KB-2)-FF(KB-1))/ &
                (C(KB-1)* (1.-VH(KB-2))-1.)
   99 CONTINUE
      DO K = 2,KB - 1
          KI = KB - K
          FF(KI) = VH(KI)*FF(KI+1) + VHP(KI)
      END DO
      DO K = 1,KB
          VH(K) = 0.0
          VHP(K) = 0.0
          A(K) = 0.0
          C(K) = 0.0
      END DO
      RETURN
      end subroutine profe
!
! EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
