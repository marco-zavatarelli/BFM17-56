!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   calculateset.f90
!
! FILE
!   calculateset.f90
!
! DESCRIPTION
!   
! AUTHORS
!   
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!
      FUNCTION CalculateSet(NUTR, mode,option,input, xinput,yinput)
        USE global_mem, ONLY:RLEN,ALLOC,error_msg_prn,ZERO
        USE bennut_variables, ONLY:ns,Y2,C,nutr_seq,nn_boundaries
        USE constants
        USE bennut_constants
        USE bennut_interface,ONLY:CompleteSet, CalculateFromCondition
        IMPLICIT  NONE
        integer,intent(IN) ::NUTR ! Specification
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::option ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification
        REAL(RLEN)            :: CalculateSet

        real(RLEN),dimension(NCOEFF)  :: Y3
        real(RLEN),dimension(NCOEFF*NCOEFF):: C2
        real(RLEN),dimension(NCOEFF*NCOEFF):: C3
        integer,dimension(NCOEFF) ::iindx
        integer, external  :: PrintSet
   
        REAL(RLEN) ::xh
        INTEGER    ::ok

        if (NUTR /= nutr_seq) stop 'error: wrong use of routine'

        !calculate coefficients.......
        if (ns%nn == 0.or.nn_boundaries == 0) then
          write(0,'(''ns%nn='',i4,'' equations='',i4)') &
                                         ns%nn,nn_boundaries
          stop 'ns%nn=0 or nn_boundaries equation=0'
        elseif (ns%nn.ne.nn_boundaries) then
          write(0,'(''termsn='',i4,'' equations='',i4)')    &
                                                 ns%nn,nn_boundaries
        endif

        if ( mode == ADD ) then
          call CompleteSet(NUTR,ADD,0,0,ZERO,value=yinput)
          ns%factor(1:nn_boundaries)=Y2(1:nn_boundaries)
          if (ns%imethod == 0) then
            call ludcmp(nn_boundaries,C,iindx,xh,ok)
            if ( ok.eq.0 ) goto 100
            call lubksb(nn_boundaries,C,iindx,ns%factor)
          else
            call svdcmp(C,nn_boundaries,nn_boundaries,ns%nn,   &
                                                 ns%nn,Y3,C3)
            call set_max_sing(NUTR,Y3,ns%nn)
            call svbksb(C,Y3,C3,nn_boundaries,nn_boundaries,   &
                                          ns%nn,ns%nn,Y2,ns%factor)
          endif
          ns%status=READY
          CalculateSet= yinput
        else 
          call re_Store(nn_boundaries,C,C2,ns%nn,ns%nn)
          ns%factor(1:nn_boundaries)=Y2(1:nn_boundaries)
          if (ns%imethod == 0) then
            call ludcmp(nn_boundaries,C2,iindx,xh,ok)
            if ( ok.eq.0 ) goto 100
            call lubksb(nn_boundaries,C2,iindx,ns%factor)
          else
            call svdcmp(C2,nn_boundaries,nn_boundaries,ns%nn,ns%nn,Y3,C3)
            call set_max_sing(NUTR,Y3,ns%nn)
            call svbksb(C2,Y3,C3,nn_boundaries,nn_boundaries,   &
                                           ns%nn,ns%nn,Y2,ns%factor)
          endif
          if ( mode ==0 ) then
            ns%status=READY
            CalculateSet=ZERO
          else
            !Set n'th row on zero:
            Y2(nn_boundaries:nn_boundaries)=ZERO
            nn_boundaries=ns%nn-1
            call re_Store(-nn_boundaries,C,C,ns%nn,ns%nn)
            call CompleteSet(NUTR,mode,option,input,xinput,value=yinput)
            xh=CalculateFromCondition(ns%nn,C,ns%factor,ns%nn)
            xh=max(yinput,xh)
            CalculateSet=xh
            ns%status=ADD
          endif
          return
        endif
        return
  100   ok= PrintSet( NUTR,"Singular matrix in ludcmp")
        stop
      end

      FUNCTION CalculateFromCondition (irow,C,Y,nn)
        USE global_mem, ONLY: RLEN
        IMPLICIT  NONE
        integer,intent(IN) ::irow ! Specification
        integer,intent(IN) ::nn ! Specification
        REAL(RLEN),intent(IN) ::c(nn,nn) ! Specification
        REAL(RLEN),intent(IN) ::y(nn) ! Specification
        REAL(RLEN)            :: CalculateFromCondition
 
        CalculateFromCondition=dot_product(C(irow,1:nn),Y(1:nn))
        return
      end
