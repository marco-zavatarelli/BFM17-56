!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   CompleteSet
!   filly
!
! FILE
!   completeset.f90
!
! DESCRIPTION
!   
!  
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij. 
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
      SUBROUTINE CompleteSet(NUTR,mode,option,input,xinput,value,mfac)
        USE global_mem, ONLY:RLEN,ALLOC,error_msg_prn
        USE bennut_variables
        USE constants
        USE bennut_constants
        USE bennut_interface,ONLY:kfind, AddEquation, transfer
        IMPLICIT  NONE
        integer,intent(IN)            ::NUTR ! Specification ...dummy
        integer,intent(IN)            ::mode ! Specification
        integer,intent(IN)            ::input ! Specification
        integer,intent(IN)            ::option ! Specification
        REAL(RLEN),intent(IN)         ::xinput ! Specification
        REAL(RLEN),intent(IN),optional::value ! Specification
        REAL(RLEN),intent(IN),optional::mfac ! Specification

        integer ::i
        integer ::j
        integer ::k
        integer ::l
        integer ::m
        integer ::n
        integer ::option_local
        REAL(RLEN) ::yinput
        REAL(RLEN) ::multi
        REAL(RLEN) ::r
        REAL(RLEN) ::s
        REAL(RLEN) ::t
        REAL(RLEN) ::u
        logical ::control

        
        yinput=0.0D+00
        if (present(value)) yinput=value;
        multi=1.0D+00
        if (present(mfac)) multi=mfac;
        if (NUTR /= nutr_seq) stop 'error: wrong use of routine'

        select case (ns%status)
          case(DEFINE)
              nn_boundaries=0
              call DefineSet(NUTR,SET_BOUNDARY,0,0,0.0D+00,0.0D+00)
              C=0.0D+00
              Y2=0.0D+00
              ns%status=SET_BOUNDARY
              ModeMass=.false.
              ns%imethod=0
          case (SET_BOUNDARY,ADD) 
          case default 
             STOP 'equation only defined after finishing definition'
        end select 

        select case (mode)
          case ( ADD,SET_BOUNDARY)
            if (mode == SET_BOUNDARY) nn_boundaries=nn_boundaries+1
            s=multi
            if ( option > 0 ) &
              call AddEquation(nn_boundaries,option,input,ns,s,xinput)
             call filly(nn_boundaries, s* yinput ,Y2)
          case (INPUT_TERM,START_ADD_TERM,INPUT_ADD_TERM)
            !special continuity equations.........
            if (mode == INPUT_TERM.or. &
              mode == START_ADD_TERM) nn_boundaries=nn_boundaries+1
            s=multi;
            r=1.D+00
            t=yinput
            j=kfind(option,ns%coeffs,ns%nn)
            if (j > 100000) stop 'CompleteSet=40??'
            i=ns%coeffs(j)%il/10
            if (input == PARAMETER) then
              ! use COEFF->PARA because rX=Y is calculated
              ! instead of X=Y/r
              r=transfer(COEFF2PARA,ns%coeffs(j),r,ns%diff(i))
            elseif(input /= STANDARD) then
              stop 'CompleteSet:error mode=input_term'
            endif
            call calcadd(nn_boundaries,j,C,ns%nn,s*r)
            call filly(nn_boundaries,s*t,Y2)
          case (FLAG)
            if (option == MASS) ModeMass=.true.
          case (SET_CONTINUITY)
            if (option == FLAG) then
              if (input == MASS) ModeMass=.true.
            endif
            s=-1.D+00
            !For all lsts:
            if (ns%equa > 1) then
              do m=2,ns%equa
                !For the first derivative and the equation:
                l=1;if ( option ==LAYERS) then
                  if ( input.ne.m) l=0
                endif
                if ( l==1) then
                  do l=EQUATION,DERIVATIVE,DERIVATIVE
                    if (ns%lst(m-1) /= ns%lst(m)) then
                      nn_boundaries=nn_boundaries+1
                      !Make sum F1(x)=F2(X)
                      do k=m-1,m
                        t=1.D+00
                        if (l == DERIVATIVE ) t=ns%diff(k)*ns%poro(k)
                        s=-s
                        call AddEquation(nn_boundaries,k,l,ns,s*t,ns%b(m))
                      enddo
                    endif
                  enddo
                endif
              enddo
            endif
          case (SET_LAYER_INTEGRAL,SET_LAYER_INTEGRAL_UNTIL )
            option_local=option
            t=1.0D+00
            s=multi
            if ( option_local > 0 ) then
               nn_boundaries=nn_boundaries+1
            else
              option_local=-option
            endif
            do m=option_local,input
              if (ModeMass) t=ns%poro(m)*(ns%ads(m)+1.D+00)
              r=ns%b(m)
              u=ns%b(m+1)
              if (mode == SET_LAYER_INTEGRAL_UNTIL) then
                if (m == input) u=xinput
              endif
              do k=m,m+1
                s=-s
                call AddEquation(nn_boundaries,m,INTEGRAL,ns,s*t,r)
                r=u
              enddo
            enddo
            call filly(nn_boundaries,yinput,Y2)
        end select
        return
      end
 
!
      SUBROUTINE filly(nr,yinput,y)
        USE global_mem, ONLY:RLEN
        IMPLICIT  NONE
        integer,intent(IN) ::nr ! Specification
        REAL(RLEN),intent(IN)    ::yinput ! Specification
        REAL(RLEN),intent(INOUT) ::y(nr) ! Specification

        y(nr)= y(nr)+yinput

        return
      end

      SUBROUTINE calcadd(irow,k,C,nn,fc)
       USE global_mem, ONLY:RLEN
       IMPLICIT  NONE
        integer,intent(IN) ::irow ! Specification
        integer,intent(IN) ::nn ! Specification
        integer,intent(IN) ::k ! Specification
        REAL(RLEN),intent(INOUT) ::c(nn,nn) ! Specification
        REAL(RLEN),intent(IN) ::fc ! Specification

        C(irow,k)=C(irow,k)+fc

        return
      end
      SUBROUTINE AddEquation(nn_bound,layer,mode,nt,s, xinput)
        USE global_mem, ONLY:RLEN
        USE bennut_interface, ONLY:funcalc
        USE bennut_variables
        USE bennut_type
        IMPLICIT  NONE

        integer, intent(IN) :: nn_bound
        integer, intent(IN) :: layer
        integer, intent(IN) :: mode
        type (ty_set), intent(IN) :: nt
        real(RLEN), intent(IN) :: s
        real(RLEN), intent(IN) :: xinput

        integer ::j
        real(RLEN) ::r

        do j=nt%lst(layer),nt%lfi(layer)
           r=s*funcalc(mode,0,nt%coeffs(j),nt%b(layer),xinput)
           call calcadd(nn_bound,j,C,nt%nn,r)
        enddo

      end   
