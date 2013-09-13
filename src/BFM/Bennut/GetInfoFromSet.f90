!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   GetInfoFromSet
!   calculate_one_term
!
! FILE
!   GetInfoFromSet.f90
!
! DESCRIPTION
!
! AUTHORS
! CHANGE_LOG
! COPYING
!
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2004 P. Ruardij, M. Vichi
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
      FUNCTION GetInfoFromSet(NUTR,option,input,termnr,at_x,to_x)
        USE global_mem, ONLY:ONE
        USE constants, ONLY:RLEN,GET,COEFFICIENT,LABDA_1,LABDA_2, &
            INTEGRAL,EXPONENTIAL_INTEGRAL,RFLUX,MASS,DERIVATIVE,PARAMETER,COEFF2PARA
        USE bennut_interface,ONLY: kfind, transfer,funcalc
        USE bennut_variables, ONLY: sets
        IMPLICIT  NONE
        integer,intent(IN)             ::nutr ! Specification
        integer,intent(IN)             ::termnr ! Specification
        integer,intent(IN)             ::option ! Specification
        integer,intent(IN)             ::input ! Specification
        REAL(RLEN),intent(IN),optional ::at_x ! Specification
        REAL(RLEN),intent(IN),optional ::to_x ! Specification
        REAL(RLEN)                     ::GetInfoFromSet

        integer      ::seqnr,j,layer
        real(RLEN)   ::bC,r,s

        seqnr=kfind(termnr,sets(NUTR)%coeffs,sets(NUTR)%nn)
        if (seqnr >1000) stop 'GetInfoFromSet mode=GET??'
        if( option == GET ) then
          select case (input)
            case (LABDA_1,LABDA_2)
              j=input-LABDA_1+1
              GetInfoFromSet=sets(NUTR)%coeffs(seqnr)%labda(j)
            case (COEFFICIENT)
              if (seqnr.gt.sets(NUTR)%nn) then
                GetInfoFromSet=-1.E31_RLEN
              else
                GetInfoFromSet=sets(NUTR)%factor(seqnr)
              endif
            case (COEFF2PARA)
               j=sets(NUTR)%coeffs(seqnr)%il/10
               GetInfoFromSet =transfer(COEFF2PARA, &
                    sets(NUTR)%coeffs(seqnr),ONE,sets(NUTR)%diff(j))
            case default
              stop 'GetInfoFromSet mode=????'
          end select
        elseif (option < GET) then
          if ( .not.present(at_x)) then
             stop 'GetInfoFromSet optional at_x NOT defined'
          endif
          layer=termnr/10
          bC=sets(NUTR)%b(layer)
          s=sets(NUTR)%factor(seqnr)
          if (input == PARAMETER)  then
            j = -2
          else
            j = 0
          end if
          if ( s/=0.0 ) then
            r= funcalc(option,j,sets(NUTR)%coeffs(seqnr),bC,at_x)
            !calculate the result at the upper border and subtract it from the
            ! lower layer
            if (option == INTEGRAL.or. &
              option == EXPONENTIAL_INTEGRAL) then
              if ( .not.present(to_x)) then
               stop 'GetInfoFromSet optional to_x NOT defined'
              endif
              r= funcalc(option,j,sets(NUTR)%coeffs(seqnr),bC,to_x)-r
              if (input == RFLUX.or.input == MASS) then
                r=r*sets(NUTR)%poro(layer)
                if(input == MASS)r=r*(sets(NUTR)%ads(layer)+1.D+00)
              endif
            elseif (option == DERIVATIVE) then
              if (input == RFLUX) &
              r=r*sets(NUTR)%diff(layer)*sets(NUTR)%poro(layer)
            endif
            if (input == PARAMETER) then
              r=transfer(COEFF2PARA,sets(NUTR)%coeffs(seqnr),r,  &
                         sets(NUTR)%diff(layer)) *sets(NUTR)%poro(layer)
            endif
          endif
          GetInfoFromSet=s*r
        endif

        return
      end function GetInfoFromSet
