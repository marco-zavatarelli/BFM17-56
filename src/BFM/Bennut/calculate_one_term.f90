!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   calculate_one_term.f90
!
! FILE
!   calculate_one_term.f90
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
      FUNCTION calculate_one_term(mode,option,xinput,coeff,b,factor)
        use global_mem, ONLY:RLEN,ONE,ZERO
        use bennut_type
        use constants
        use bennut_interface,ONLY:funcalc
        IMPLICIT  NONE
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::option ! Specification
        type (ty_coeff),intent(IN) ::coeff ! Specification
        real(RLEN),intent(IN) ::xinput ! Specification
        real(RLEN),intent(IN) ::b ! Specification
        real(RLEN),intent(IN) ::factor ! Specification
        real(RLEN)            ::calculate_one_term

        real(RLEN) ::r 
        integer ::i

        r=factor;
        if (r /= ZERO) then
          if (mode == PARAMETER) then
             i = -2
          else
             i = 0
          end if
          r=r*funcalc(option,i,coeff,b,xinput)
        endif
        calculate_one_term=r

        return
      end

