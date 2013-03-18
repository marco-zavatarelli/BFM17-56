!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   LimitChange.f90
!
! FILE
!   LimitChange.f90
!
! DESCRIPTION
!    Function to limit the rate of change of a variable
!    according to its value. Uses a saturation function:
!                      M
!    limF = F * ---------------
!                  |F|/C + M
!  
! !INTERFACE
        subroutine LimitChange(mode,jK0x,K0x,max_change)
!
! !AUTHORS
!   Piet Ruardij   
!
! !USES:
        USE global_mem, ONLY:RLEN,ZERO
        USE mem_Param,  ONLY: p_small
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2007 P. Ruardij, the BFM team
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
        IMPLICIT  NONE
        INTEGER,intent(IN)        :: mode
        REAL(RLEN),intent(INOUT)  :: jK0x ! limited flux
        REAL(RLEN),intent(IN)     :: K0x  ! variable value
        REAL(RLEN),intent(IN)     :: max_change ! half-saturation constant
 
        if ( mode ==1 .or. ( mode==2.and.jK0x < ZERO) ) &
           jK0x=jK0x*max_change/(abs(jK0x)/(max(p_small,K0x))+max_change)

        return
      end subroutine LimitChange

