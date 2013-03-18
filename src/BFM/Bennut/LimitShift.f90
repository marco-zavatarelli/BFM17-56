!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   LimitShift.f90
!
! FILE
!   LimitShift.f90
!
! DESCRIPTION
!	function to calculate the limit of  the layer shift
!	input: 
!  
! !INTERFACE
        subroutine LimitShift(jK10K0x,K0x,K10x,max_shift)
!
! !AUTHORS
!   Piet Ruardij   
!
! !USES:
        USE global_mem,      ONLY:RLEN,ZERO
        use mem_globalfun,   ONLY: insw
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
        IMPLICIT  NONE
        REAL(RLEN),intent(INOUT) ::jK10K0x ! Specification
        REAL(RLEN),intent(IN) ::K0x ! Specification
        REAL(RLEN),intent(IN) ::K10x ! Specification
        REAL(RLEN),intent(IN) ::max_shift ! Specification
        REAL(RLEN)            :: r 
 
       
        r= 1.0D-80+max(ZERO,insw(jK10K0x)* K10x+insw(-jK10K0x)* K0x)
        jK10K0x=jK10K0x*max_shift/(abs(jK10K0x/r)+max_shift);

        return
      end

