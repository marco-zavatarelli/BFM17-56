!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   calculatetau.f90
!
! FILE
!   calculatetau.f90
!
! DESCRIPTION
!   
!	function calculateTau
!	function to calculate the adaptation time
!	input: sMI             :first order process (/day)
!	xdiffMI         :diffusion constant (m2/day)
!	ptMI            :proportion between dissolved and dissolved in layer
!	0-DXm (-)
!	DXm             :underside of layer DXm (m)
!	local: pip2            : pi *pi 
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
!
!
      FUNCTION CalculateTau(sMI,xdiffMI,ptMI,DXm)
        USE global_mem, ONLY:RLEN,ONE
        IMPLICIT  NONE
        REAL(RLEN),intent(IN) ::smi ! Specification
        REAL(RLEN),intent(IN) ::xdiffmi ! Specification
        REAL(RLEN),intent(IN) ::ptmi ! Specification
        REAL(RLEN),intent(IN) ::dxm ! Specification
        REAL(RLEN) ::pip2
        REAL(RLEN) ::CalculateTau
        parameter (pip2=9.8696044D+00)

        CalculateTau= (ONE+ptMI)/(sMI+pip2*xdiffMI/(DXm**2))

        return
      end

