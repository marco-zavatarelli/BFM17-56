#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   LimitShift.f90
!
! FILE
!   LimitChange.f90
!
! DESCRIPTION
!   
!	function FixProportionCoeff
!	function to calculate the limit  the shift
!	input: 
!  
! !INTERFACE
        subroutine FixProportionCoeff(NUTR,coeff1,coeff2,value1,value2)
!
! !AUTHORS
!   Piet Ruardij   
!
! !USES:
        USE global_mem,      ONLY:RLEN
        use constants, ONLY: PARAMETER, INPUT_TERM, START_ADD_TERM, INPUT_ADD_TERM 
        use bennut_interface, ONLY: CompleteSet
        use mem,ONLY:dummy
!
! CHANGE_LOG
!   
!
! COPYING
!   
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
        IMPLICIT  NONE
        INTEGER,intent(IN)        :: NUTR
        INTEGER,intent(IN)        :: coeff1
        INTEGER,intent(IN)        :: coeff2
        REAL(RLEN),intent(IN) ::value1
        REAL(RLEN),intent(IN) ::value2

        select case ( abs(value2)> 1.0D-20)
           case( .FALSE. )
              call CompleteSet( NUTR, INPUT_TERM, coeff2, PARAMETER, &
                                                          dummy, value=0.0D+00)
           case( .TRUE. )
             call CompleteSet( NUTR, START_ADD_TERM, coeff2, PARAMETER, &
                                                           dummy, mfac=1.0D+00/value2)
             call CompleteSet( NUTR, INPUT_ADD_TERM, coeff1, PARAMETER, &
                                                           dummy, mfac=-1.0/value1)
       end select
 
       return
      end

