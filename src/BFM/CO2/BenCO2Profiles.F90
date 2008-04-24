#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   BenCO2Profiles.f90
!
! FILE
!   LimitChange.f90
!
! DESCRIPTION
!   
!	function BenCO2Profiles
!  
! !INTERFACE
      subroutine BenCO2Profiles
!
! !AUTHORS
!   Piet Ruardij   
!
! !USES:
        USE global_mem,      ONLY:RLEN,LOGUNIT,ZERO
#ifdef INCLUDE_BENPROFILES
#ifdef INCLUDE_BENCO2
        USE mem, ONLY:BoxNumberZ, NO_BOXES_Z, BoxNumberX, NO_BOXES_X, &
           BoxNumberY,NO_BOXES_Y,BoxNumber,BoxNumberXY, seddepth
        use mem,ONLY:PrDIC,PrAc,KALK,KCO2,PrM1p,PrM5s,PrpH, &
           ESW_Ben,ETW_Ben,ERHO_Ben
        USE constants, ONLY: INTEGRAL,STANDARD,MW_C
        use bennut_interface, ONLY:CalculateFromSet
        use CO2System,ONLY: CalcCO2System
        use mem_CO2,ONLY: MethodCalcCO2
        USE BFM_ERROR_MSG, ONLY: BFM_ERROR

        
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
        integer    ::error
        REAL(RLEN)    ::r,s,dumCO2,dumHCO3,dumCO3,dumpCO2,dummy
        logical,save :: start=.TRUE.
        logical,save :: llDIC,llAc,llpH
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! user defined external functions
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        integer, external  :: D3toD1
        integer, external  :: D2toD1
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


        if ( start) then
          start=.FALSE.
          call check_if_in_output('PrDIC',llDIC)
          call check_if_in_output('PrAc',llAc)
          call check_if_in_output('PrpH',llPh)
          llDIC=llDIC.or.llPh
          llAc=llAc.or.llPh
          PrpH(:)=-1.0
        endif


        s=0.0;
        DO BoxNumberZ= 1,NO_BOXES_Z
          DO BoxNumberY=1,NO_BOXES_Y
            DO BoxNumberX=1,NO_BOXES_X
             BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
             BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)
 
             r=abs(seddepth(BoxNumber))
             if (llDIC) PrDIC(BoxNumber)= CalculateFromSet &
                  (KCO2(BoxNumberXY),INTEGRAL,STANDARD,s,r)/(r-s)
             if (llAc)  PrAc(BoxNumber)= CalculateFromSet &
                  (KALK(BoxNumberXY),INTEGRAL,STANDARD,s,r)/(r-s)
             if (llpH) then
               error= CalcCO2System(MethodCalcCO2,ESW_Ben(BoxNumberXY), &
                   ETW_Ben(BoxNumberXY),ERHO_Ben(BoxNumberXY),&
                   PrM1p(BoxNumber),max(ZERO,PrM5s(BoxNumber)),PrAc(BoxNumber),&
                   dumCO2,dumHCO3,dumCO3,PrpH(BoxNumber),&
                   DIC_in=PrDIC(BoxNumber),pCO2_out=dumpCO2)
               if ( error > 0 ) then
!                 write(LOGUNIT,*)"Warning: Ph outside range"
!                 write(LOGUNIT,'(A,'' ='',G12.6)') 'Prn1p',PrM1p(BoxNumber)
!                 write(LOGUNIT,'(A,'' ='',G12.6)') 'PrN5s',PrM5s(BoxNumber)
!                 write(LOGUNIT,'(A,'' ='',G12.6)') 'PrDIC',PrDIC(BoxNumber)
!                 write(LOGUNIT,'(A,'' ='',G12.6)') 'PRAc',PrAc(BoxNumber)
!                 write(LOGUNIT,'(''layer:'',I4,'' pH='',G12.6)') &
!                          BoxNumberXY,PrpH(BoxNumber)
                  Prph(BoxNumber)=-1.0
!                 call BFM_ERROR("BenCO2Profiles",&
!                          "pH outside range 2-11")
               endif
             endif
            s=r;
            ENDDO
          ENDDO
        ENDDO

        return
#endif
#endif
      end

