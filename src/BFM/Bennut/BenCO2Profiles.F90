#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.3
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
        USE global_mem,      ONLY:RLEN,LOGUNIT
#ifdef INCLUDE_BENPROFILES
#ifdef INCLUDE_BENCO2
#ifdef NOPOINTERS
        use mem
#else
        USE mem, ONLY:BoxNumberZ, NO_BOXES_Z, BoxNumberX_ben, NO_BOXES_X_BEN, &
           BoxNumberY_ben,NO_BOXES_Y_BEN,BoxNumber_ben,BoxNumberXY_ben, seddepth
        use mem,ONLY:PrDIC,PrAc,KHplus,KCO2,PrM1p,PrM5s,PrpH, &
           ESW_Ben,ETW_Ben,ERHO_Ben
#endif
        USE constants, ONLY: INTEGRAL,STANDARD
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
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi
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
        REAL(RLEN)    ::r,s,dumCO2,dumHCO3,dumCO3,dumpCO2
        logical,static :: start=.TRUE.
        logical,static :: llDIC,llAc,llpH
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
        endif


        s=0.0;
        DO BoxNumberZ= 1,NO_BOXES_Z
          DO BoxNumberY_ben=1,NO_BOXES_Y_BEN
            DO BoxNumberX_ben=1,NO_BOXES_X_BEN
             BoxNumber_ben=D3toD1(BoxNumberX_ben,BoxNumberY_ben,BoxNumberZ)
             BoxNumberXY_ben=D2toD1(BoxNumberX_ben,BoxNumberY_ben)
 
             r=abs(seddepth(BoxNumber_ben))
             if (llDIC) PrDIC(BoxNumber_ben)= CalculateFromSet &
                  (KCO2(BoxNumberXY_ben),INTEGRAL,STANDARD,s,r)/(r-s)/12.0
             if (llAc)  PrAc(BoxNumber_ben)= CalculateFromSet&
                  (KHplus(BoxNumberXY_ben),INTEGRAL,STANDARD,s,r)/(r-s)+2225.0D+00
             if (llpH) then
               error= CalcCO2System(MethodCalcCO2,ESW_Ben(BoxNumberXY_ben), &
                   ETW_Ben(BoxNumberXY_ben),ERHO_Ben(BoxNumberXY_ben),&
                   PrM1p(BoxNumber_ben),PrM5s(BoxNumber_ben),PrAc(BoxNumber_ben),&
                   dumCO2,dumHCO3,dumCO3,PrpH(BoxNumber_ben),&
                   DIC_in=PrDIC(BoxNumber_ben),pCO2_out=dumpCO2)
               if ( error > 0 ) then
                  write(LOGUNIT,*)" Ph outside range"
                  write(LOGUNIT,'(A,'' ='',G12.6)') 'PrN1p',PrM1p(BoxNumber_ben)
                  write(LOGUNIT,'(A,'' ='',G12.6)') 'PrN5s',PrM5s(BoxNumber_ben)
                  write(LOGUNIT,'(A,'' ='',G12.6)') 'PrDIC',PrDIC(BoxNumber_ben)
                  write(LOGUNIT,'(A,'' ='',G12.6)') 'PRAc',PrAc(BoxNumber_ben)
                  write(LOGUNIT,'(''layer:'',I4,'' pH='',G12.6)') &
                           BoxNumberXY_ben,PrpH(BoxNumber_ben)
                  call BFM_ERROR("BenCO2Profiles",&
                           "pH outside range 2-11")
               endif
             endif
            s=r;
            ENDDO
          ENDDO
        ENDDO

        return
#endif
#endif
      end subroutine BenCO2Profiles

