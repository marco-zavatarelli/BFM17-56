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
#if defined INCLUDE_BENPROFILES && defined INCLUDE_BENCO2
! !USES:
        use global_mem,      ONLY:RLEN,LOGUNIT,ZERO
#ifdef NOPOINTERS
        use mem
#else
        use mem, ONLY:BoxNumberZ, NO_BOXES_Z, BoxNumberX_ben, NO_BOXES_X, &
           BoxNumberY_ben,NO_BOXES_Y,BoxNumber_ben,BoxNumberXY_ben, seddepth
        use mem,ONLY:PrDIC,PrAc,KALK,KCO2,PrM1p,PrM5s,PrpH, &
           ESW_Ben,ETW_Ben,ERHO_Ben
#endif
        use constants, ONLY: INTEGRAL,STANDARD,MW_C
        use bennut_interface, ONLY:CalculateFromSet
        use CO2System,ONLY: CalcCO2System
        use mem_CO2,ONLY: MethodCalcCO2
        use bfm_error_msg, ONLY: bfm_error
!
! CHANGE_LOG
!   
!
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
!
!

        IMPLICIT  NONE
        integer    ::error
        REAL(RLEN)    ::r,s,dumCO2,dumHCO3,dumCO3,dumpCO2
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
          DO BoxNumberY_ben=1,NO_BOXES_Y
            DO BoxNumberX_ben=1,NO_BOXES_X
             BoxNumber_ben=D3toD1(BoxNumberX_ben,BoxNumberY_ben,BoxNumberZ)
             BoxNumberXY_ben=D2toD1(BoxNumberX_ben,BoxNumberY_ben)
 
             r=abs(seddepth(BoxNumber_ben))
             if (llDIC) PrDIC(BoxNumber_ben)= CalculateFromSet &
                  (KCO2(BoxNumberXY_ben),INTEGRAL,STANDARD,s,r)/(r-s)/ &
                  MW_C/ERHO_Ben(BoxNumberXY_ben)*1000._RLEN
             if (llAc)  PrAc(BoxNumber_ben)= CalculateFromSet &
                  (KALK(BoxNumberXY_ben),INTEGRAL,STANDARD,s,r)/(r-s)/ &
                  ERHO_Ben(BoxNumberXY_ben)*1000._RLEN
             if (llpH) then
               error= CalcCO2System(MethodCalcCO2,ESW_Ben(BoxNumberXY_ben), &
                   ETW_Ben(BoxNumberXY_ben),ERHO_Ben(BoxNumberXY_ben),&
                   PrM1p(BoxNumber_ben),max(ZERO,PrM5s(BoxNumber_ben)),PrAc(BoxNumber_ben),&
                   dumCO2,dumHCO3,dumCO3,PrpH(BoxNumber_ben),&
                   DIC_in=PrDIC(BoxNumber_ben),pCO2_out=dumpCO2)
               if ( error > 0 ) then
                 write(LOGUNIT,*)"Warning: Ph outside range"
                 write(LOGUNIT,'(A,'' ='',G12.6)') 'PrM1p',PrM1p(BoxNumber_ben)
                 write(LOGUNIT,'(A,'' ='',G12.6)') 'PrM5s',PrM5s(BoxNumber_ben)
                 write(LOGUNIT,'(A,'' ='',G12.6)') 'PrDIC',PrDIC(BoxNumber_ben)
                 write(LOGUNIT,'(A,'' ='',G12.6)') 'PrAc',PrAc(BoxNumber_ben)
                 write(LOGUNIT,'(''layer:'',I4,'' pH='',G12.6)') &
                          BoxNumberXY_ben,PrpH(BoxNumber_ben)
                  Prph(BoxNumber_ben)=-1.0
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
      end

