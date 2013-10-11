#include "INCLUDE.h"
#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenpH
!
! DESCRIPTION
!   Computation of ph in sediments according to the carbonate system 
!   equations (see ModuleCO2_System.F90)
!
! !INTERFACE
  subroutine BenpHDynamics
!

#ifdef INCLUDE_BENCO2
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,LOGUNIT
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: G13c, G3c, D1m, Q1c, D6m, D2m, D2STATE_BEN
  use mem, ONLY: ppG13c, ppG3c, ppD1m, ppQ1c, ppD6m, ppD2m, &
     NO_BOXES_XY,   &
     BoxNumberXY_ben, DICae,  pHAe, pCO2ae, DICan,  pHan, pCO2an,  ETW_Ben, &
    ESW_Ben, ERHO_Ben, M1p, M5s,AcAe, AcAn,M11p,M21p,D1m,D2m
#endif
  use bfm_error_msg, ONLY: bfm_error
  use CO2System,ONLY: CalcCO2System
  use mem_Param,  ONLY: p_d_tot 
  use mem_CO2, ONLY: DYNAMIC ! force dynamic method for CO2 equilibria computation
  IMPLICIT NONE
!  
! !LOCAL VARIABLES
  real(RLEN)  :: r
  real(RLEN)  :: m1
  real(RLEN)  :: dummy
  integer     :: error
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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
!EOP
!-------------------------------------------------------------------------!
!BOC

  do BoxNumberXY_ben=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pH value in oxic sediments
      ! Only the iterative solution of the carbonate system can be
      ! used
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       error= CalcCO2System(DYNAMIC, &
            ESW_Ben(BoxNumberXY_ben), & 
            ETW_Ben(BoxNumberXY_ben), &
            ERHO_Ben(BoxNumberXY_ben), &
            M1p(BoxNumberXY_ben), &
            M5s(BoxNumberXY_ben), &
            Acae(BoxNumberXY_ben),&
            dummy,dummy,dummy,pHae(BoxNumberXY_ben),&
            DIC_in=DICae(BoxNumberXY_ben), &
            pCO2_out=pCO2ae(BoxNumberXY_ben))
       if ( error > 0 ) then
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ESW_Ben',ESW_Ben(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ETW_Ben',ETW_Ben(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ERHO_Ben',ERHO_Ben(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICae',DICae(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M1p',M1p(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M5s',M5s(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acae',Acae(BoxNumberXY_ben)
            write(LOGUNIT,'('' pHae='',G12.6)') pHae(BoxNumberXY_ben)
            write(LOGUNIT,*) "BenpHDynamics pHae outside range 2-11"
            pHae(BoxNumberXY_ben)=-1
           call BFM_ERROR("BenpHDynamics","pHae outside range 2-11")
       endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pH value in anoxic sediments
      ! Only the iterative solution of the carbonate system can be
      ! used
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       m1=M11p(BoxNumberXY_ben)* &
            (D2m(BoxNumberXY_ben)- &
            D1m(BoxNumberXY_ben))+ &
            M21p(BoxNumberXY_ben)* &
            (p_d_tot-D2m(BoxNumberXY_ben))/ &
            (p_d_tot-D1m(BoxNumberXY_ben))
       error= CalcCO2System(DYNAMIC, &
            ESW_Ben(BoxNumberXY_ben), & 
            ETW_Ben(BoxNumberXY_ben), &
            ERHO_Ben(BoxNumberXY_ben), &
            m1,M5s(BoxNumberXY_ben), &
            Acan(BoxNumberXY_ben), &
            dummy,dummy,dummy,pHan(BoxNumberXY_ben), &
            DIC_in=DICan(BoxNumberXY_ben), &
            pCO2_out=pCO2an(BoxNumberXY_ben))
       if ( error > 0 ) then
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICan',DICan(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acan',Acan(BoxNumberXY_ben)
            write(LOGUNIT,'('' pHan='',G12.6)') pHan(BoxNumberXY_ben)
            write(LOGUNIT,*) "BenpHDynamics:pHan outside range 2-11"
            pHan(BoxNumberXY_ben)=-1
           call BFM_ERROR("BenpHDynamics","pHan outside range 2-11")
       endif
#ifdef DEBUG
            write(LOGUNIT,*) "in BenpH:"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ESW_Ben',ESW_Ben(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ETW_Ben',ETW_Ben(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ERHO_Ben',ERHO_Ben(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M1p',M1p(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M5s',M5s(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICae', DICae(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acae',  Acae(BoxNumberXY_ben)
            write(LOGUNIT,'('' pHae='',G12.6)') pHae(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICan',DICan(BoxNumberXY_ben)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acan',Acan(BoxNumberXY_ben)
            write(LOGUNIT,'('' pHan='',G12.6)') pHan(BoxNumberXY_ben)
#endif
  end do
#endif

  end subroutine BenpHDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
