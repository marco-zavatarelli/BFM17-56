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
  use mem,  ONLY: G13c, G3c, D1m, Q1c, D6m, D2m, D2STATE
  use mem, ONLY: ppG13c, ppG3c, ppD1m, ppQ1c, ppD6m, ppD2m, &
     NO_BOXES_XY,   &
     BoxNumberXY, DICae,  pHAe, pCO2ae, DICan,  pHan, pCO2an,  ETW_Ben, &
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

  do BoxNumberXY=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pH value in oxic sediments
      ! Only the iterative solution of the carbonate system can be
      ! used
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       error= CalcCO2System(DYNAMIC,ESW_Ben(BoxNumberXY),& 
                   ETW_Ben(BoxNumberXY),ERHO_Ben(BoxNumberXY),&
                   M1p(BoxNumberXY),M5s(BoxNumberXY),Acae(BoxNumberXY),&
                   dummy,dummy,dummy,pHae(BoxNumberXY),&
                   DIC_in=DICae(BoxNumberXY),pCO2_out=pCO2ae(BoxNumberXY))
       if ( error > 0 ) then
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ESW_Ben',ESW_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ETW_Ben',ETW_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ERHO_Ben',ERHO_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICae',DICae(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M1p',M1p(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M5s',M5s(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acae',Acae(BoxNumberXY)
            write(LOGUNIT,'('' pHae='',G12.6)') pHae(BoxNumberXY)
            write(LOGUNIT,*) "BenpHDynamics pHae outside range 2-11"
            pHae(BoxNumberXY)=-1
           call BFM_ERROR("BenpHDynamics","pHae outside range 2-11")
       endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pH value in anoxic sediments
      ! Only the iterative solution of the carbonate system can be
      ! used
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       m1=M11p(BoxNumberXY)*(D2m(BoxNumberXY)-D1m(BoxNumberXY)) + &
       M21p(BoxNumberXY)*(p_d_tot-D2m(BoxNumberXY))/ (p_d_tot-D1m(BoxNumberXY))
       error= CalcCO2System(DYNAMIC,ESW_Ben(BoxNumberXY),& 
                   ETW_Ben(BoxNumberXY),ERHO_Ben(BoxNumberXY),&
                   m1,M5s(BoxNumberXY),Acan(BoxNumberXY),&
                   dummy,dummy,dummy,pHan(BoxNumberXY),&
                   DIC_in=DICan(BoxNumberXY),pCO2_out=pCO2an(BoxNumberXY))
       if ( error > 0 ) then
            write(LOGUNIT,*)" Ph outside range"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICan',DICan(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acan',Acan(BoxNumberXY)
            write(LOGUNIT,'('' pHan='',G12.6)') pHan(BoxNumberXY)
            write(LOGUNIT,*) "BenpHDynamics:pHan outside range 2-11"
            pHan(BoxNumberXY)=-1
           call BFM_ERROR("BenpHDynamics","pHan outside range 2-11")
       endif
#ifdef DEBUG
            write(LOGUNIT,*) "in BenpH:"
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ESW_Ben',ESW_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ETW_Ben',ETW_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'ERHO_Ben',ERHO_Ben(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M1p',M1p(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'M5s',M5s(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICae', DICae(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acae',  Acae(BoxNumberXY)
            write(LOGUNIT,'('' pHae='',G12.6)') pHae(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'DICan',DICan(BoxNumberXY)
            write(LOGUNIT,'(A,'' ='',G12.6)') 'Acan',Acan(BoxNumberXY)
            write(LOGUNIT,'('' pHan='',G12.6)') pHan(BoxNumberXY)
#endif
  end do
#endif

  end subroutine BenpHDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
