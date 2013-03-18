!-----------------------------------------------------------------------
!BOP
!
! !PROGRAM: main
! 
! !INTERFACE:
   PROGRAM main
!
! !DESCRIPTION: 
!
! !USES:
   use standalone
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Author(s): Momme Butenschoen (UNIBO) and Marcello Vichi (INGV)
!
!
! COPYING
!
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
! !LOCAL VARIABLES:
#ifdef IFORT
   real                      :: t1=-1,t2=-1
#endif
! 
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef IFORT
   call CPU_Time(t1)
#endif
      call init_standalone
      call timestepping
      call end_standalone

#ifdef IFORT
   call CPU_Time(t2)
   write(*,*) 'CPU-time was in loop:  ',t2-t1,' seconds'
#endif

      END PROGRAM main
!EOC
