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
