!$Id: getm_error.F90,v 1.3 2004-04-06 16:54:33 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: error_functions ---
!
! !INTERFACE:
   MODULE bfm_error_msg
!
! !DESCRIPTION:
!   Aim: to have a module and a call to routine wghich has the same
!   structure as the gotm_error_msg
!   This routine in used in standalone and other application than gotm.
!   
! !USES:
   use global_mem
   IMPLICIT NONE

!EOP
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: bfm_error() - global error reporting routine
!
! !INTERFACE:
   subroutine bfm_error(sub,msg)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*),intent(IN)         :: sub,msg
!
!
!-----------------------------------------------------------------------
!BOC
     FATAL "bfm_error: Called from: ",trim(sub)
     FATAL "bfm_error: Message:     ",trim(msg)
     stop  "bfm_error( see bfm.log !! )"

   return
   end subroutine bfm_error
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_warning_for_getm()
!
! !INTERFACE:
   subroutine set_warning_for_getm()
!
! !DESCRIPTION:
!       empty function  if bfm_error_msg is replaced by the gotm_version
!       ( with #define this version can be replaced by the one in gotm_error_message.
!
! !USES:
   IMPLICIT NONE
!
!
!-----------------------------------------------------------------------
!BOC

  return
   end subroutine set_warning_for_getm
!-----------------------------------------------------------------------


end module
!-----------------------------------------------------------------------
! Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
! Copyright (C) 2006 -  Piet Ruardij and Marcello Vichi         !
!-----------------------------------------------------------------------
