#include "DEBUG.h"


SUBROUTINE PrintSet(NUTR,msg)
     USE mem, ONLY:BoxNumberX_ben,BoxNumberY_ben
     USE mem_BenthicNutrient3, ONLY:bennut_messages
     USE global_mem, ONLY: LOGUNIT
     USE bennut_variables,ONLY:sets
     USE BFM_ERROR_MSG, ONLY: set_warning_for_getm
     implicit none
     integer,intent(IN)             ::NUTR    !Specification`
     character(len=*),intent(IN)    ::msg

     integer                      ::i,k,l,n
    
     write(LOGUNIT,'(''Warnings benthic nutrient model'')')
     write(LOGUNIT,'(''Message;'',A)') msg(1:len_trim(msg))
     call set_warning_for_getm
     if ( bennut_messages == 0 ) return
     write(LOGUNIT,'(''Nutrient Sequence number(NUTR):'',I3)') NUTR
     write(LOGUNIT,'(''Point BoxNumberX_ben,BoxNumberY_ben:'',I5,'','',I5 )') BoxNumberX_ben,BoxNumberY_ben
     write(LOGUNIT,'('' Layer Definition'')')
     write(LOGUNIT, &
       '('' Layers/Equations:'',I5,'' total number of terms:'',I5)') &
          sets(NUTR)%equa, sets(NUTR)%nn

     write( LOGUNIT,'(/,7X,''nr'',6X,''diffusion'',7X,''porosity'',5X,''adsorption   Depth'')') 
     do i= 1,sets(NUTR)%equa
       write( LOGUNIT,'(''layer'',I6,1X,G16.6,1X,G12.6,1X,G12.6,1X,F8.4)') i,  &
              sets(NUTR)%diff(i), sets(NUTR)%poro(i), sets(NUTR)%ads(i),sets(NUTR)%b(i)
     enddo
     write(LOGUNIT,*)

     write(LOGUNIT,'(''Set Definition:'')')
     write(LOGUNIT,'(''nr  term type      coeff'',9x,''lambdas'')') 
     do i= 1,sets(NUTR)%nn
           k=sets(NUTR)%coeffs(i)%il
           l=sets(NUTR)%coeffs(i)%ia
           if ( l>=0)  then
              write(LOGUNIT,'(I2,4x,I2,3x, i2,G16.6)') &
                   i,k,l,sets(NUTR)%factor(i)
           else 
             write(LOGUNIT,'(I2,4x,I2,3x, i2,G16.6,1X,2(G12.6,1X))') &
                 i,k,l,sets(NUTR)%factor(i),  &
                 sets(NUTR)%coeffs(i)%labda(1),sets(NUTR)%coeffs(i)%labda(2)
           endif
     enddo
     write(LOGUNIT,*)
END 
