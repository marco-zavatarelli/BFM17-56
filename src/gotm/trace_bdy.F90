!$Id: $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: trace_bdy --- BFM bio model \label{sec:trace_bdy}
!
! !INTERFACE:
   module trace_bdy
!
! !DESCRIPTION:
!  
!
! !USES:
!  default: all is private.
   use bio_var
   private
   integer                        :: ivar,err,rriver
   integer                        :: bdys(4)
   character(len=64),public       :: tracer_type
   character(len=64),public       :: used_for_missing_tracer
! !PUBLIC MEMBER FUNCTIONS:
   public init_trace_bdy, &
          init_var_trace, &
          end_trace_bdy
!

!
! !REVISION HISTORY:
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the template bio module
!
! !INTERFACE:
   subroutine init_trace_bdy(fname,nlev)
!
! !DESCRIPTION:
!  Here, the main communication of array dimensions between GOTM
!  and BFM is done.
!  
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)   :: fname
   integer,          intent(in)   :: nlev
!

   integer             :: n,k,i,rc,ir,jr,lltype,ncid,err
   integer             :: river
   integer             :: unit=410
   real                :: w,wfj,wlj
   logical             :: exist
 

   character(len=64)   :: first_name,second_name

    namelist /trace_nml/ tracer_type,used_for_missing_tracer  

!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_trace_bdy'
    
   inquire(file=fname,exist=exist)

   tracer_type=''
   used_for_missing_tracer=''
   if ( exist) then

     open(unit,file=fname,action='read',status='old',err=99)
     read(unit,nml=trace_nml,err=98)
     close(unit)
   endif

     open(unit,file='riverinfo.dat',action='read',status='old',err=90)
     read(unit,*) i

     ivar=1
     read(unit,*) ir,jr,first_name
     first_name = trim(first_name)
     rriver=1
     do n=2,i
         read(unit,*) ir,jr,second_name
         second_name = trim(second_name)
         if (first_name .ne. second_name) then 
            rriver=rriver+1
            first_name=second_name
         endif
      enddo

     open(unit,file='bdyinfo.dat',status='unknown',ERR=90)
     bdys=0
     do k=1,4 
       read(unit,*,END=91,ERR=92) i
       if (i .ge. 1) then
          bdys(k)=i
          do n = 1,i
            read(unit,*,END=91,ERR=92) w,wfj,wlj
          end do
       end if
     enddo
  91 close(unit)
     ivar=rriver+bdys(1)+bdys(2)+bdys(3)+bdys(4)+1
      numc=ivar
      numcc=numc
      numc_diag=0
      numc_flux=0
      numbc=0
      numbc_diag=0
      numbc_flux=0
      stPelStateS=1
      stPelStateE=numc
      return


  90 stop ' error: when  opening bdyinfo.dat'
  92 stop ' error when reading bdyinfo.dat'
  98 stop ' error when reading trace.nml'
  99  FATAL 'I could not read trace.nml'
      stop 'int_bdy'
 
   end subroutine init_trace_bdy
!EOC

!-----------------------------------------------------------------------
!BOP
!

! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_trace(bio_setup)
!
! !DESCRIPTION:
!  Allocation of BFM variables and initialisation of
!  parameters  and state variables
!
! !USES:
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: bio_setup
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
   integer              :: icontrol,i,ll,end_bdy,j
   character(len=6)    :: bdys_name(4)
   character(len=64)   :: name
   integer             :: unit=410
   integer, allocatable:: r_ids(:)
!!
   include "netcdf.inc"


!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL2 'init_var_trace'
      open(unit,file='riverinfo.dat',action='read',status='old',err=90)
      read(unit,*) i

     ivar=1
     read(unit,*) ir,jr,var_names(ivar)
     var_names(ivar) = trim(var_names(ivar))
     ivar=ivar+1
     do n=2,i
         read(unit,*) ir,jr,var_names(ivar)
         var_names(ivar) = trim(var_names(ivar))
         if (var_names(ivar) .ne. var_names(ivar-1)) ivar=ivar+1
      enddo
      close(unit)

      allocate(r_ids(rriver),stat=err)
      if (err /= 0) stop 'ncdf_river: Error allocating memory (r_ids)'
      err = nf_open('rivers.nc',NCNOWRIT,ncid)
      if (err .ne. NF_NOERR) go to 10
      do n=1,rriver
        err = nf_inq_varid(ncid,var_names(n),r_ids(n))
        if (err .ne. NF_NOERR) go to 10
        err=nf_inq_att(ncid,r_ids(n),"long_name",lltype,ll)
        if (err .ne. NF_NOERR) go to 10
        if ( lltype == NF_CHAR ) then
          err=nf_get_att_text(ncid,r_ids(n),"long_name",name)
          i=index(name,' ')
          if ( i == 0 ) i=ll
          i=min(i,ll)
          write(var_names(n),'(''@riv_'',a)'),name(1:i)
          write(var_long(n),'(''tracer '',A)') name(1:ll)
          write(var_units(n),'(''-'')') 
          if (err .ne. NF_NOERR) go to 10
        endif
      enddo

     err= nf_close(ncid)
     sfl=_ZERO_
     sfl_read=_ZERO_

      cc(1:numc,:)=0.0;
      ivar=rriver+1
      var_names(ivar)="@old"
      write(var_long(ivar),'(''tracer of initial water body'')') 
      write(var_units(ivar),'(''-'')') 
      cc(ivar,:)=1.0
      bdys_name(1)='W'
      bdys_name(2)='N'
      bdys_name(3)='E'
      bdys_name(4)='S'
      do i=1,4
        do n=1,bdys(i)
          ivar=ivar+1
          write(var_names(ivar),'(''@bndy_'',A1,''_'',I2.2)') bdys_name(i),n
          write(var_long(ivar),'(''tracer_'',A1,''_'',I2.2)') bdys_name(i),n
          write(var_units(ivar),'(''-'')') 
        enddo
      enddo

      var_ids=-1
      var_ave=0

    
     return
  90 stop ' error: when  opening bdyinfo.dat'
  92 stop ' error when reading bdyinfo.dat'
  10 stop ' error when reading nc-file'

   end subroutine init_var_trace
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_trace_bdy
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_trace_bdy
!EOC

!-----------------------------------------------------------------------

   end module trace_bdy

!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License 
! www.gnu.org
!-----------------------------------------------------------------------
