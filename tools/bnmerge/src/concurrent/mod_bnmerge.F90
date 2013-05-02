! BFM_NEMO-MERGE bnmerge
!    Copyright (C) 2009-2011 Marcello Vichi (marcello.vichi@bo.ingv.it)
!    UPDATES:   2012 - Tomas Lovato (tomas.lovato@cmcc.it)
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

module mod_bnmerge

  integer,    parameter, public :: RLEN=selected_real_kind(12,307)
  real(RLEN), parameter, public :: ZERO=0.0_RLEN
  real(RLEN), parameter, public :: ONE=1.0_RLEN

   integer, public ::               &  
      jpi       ,                   &  !: number of grid points along i (proc)
      jpj       ,                   &  !: number of grid points along j (proc)
      jpk       ,                   &  !: number of grid points along k (proc)
      jpiglo    ,                   &  !: number of grid points along i (global domain)
      jpjglo    ,                   &  !: number of grid points along j (global domain)
      jpni      ,                   &  !: number of processors following i 
      jpnj      ,                   &  !: number of processors following j
      jpnij                            !: nb of local domain = nb of processors 
      !                                !  ( <= jpni x jpnj )

   integer, public, allocatable, dimension(:) ::   &  !:
      nimppt, njmppt,  &  !: i-, j-indexes for each processor
      ibonit, ibonjt,  &  !: i-, j- processor neighbour existence
      nlcit , nlcjt,   &  !: dimensions of every subdomain
      nldit , nldjt,   &  !: first, last indoor index for each i-domain
      nleit , nlejt       !: first, last indoor index for each j-domain

   ! masks
   real, public, allocatable, dimension(:,:)                 :: latglo, longlo ! FLOAT
   real(RLEN), public, allocatable, target, dimension(:,:,:) :: maskglo

   ! no. of ocean points in the whole domain (volume and surface)
   integer,public           :: ocepoints, srfpoints

   ! NetCDF IDs  of variables to be merged
   integer,public                                  :: n_bfmvar
   integer,allocatable,dimension(:),public         :: bfmvarid

   ! namelist variables
   character(LEN=400) :: cf_nml_bnmerge='bnmerge.nml'     ! namelist name
   character(LEN=100) :: inp_dir, out_dir, chunk_fname
   integer,parameter  :: NSAVE=120      ! Maximum no variables which can be saved
   character(len=64),dimension(NSAVE):: var_save
   logical :: ln_mask=.FALSE.
   

   public 
   contains 
!
!   ------------------------------------------------------------------------------    
!    Handle errors of NetCDF operations
!   ------------------------------------------------------------------------------
!
   subroutine handle_err(iret,errstring)
   use netcdf
   implicit none
     integer,intent(in)  :: iret
     character(len=*),optional,intent(in) :: errstring
     if (iret .ne. NF90_NOERR) then
       write(*,*) "====== NetCDF Error ======"
       if (present(errstring)) write(*,*) errstring
       write(*,*) NF90_STRERROR(iret)
       stop "stop in function handle_err"
     endif
   end subroutine handle_err
!
!   ------------------------------------------------------------------------------    
!    Retrieve the namelist file name and location
!   ------------------------------------------------------------------------------
!
   subroutine GET_ARGUMENTS()
    
   INTEGER            :: iargc, jarg
   CHARACTER(len=400) :: cr 
   CHARACTER(LEN=2), DIMENSION(1), PARAMETER :: &
        clist_opt = (/ '-f' /)     

   ! --- get arguments
   jarg = 0

   ! set the default name
   if (jarg == iargc()) then
          PRINT *, ' Input namelist not specified !!! '
          PRINT *, ' Set the name to the  default (bnmerge.nml). '
          PRINT *, ' -----'
          PRINT *, ' To specify a different filename use the option -f, e.g. '
          PRINT *, '       bnmerge -f bnmerge2000.nml'
          PRINT *, ''
          cr=trim(cf_nml_bnmerge)
   endif

   DO WHILE ( jarg < iargc() )
      !!
      jarg = jarg + 1
      CALL getarg(jarg,cr)
      !!
      SELECT CASE (trim(cr))
 
      CASE('-f')
        
         IF ( jarg + 1 > iargc() ) THEN
            PRINT *, 'ERROR: Missing namelist name!'; goto 123; 
         ELSE
            !!
            jarg = jarg + 1
            CALL getarg(jarg,cr)
            !!
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: ', trim(cr), ' is definitively not the name of the namelist!'
               goto 123;
            END IF
            !!
            !!
         END IF

      CASE DEFAULT
         PRINT *, ' Unrecognized input option !!! '
         PRINT *, ' Set the name to the  default (bnmerge.nml). ' 
         PRINT *, ' -----'
         PRINT *, ' To specify a the namelis the option -f, e.g. '
         PRINT *, '       bnmerge -f bnmerge2000.nml'
         PRINT *, '' 
         cr=trim(cf_nml_bnmerge)
      END SELECT       
 
   END DO

    !
    ! set the namelist name
    cf_nml_bnmerge=trim(cr)

    PRINT *, ''
    PRINT *, 'BNMERGE Namelist is: ', trim(cf_nml_bnmerge)   
    PRINT *, ''

   return

123 PRINT*, 'IF you see this something is wrong with the specified namelist! ' 
   STOP
   end subroutine GET_ARGUMENTS


end module mod_bnmerge
