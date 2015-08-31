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

  use netcdf
  implicit none

#ifdef PARALLEL
  include 'mpif.h'
! #define INTEGER integer(kind=MPI_OFFSET_KIND)
#endif


  integer,    parameter, public :: RLEN=selected_real_kind(12,307)
  real(RLEN), parameter, public :: ZERO=0.0_RLEN
  real(RLEN), parameter, public :: ONE=1.0_RLEN

  integer, public :: &  
       jpkglo    ,                   &  !: number of grid points along k (proc)
       jpiglo    ,                   &  !: number of grid points along i (global domain)
       jpjglo    ,                   &  !: number of grid points along j (global domain)
       jpnij                            !: nb of local domain = nb of processors 
  !                                !  ( <= jpni x jpnj )

  integer, public, allocatable, dimension(:) ::   &  !:
       nimppt, njmppt,  &  !: i-, j-indexes for each processor
       nlcit , nlcjt       !: dimensions of every subdomain

  ! masks
  real, public, allocatable, target, dimension(:,:)   :: latglo, longlo ! FLOAT
  real, public, allocatable, target, dimension(:,:,:) :: maskglo

  ! NetCDF IDs, types and dimensions  of variables to be merged
  integer,public                          :: n_bfmvar_out, n_bfmvar_res
  integer,allocatable,dimension(:),public :: bfmvarid_out, bfmvarid_res         ! id of variables in input
  integer,allocatable,dimension(:),public :: bfmvartype_out, bfmvartype_res     ! type of var dimension
  integer,allocatable,dimension(:),public :: bfmvartarget_out, bfmvartarget_res ! id of variables in output
  integer, parameter, public              :: TYPE_OCE=1, TYPE_BTN=2, TYPE_SRF=3, TYPE_PH=4 ! type of input dimensions

  ! namelist variables
  character(LEN=400) :: cf_nml_bnmerge='bnmerge.nml'     ! namelist name
  character(LEN=NF90_MAX_NAME)   :: inp_dir, out_dir, chunk_fname='', bfm_restart='', out_fname=''
  logical :: do_restart, do_output
  integer,parameter  :: NSAVE=120      ! Maximum no variables which can be saved
  character(len=64),dimension(NSAVE):: var_save="NotVar"
  logical :: ln_mask=.FALSE.
  logical :: nc_compres=.FALSE.
  integer,public :: nc_shuffle=0,nc_deflate=0,nc_defllev=0

  public 
contains 
  !
  !   ------------------------------------------------------------------------------    
  !    Handle errors of NetCDF operations
  !   ------------------------------------------------------------------------------
!
  subroutine handle_err(iret,errstring)
    implicit none
    integer,intent(in)  :: iret
    character(len=*),optional,intent(in) :: errstring
#ifdef PARALLEL
    if (iret .ne. NF90_NOERR) then
       nf90_set_log_level(6)
       write(*,*) "====== NetCDF Error ======"
       if (present(errstring)) write(*,*) errstring
       write(*,*) NF90_STRERROR(iret)
       call MPI_ABORT(MPI_COMM_WORLD, 1, iret)
#else
    if (iret .ne. NF90_NOERR) then
       write(*,*) "====== NetCDF Error ======"
       if (present(errstring)) write(*,*) errstring
       write(*,*) NF90_STRERROR(iret)
       stop "stop in function handle_err"
#endif
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
  
  subroutine tick(t)
    implicit none
    integer, intent(OUT) :: t
    call system_clock(t)
  end subroutine tick

  ! returns time in seconds from now to time described by t
  real(8) function tock(t)
    implicit none
    integer, intent(in) :: t
    integer :: now, clock_rate

    call system_clock(now,clock_rate)

    tock = real(now - t)/real(clock_rate)
  end function tock


  subroutine replace_char(str,tar,rep)
    implicit none
    character(LEN=*), intent(INOUT) :: str
    character(LEN=*), intent(IN)    :: tar, rep

    integer                 :: times

    times = scan(str, tar)
    do while ( times .ne. 0 )
       str(times:times) = rep
       times = scan(str, tar)
    end do
  end subroutine replace_char


end module mod_bnmerge
