! BFM_NEMO-MERGE bnmerge
!    Copyright (C) 2009-2011 Marcello Vichi (marcello.vichi@bo.ingv.it)
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
! --------------------------------------------------------------------------
!    Notes 
!    ncbfmid   : BFM input file identifier
!    ncid      : BFM output merged file identifier 
! --------------------------------------------------------------------------
subroutine create_outputfile

  use pnetcdf
  use mod_bnmerge
  use mpi

  implicit none
  integer           :: ncid, ncbfmid, status
  integer           :: p,n,d,t
  character(LEN=172) :: fname
  integer           :: ndims, nVars, nGlobalAtts
  integer           :: IDtime,IDtimetmp,IDunlimdim,IDvartime
  integer           :: IDx, IDy, IDz, IDvar, IDtarget, IDatt
  integer(kind=MPI_OFFSET_KIND) :: ntime, length
  real(4), allocatable, dimension(:) :: time
  integer           :: npoints
  integer :: vartype,dimids(4),dimlen(4), natts
  character(len = NF_MAX_NAME) :: DimName,varname,attname
  integer(kind=MPI_OFFSET_KIND) :: time_start(1), time_count(1)
  real(4) :: fillValue(1)

  real(4), allocatable, dimension(:,:,:,:) :: tmpfillvar3d
  real(4), allocatable, dimension(:,:,:)   :: tmpfillvar2d


  ! Define missing value
  fillvalue=(/ NF_FILL_REAL /)

  call handle_err( nfmpi_create(MPI_COMM_WORLD, trim(out_dir)//"/"//trim(chunk_fname)//".nc", NF_NOCLOBBER, MPI_INFO_NULL, ncid), &
       errstring="A file named "//trim(chunk_fname)//".nc already exists!" )

  ! Define the dimensions
  call handle_err( nfmpi_def_dim(ncid, 'time', NF_UNLIMITED, IDtime) )
  call handle_err( nfmpi_def_dim(ncid, "x", jpiglo, IDx) )
  call handle_err( nfmpi_def_dim(ncid, "y", jpjglo, IDy) )
  call handle_err( nfmpi_def_dim(ncid, "depth", jpk, IDz) )

  ! read variables from domain 0000 and copy attributes
  fname = trim(inp_dir)//"/"//trim(chunk_fname)//"_0000.nc" 
  call handle_err( nfmpi_open(MPI_COMM_WORLD, fname, NF_NOWRITE, MPI_INFO_NULL, ncbfmid) )
  call handle_err( nfmpi_inq(ncbfmid, nDims, nVars, nGlobalAtts, IDunlimdim) )
  call handle_err( nfmpi_inq_dimlen(ncbfmid, IDunlimdim, ntime) )

  ! define geographic variables
  if (ln_mask) then
     call handle_err ( nfmpi_def_var(ncid, "mask", NF_REAL, 3, (/ IDx, IDy, IDz /), IDtarget) )
     length = len("lon lat")
     call handle_err ( nfmpi_put_att_text(ncid, IDtarget, "coordinates", length, "lon lat") )
  endif

  call handle_err ( nfmpi_def_var(ncid, "lat", NF_REAL, 2, (/ IDx, IDy /), IDtarget) )
  call handle_err( nfmpi_put_att_real(ncid, IDtarget, "_FillValue", NF_REAL, 1, fillvalue) )
  length = len("degrees_north")
  call handle_err ( nfmpi_put_att_text(ncid, IDtarget, "units", length, "degrees_north") )
  call handle_err ( nfmpi_def_var(ncid, "lon", NF_REAL, 2, (/ IDx, IDy /), IDtarget) )
  call handle_err( nfmpi_put_att_real(ncid, IDtarget, "_FillValue", NF_REAL, 1, fillvalue) )
  length = len("degrees_east")
  call handle_err ( nfmpi_put_att_text(ncid, IDtarget, "units", length, "degrees_east") )
  call handle_err ( nfmpi_def_var(ncid, "depth", NF_REAL, 1, (/ IDz /), IDtarget) )
  call handle_err( nfmpi_put_att_real(ncid, IDtarget, "_FillValue", NF_REAL, 1, fillvalue) )
  length = len("depth_below_sea")
  call handle_err ( nfmpi_put_att_text(ncid, IDtarget, "long_name", length, "depth_below_sea") )
  length = len("m")
  call handle_err ( nfmpi_put_att_text(ncid, IDtarget, "units", length, "m") )
  length = len("down")
  call handle_err ( nfmpi_put_att_text(ncid, IDtarget, "positive", length, "down") )
  length = len("Z")
  call handle_err ( nfmpi_put_att_text(ncid, IDtarget, "axis", length, "Z") )

  ! copy global attributes
#ifdef DEBUG
  write(*,*) "Creating file:",trim(chunk_fname)//".nc"," containing",ntime,"time frames"
  write(*,*) "Start copying global attributes ..."
#endif
  do IDatt=1,nGlobalAtts
     call handle_err( nfmpi_inq_attname(ncbfmid, NF_GLOBAL, IDatt, name=attname), &
          errstring="inquiring attribute "//trim(attname))
     call handle_err( nfmpi_copy_att(ncbfmid, NF_GLOBAL, trim(attname), ncid, NF_GLOBAL), &
          errstring="copying attribute "//trim(attname))
  end do
  ! Tracks of the variables that have to be stored
  allocate(bfmvarid(nVars))
  n_bfmvar=0
  do IDvar=1,nVars
     call handle_err( nfmpi_inq_var(ncbfmid, IDvar, datatype=vartype, ndims=ndims, dimids=dimids, name=varname, natts=natts), &
          errstring="variable: "//trim(varname) )
     do n = 1 , NSAVE
        if ( trim(var_save(n)) == trim(varname) ) then
           n_bfmvar = n_bfmvar + 1
           bfmvarid(n_bfmvar)= IDvar
#ifdef DEBUG
           write(*,*) "Assigned output ",trim(varname), " with ID:", IDvar
#endif
        endif
     enddo
  enddo
#ifdef DEBUG
  write(*,*) "Total Output Variables ", n_bfmvar
  write(*,*)
#endif
  if (n_bfmvar == 0) then
     write(*,*) "Selected output variables do not match the content of the input files.", n_bfmvar
     stop
  endif


  ! Assign dimensions and attributes to variables 
  do n = 1 , n_bfmvar
     IDvar = bfmvarid(n)

     ! call handle_err( nfmpi_redef(ncid), &
     !      errstring="enter def mode in variable" )

     ! inquire variable
     call handle_err( nfmpi_inq_var(ncbfmid, IDvar, datatype=vartype, ndims=ndims, dimids=dimids(1:ndims), name=varname, natts=natts), &
          errstring="variable: "//trim(varname) )
     call handle_err( nfmpi_inq_dimname(ncbfmid, dimids(ndims), name=DimName) )

#ifdef DEBUG
     write(*,*) "Define variable: ",trim(varname)," with ID: ",IDvar
     write(*,*) "from file ",trim(fname)
     write(*,*) "last dimension name ",trim(DimName), " nun dims: ", ndims
#endif

     if (DimName /= "time" .OR. ndims == 1) cycle ! enter only with time-varying variables

     ! check the dimension of the variable
     call handle_err( nfmpi_inq_dimname(ncbfmid, dimids(1), name=DimName) )

     if (DimName == "oceanpoint") then
        ! 3D array
        call handle_err( nfmpi_def_var(ncid, trim(varname), NF_REAL, 4, (/ IDx, IDy, IDz, IDtime /), IDtarget),errstring="3D var" )
     else
        ! 2D array
        call handle_err( nfmpi_def_var(ncid, trim(varname), NF_REAL, 3, (/ IDx, IDy, IDtime /), IDtarget),errstring="2D var" )
     end if

     ! copy attributes
     call handle_err( nfmpi_copy_att(ncbfmid, IDvar, "units", ncid, IDtarget) )
     call handle_err( nfmpi_copy_att(ncbfmid, IDvar, "long_name", ncid, IDtarget) )

     ! Add fill value
     call handle_err( nfmpi_put_att_real(ncid, IDtarget, "_FillValue", NF_REAL, 1, fillvalue) )

     ! Add reference coordinate names (needed with CDO)
     length = len("lon lat")
     call handle_err( nfmpi_put_att_text(ncid, IDtarget, "coordinates", length, "lon lat") )

  end do ! n_bfmvar 

  ! copy time variable and attributes
  allocate(time(ntime))
  time_start(1) = 1
  time_count(1) = ntime
  call handle_err( nfmpi_inq_varid(ncbfmid, "time", IDtimetmp), &
       errstring="inquiring time var in "//fname)
  call handle_err( nfmpi_begin_indep_data(ncbfmid) )
  call handle_err( nfmpi_get_vara_real(ncbfmid, IDtimetmp, time_start, time_count, time), &
       errstring="reading time values from"//fname )
  call handle_err( nfmpi_end_indep_data(ncbfmid) )
  call handle_err( nfmpi_def_var(ncid, "time", NF_REAL, 1, (/ IDtime /), IDvartime) )
  call handle_err( nfmpi_copy_att(ncbfmid, IDtimetmp, "units", ncid, IDvartime), &
       errstring="copying time units" )
  ! exit definition mode
  call handle_err( nfmpi_enddef(ncid), &
       errstring="while exiting definition mode" )
  ! write time values
  call handle_err( nfmpi_begin_indep_data(ncid) )
  call handle_err( nfmpi_put_vara_real(ncid, IDvartime, time_start, time_count, time), &
       errstring="variable: time")
  call handle_err( nfmpi_end_indep_data(ncid) )


  !write fill value
  allocate(tmpfillvar3d(jpiglo,jpjglo,jpk,ntime),tmpfillvar2d(jpiglo,jpjglo,ntime))
  tmpfillvar3d = NF_FILL_REAL
  tmpfillvar2d = NF_FILL_REAL
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP PRIVATE ( IDvar, varname, IDtarget, ndims) &
  !$OMP SHARED   ( n_bfmvar, bfmvarid, tmpfillvar3d, tmpfillvar2d, ncbfmid, ncid )
  do n = 1 , n_bfmvar
     IDvar = bfmvarid(n)
     ! inquire ID in the output file
     call handle_err( nfmpi_inq_varname(ncbfmid, IDvar, varname) )
     call handle_err( nfmpi_inq_varid(ncid, varname, IDtarget) )
     ! check the dimension of the variable
     call handle_err( nfmpi_inq_varndims(ncid, IDtarget, ndims) )
     ! write empty variable values
     SELECT CASE (ndims)
     case (4)
        ! 3D array
        !$OMP CRITICAL
        call handle_err( nfmpi_begin_indep_data(ncid) )
        call handle_err( nfmpi_put_var_real(ncid, IDtarget, tmpfillvar3d), &
             errstring="variable:"//trim(varname))
        call handle_err( nfmpi_end_indep_data(ncid) )
        !$OMP END CRITICAL
     case (3)
        ! 2D array
        !$OMP CRITICAL
        call handle_err( nfmpi_begin_indep_data(ncid) )
        call handle_err( nfmpi_put_var_real(ncid, IDtarget, tmpfillvar2d), &
             errstring="variable:"//trim(varname))
        call handle_err( nfmpi_end_indep_data(ncid) )
        !$OMP END CRITICAL
     case default
        write(*,'(A,I4)') "invalid dimension size: ",ndims
     end select
  end do
  !$OMP END PARALLEL DO


  deallocate(tmpfillvar3d,tmpfillvar2d)

  ! close the 0000 BFM netcdf file
  call handle_err( nfmpi_close(ncbfmid), errstring="while closing BFM input file" )
  call handle_err( nfmpi_close(ncid),    errstring="while closing BFM output file" )



#ifdef DEBUG
  write(*,*) "Output file created with ",ntime,"time frames"
#endif

  return

end subroutine create_outputfile
