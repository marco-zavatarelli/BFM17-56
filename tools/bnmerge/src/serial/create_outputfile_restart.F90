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
!    ncresid   : BFM input file identifier
!    ncid      : BFM output merged file identifier 
! --------------------------------------------------------------------------
subroutine create_outputfile_restart

  use netcdf
  use mod_bnmerge
  implicit none
  integer                                :: ncid, ncresid, status
  integer                                :: n
  character(LEN=172)                     :: fname, fname_res
  integer                                :: ndims, nGlobalAtts, nChars
  integer                                :: IDtime,IDtimetmp,IDunlimdim,IDvartime
  integer                                :: IDx, IDy, IDz, IDtarget, IDatt, IDchars
  character(len=64), allocatable         :: res_names(:),res_units(:),res_long(:)

  integer                                :: ndims3dstate, dimvalue3d
  integer                                :: ID3dvars, ID3dstate, ID3dname, ID3dunits, ID3dlong
  integer, allocatable, dimension(:)     :: dimsid3d
  character(len=NF90_MAX_NAME)           :: dimname3d

  real,dimension(1)                   :: time
  character(len = NF90_MAX_NAME)      :: varname,attname
  real, allocatable, dimension(:,:,:,:) :: tmpfillvar3d
  real, allocatable, dimension(:,:,:)   :: tmpfillvar2d  

  integer :: step_start_arr2(2), step_count_arr2(2)


  status = nf90_create(trim(out_dir)//"/"//trim(bfm_restart)//".nc", NF90_NOCLOBBER, ncid)
  if(status /= NF90_NOERR) call handle_err(status,errstring="A file named "//trim(bfm_restart)//".nc already exists!" )

  ! Define the dimensions
  status = nf90_def_dim(ncid, "time", NF90_UNLIMITED, IDtime)
  if(status /= NF90_NOERR) call handle_err(status)
  status = nf90_def_dim(ncid, "x", jpiglo, IDx)
  if (status /= NF90_NOERR) call handle_err(status)
  status = nf90_def_dim(ncid, "y", jpjglo, IDy)
  if (status /= NF90_NOERR) call handle_err(status)
  status = nf90_def_dim(ncid, "depth", jpk, IDz)
  if (status /= NF90_NOERR) call handle_err(status)

  ! read variables from BFM restart domain 0000 and copy attributes
  fname_res = trim(inp_dir)//"/"//trim(bfm_restart)//"_0000.nc" 
  status = nf90_open(path = fname_res, mode = NF90_WRITE, ncid = ncresid)
  if (status /= NF90_NOERR) call handle_err(status,errstring="opening named "//trim(fname_res)//".nc!" )
  status = nf90_inq_dimid(ncresid, "char_max", IDchars)
  if (status /= NF90_NOERR) call handle_err(status)
  status = nf90_inquire_dimension(ncresid, IDChars, len=nChars)
  if (status /= NF90_NOERR) call handle_err(status)
  status = nf90_inquire(ncresid, nDimensions=nDims, nAttributes=nGlobalAtts, unlimitedDimId=IDunlimdim)
  if (status /= NF90_NOERR) call handle_err(status)
  status = nf90_inquire_dimension(ncresid, IDunlimdim)
  if (status /= NF90_NOERR) call handle_err(status)


  ! define geographic variables
  if (ln_mask) then
     call handle_err (nf90_def_var(ncid, "mask", NF90_REAL, (/ IDx, IDy, IDz /), IDtarget))
     call handle_err (nf90_put_att(ncid, IDtarget, "coordinates", "lon lat"))
  endif
  call handle_err (nf90_def_var(ncid, "lat", NF90_REAL, (/ IDx, IDy /), IDtarget))
  call handle_err( nf90_put_att(ncid, IDtarget, "_FillValue", NF90_FILL_REAL) )
  call handle_err (nf90_put_att(ncid, IDtarget, "units", "degrees_north"))
  call handle_err (nf90_def_var(ncid, "lon", NF90_REAL, (/ IDx, IDy /), IDtarget))
  call handle_err( nf90_put_att(ncid, IDtarget, "_FillValue", NF90_FILL_REAL) )
  call handle_err (nf90_put_att(ncid, IDtarget, "units", "degrees_east"))
  call handle_err (nf90_def_var(ncid, "depth", NF90_REAL, (/ IDz /), IDtarget))
  call handle_err( nf90_put_att(ncid, IDtarget, "_FillValue", NF90_FILL_REAL) )
  call handle_err (nf90_put_att(ncid, IDtarget, "long_name", "depth_below_sea"))
  call handle_err (nf90_put_att(ncid, IDtarget, "units", "m"))
  call handle_err (nf90_put_att(ncid, IDtarget, "positive", "down"))
  call handle_err (nf90_put_att(ncid, IDtarget, "axis", "Z"))

  ! copy global attributes
#ifdef DEBUG
  write(*,*) "Creating file:",trim(bfm_restart)//".nc"
  write(*,*) "Start copying global attributes ..."
#endif
  do IDatt=1,nGlobalAtts
     status=nf90_inq_attname(ncresid, NF90_GLOBAL, IDatt, name=attname)
     status = nf90_copy_att(ncresid, NF90_GLOBAL, trim(attname), ncid, NF90_GLOBAL)
     if (status /= NF90_NOERR) call handle_err(status,errstring="copying attribute "//trim(attname))
  end do

  ! extract variable name and attributes
  n_bfmvar3d = 0
  status = nf90_inq_dimid(ncresid, "d3vars", ID3dvars)
  if (status == NF90_NOERR) then
     status = nf90_inquire_dimension(ncresid, ID3dvars, len=n_bfmvar3d)
     if (status /= NF90_NOERR) call handle_err(status)
     write(*,*) "Restart Vars: ",n_bfmvar3d

     ! for pH variable
     n_bfmvar3d = n_bfmvar3d + 1

     allocate(res_names(n_bfmvar3d))
     allocate(res_units(n_bfmvar3d))
     allocate(res_long(n_bfmvar3d))

     status = nf90_inq_varid(ncresid, "D3STATE", ID3dstate)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
     status=nf90_inquire_variable(ncresid, ID3dstate, ndims=ndims3dstate)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
     allocate(dimsid3d(ndims3dstate))

     status=nf90_inquire_variable(ncresid, ID3dstate, dimids=dimsid3d)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
     status = nf90_inquire_dimension(ncresid, dimsid3d(2), name=dimname3d, len=dimvalue3d)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
     deallocate(dimsid3d)

     status = nf90_inq_varid(ncresid, "D3STATE_NAME", ID3dname)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_NAME var in "//fname_res)
     step_start_arr2 = (/ 1, 1 /)
     step_count_arr2 = (/ LEN(res_names), n_bfmvar3d-1 /)
     status = nf90_get_var(ncresid, ID3dname, res_names, start=step_start_arr2, count=step_count_arr2)
     if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_NAME var in "//fname_res)

     status = nf90_inq_varid(ncresid, "D3STATE_UNITS", ID3dunits)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_UNITS var in "//fname_res)
     status = nf90_get_var(ncresid, ID3dunits, res_units, start=step_start_arr2, count=step_count_arr2)
     if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_UNITS var in "//fname_res)

     status = nf90_inq_varid(ncresid, "D3STATE_LONG", ID3dlong)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_LONG var in "//fname_res)
     status = nf90_get_var(ncresid, ID3dlong, res_long, start=step_start_arr2, count=step_count_arr2)
     if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_LONG var in "//fname_res)

     ! Assign dimensions and attributes to variables
     allocate(bfmvarid3d(n_bfmvar3d))

     !Define pH variable in last position in output file
     res_names(n_bfmvar3d) = "pH"
     res_units(n_bfmvar3d) = "pH"
     res_long(n_bfmvar3d) = "acidity or basicity of an aqueous solution"

     !Define the variables inside D3STATE
     do n = 1 , n_bfmvar3d
#ifdef DEBUG
        write(*,*) "Define variable: ",trim(res_names(n))
        write(*,*) "from file ",trim(fname_res)
        write(*,*) "last dimension name ",trim(dimname3d)
#endif
        ! if (dimname3d /= "time" .OR. ndims3dstate == 1) cycle ! enter only with time-varying variables
        ! check the dimension of the variable
        if (dimname3d == "oceanpoint") then
           ! 3D array
           status = nf90_def_var(ncid, trim(res_names(n)), NF90_REAL, (/ IDx, IDy, IDz, IDtime /), IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
           bfmvarid3d(n)= IDtarget
        else
           ! 2D array
           status = nf90_def_var(ncid, trim(res_names(n)), NF90_REAL, (/ IDx, IDy, IDtime /), IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
           bfmvarid3d(n)= IDtarget
        end if
        ! copy attributes
        status = nf90_put_att(ncid, IDtarget, "units", trim(res_units(n)))
        if (status /= NF90_NOERR) call handle_err(status)
        status = nf90_put_att(ncid, IDtarget, "long_name", trim(res_long(n)))
        if (status /= NF90_NOERR) call handle_err(status)
        ! Add fill value
        status = nf90_put_att(ncid, IDtarget, "_FillValue", NF90_FILL_REAL)
        if (status /= NF90_NOERR) call handle_err(status)
        ! Add reference coordinate names (needed with CDO)
        status = nf90_put_att(ncid, IDtarget, "coordinates", "lon lat")
        if (status /= NF90_NOERR) call handle_err(status)
        ! inquire variable
        write(*,*) "Variable saved: ",trim(res_names(n))," : ",n
     end do ! n_bfmvar3d

     deallocate(res_names)
     deallocate(res_units)
     deallocate(res_long)

  end if

  ! copy time variable and attributes
  status = nf90_inq_varid(ncresid, "time", IDtimetmp)
  if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring time var in "//fname)
  status = nf90_get_var(ncresid, IDtimetmp, time, start = (/ 1 /), count = (/ 1 /))
  if (status /= NF90_NOERR) call handle_err(status,errstring="reading time values from"//fname)
  status = nf90_def_var(ncid, "time", NF90_REAL, (/ IDtime /), IDvartime)
  if (status /= NF90_NOERR) call handle_err(status)
  status = nf90_copy_att(ncresid, IDtimetmp, "units", ncid, IDvartime)
  if (status /= NF90_NOERR) call handle_err(status,errstring="copying time units")     

  ! exit definition mode
  status = nf90_enddef(ncid)
  if (status /= NF90_NOERR) call handle_err(status,errstring="while exiting definition mode")

  ! write time values
  status = nf90_put_var(ncid, IDvartime, time(1))
  if (status /= NF90_NOERR) call handle_err(status,errstring="variable: time")

#ifdef DEBUG
  write(*,*) "Writting filling value "
#endif

  !write fill value
  allocate(tmpfillvar3d(jpiglo,jpjglo,jpk,1),tmpfillvar2d(jpiglo,jpjglo,1))
  tmpfillvar3d = NF90_FILL_REAL
  tmpfillvar2d = NF90_FILL_REAL
  do n = 1 , n_bfmvar3d
     IDtarget = bfmvarid3d(n)
     ! inquire ID in the output file
     call handle_err( nf90_inquire_variable(ncid, IDtarget, ndims=ndims, name=varname) )
     ! write empty variable values
     SELECT CASE (ndims)
     case (4)
        ! 3D array
        call handle_err( nf90_put_var(ncid, IDtarget, tmpfillvar3d), &
             errstring="variable:"//trim(varname))
     case (3)
        ! 2D array
        call handle_err( nf90_put_var(ncid, IDtarget, tmpfillvar2d), &
             errstring="variable:"//trim(varname))
     case default
        write(*,'(A,I4,A,A,A,I4)') "invalid dimension size: ",ndims, " for var: '",trim(varname), "'  ID: ", IDtarget
     end select
  end do

  deallocate(tmpfillvar3d,tmpfillvar2d)

#ifdef DEBUG
  write(*,*) "Closing files "
#endif


  ! close the 0000 BFM restart netcdf file
  status = nf90_close(ncresid)
  if (status /= NF90_NOERR) call handle_err(status,errstring="while closing BFM restart input file")
  ! close output file
  status = nf90_close(ncid)
  if (status /= NF90_NOERR) call handle_err(status,errstring="while closing OUTPUT file")

#ifdef DEBUG
  write(*,*) "Output file created "
#endif

end subroutine create_outputfile_restart
