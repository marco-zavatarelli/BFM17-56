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
!    ncoutid   : BFM input file identifier
!    ncresid   : BFM input restart file identifier
!    ncidout   : BFM output merged file identifier 
!    ncidres   : BFM output restart merged file identifier 
! --------------------------------------------------------------------------
subroutine create_outputfile

  use netcdf
  use mod_bnmerge, ONLY: n_bfmvar_out, n_bfmvar_res, bfmvarid_out, bfmvarid_res, &
       inp_dir, out_dir, chunk_fname, bfm_restart, &
       do_restart, do_output, &
       jpiglo, jpjglo, jpkglo, &
       ln_mask, NSAVE, &
       handle_err, var_save

  implicit none
  integer            :: ncidout, ncidres, ncoutid, status, ncresid
  integer            :: n
  character(LEN=172) :: fname_out, fname_res
  integer            :: ndims, nVars_out, nGlobalAtts_out, nGlobalAtts_res, nChars

  integer           :: IDtime,IDtimetmp,IDunlimdim_out,IDunlimdim_res,IDvartime_out,IDvartime_res
  integer           :: IDx, IDy, IDz, IDvar, IDtarget, IDatt, IDchars

  character(len=64), allocatable         :: res_names(:),res_units(:),res_long(:)
  integer                                :: ndims_state_res
  integer                                :: IDvars_res, IDstate_res, IDname_res, IDunits_res, IDlong_res
  integer, allocatable, dimension(:)     :: dimsid_res

  integer                         :: ntime_res, ntime_out
  real, allocatable, dimension(:) :: time_res,time_out

  integer :: dimids(4)
  character(len = NF90_MAX_NAME) :: DimName,dimname3d_res, varname,attname

  integer :: step_start_arr2(2), step_count_arr2(2)
  real, allocatable, dimension(:,:,:,:) :: tmpfillvar3d
  real, allocatable, dimension(:,:,:)   :: tmpfillvar2d  


  if( do_output ) then
     status = nf90_create(trim(out_dir)//"/"//trim(chunk_fname)//".nc", NF90_NOCLOBBER, ncidout)
     if(status /= NF90_NOERR) call handle_err(status,errstring="A file named "//trim(chunk_fname)//".nc already exists!" )
  end if
  if( do_restart ) then
     status = nf90_create(trim(out_dir)//"/"//trim(bfm_restart)//".nc", NF90_NOCLOBBER, ncidres)
     if(status /= NF90_NOERR) call handle_err(status,errstring="A file named "//trim(bfm_restart)//".nc already exists!" )
  end if


  ! Define the dimensions
  if( do_output ) then
     write(*,*) "Preparing output"
     status = nf90_def_dim(ncidout, "time", NF90_UNLIMITED, IDtime)
     if(status /= NF90_NOERR) call handle_err(status,errstring="reading time in out")
     status = nf90_def_dim(ncidout, "x", jpiglo, IDx)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading X in out")
     status = nf90_def_dim(ncidout, "y", jpjglo, IDy)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading Y in out")
     status = nf90_def_dim(ncidout, "depth", jpkglo, IDz)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading depth in out")
  endif
  if( do_restart ) then
     write(*,*) "Preparing restart..."
     status = nf90_def_dim(ncidres, "time", NF90_UNLIMITED, IDtime)
     if(status /= NF90_NOERR) call handle_err(status,errstring="reading time in restart out")
     status = nf90_def_dim(ncidres, "x", jpiglo, IDx)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading X in restart out")
     status = nf90_def_dim(ncidres, "y", jpjglo, IDy)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading Y in restart out")
     status = nf90_def_dim(ncidres, "depth", jpkglo, IDz)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading depth in restart out")
  endif


  if( do_output ) then
     ! read variables from domain 0000 and copy attributes
     fname_out = trim(inp_dir)//"/"//trim(chunk_fname)//"_0000.nc" 
     status = nf90_open(path = fname_out, mode = NF90_WRITE, ncid = ncoutid)
     if (status /= NF90_NOERR) call handle_err(status,errstring="opening named "//trim(chunk_fname)//".nc!" )
     status = nf90_inquire(ncoutid, nVariables=nVars_out, nAttributes=nGlobalAtts_out, unlimitedDimId=IDunlimdim_out)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading info")
     status = nf90_inquire_dimension(ncoutid, IDunlimdim_out, len = ntime_out)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading time")
  end if
  if( do_restart) then
     ! read variables from BFM restart domain 0000 and copy attributes
     fname_res = trim(inp_dir)//"/"//trim(bfm_restart)//"_0000.nc" 
     status = nf90_open(path = fname_res, mode = NF90_WRITE, ncid = ncresid)
     if (status /= NF90_NOERR) call handle_err(status,errstring="opening named "//trim(fname_res)//".nc!" )
     status = nf90_inq_dimid(ncresid, "char_max", IDchars)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading chars")
     status = nf90_inquire_dimension(ncresid, IDChars, len=nChars)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading len chars")
     status = nf90_inquire(ncresid, nAttributes=nGlobalAtts_res, unlimitedDimId=IDunlimdim_res)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading info")
     status = nf90_inquire_dimension(ncresid, IDunlimdim_res, len = ntime_res)
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading time")
  end if

  ! define geographic variables
  if( do_output ) then
     if ( ln_mask ) then
        call handle_err (nf90_def_var(ncidout, "mask", NF90_REAL, (/ IDx, IDy, IDz /), IDtarget))
        call handle_err (nf90_put_att(ncidout, IDtarget, "coordinates", "lon lat"))
     end if
     call handle_err (nf90_def_var(ncidout, "lat", NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err( nf90_put_att(ncidout, IDtarget, "_FillValue", NF90_FILL_REAL) )
     call handle_err (nf90_put_att(ncidout, IDtarget, "units", "degrees_north"))
     call handle_err (nf90_def_var(ncidout, "lon", NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err( nf90_put_att(ncidout, IDtarget, "_FillValue", NF90_FILL_REAL) )
     call handle_err (nf90_put_att(ncidout, IDtarget, "units", "degrees_east"))
     call handle_err (nf90_def_var(ncidout, "depth", NF90_REAL, (/ IDz /), IDtarget))
     call handle_err( nf90_put_att(ncidout, IDtarget, "_FillValue", NF90_FILL_REAL) )
     call handle_err (nf90_put_att(ncidout, IDtarget, "long_name", "depth_below_sea"))
     call handle_err (nf90_put_att(ncidout, IDtarget, "units", "m"))
     call handle_err (nf90_put_att(ncidout, IDtarget, "positive", "down"))
     call handle_err (nf90_put_att(ncidout, IDtarget, "axis", "Z"))
  end if
  if ( do_restart ) then
     if ( ln_mask ) then
        call handle_err (nf90_def_var(ncidres, "mask", NF90_REAL, (/ IDx, IDy, IDz /), IDtarget))
        call handle_err (nf90_put_att(ncidres, IDtarget, "coordinates", "lon lat"))
     end if
     call handle_err (nf90_def_var(ncidres, "lat", NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err( nf90_put_att(ncidres, IDtarget, "_FillValue", NF90_FILL_REAL) )
     call handle_err (nf90_put_att(ncidres, IDtarget, "units", "degrees_north"))
     call handle_err (nf90_def_var(ncidres, "lon", NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err( nf90_put_att(ncidres, IDtarget, "_FillValue", NF90_FILL_REAL) )
     call handle_err (nf90_put_att(ncidres, IDtarget, "units", "degrees_east"))
     call handle_err (nf90_def_var(ncidres, "depth", NF90_REAL, (/ IDz /), IDtarget))
     call handle_err( nf90_put_att(ncidres, IDtarget, "_FillValue", NF90_FILL_REAL) )
     call handle_err (nf90_put_att(ncidres, IDtarget, "long_name", "depth_below_sea"))
     call handle_err (nf90_put_att(ncidres, IDtarget, "units", "m"))
     call handle_err (nf90_put_att(ncidres, IDtarget, "positive", "down"))
     call handle_err (nf90_put_att(ncidres, IDtarget, "axis", "Z"))
  end if

  ! copy global attributes
  if ( do_output ) then
#ifdef DEBUG
     write(*,*) "Creating file:",trim(chunk_fname)//".nc"," containing",ntime_out,"time frames"
     write(*,*) "Start copying global attributes ..."
#endif
     do IDatt=1,nGlobalAtts_out
        status=nf90_inq_attname(ncoutid, NF90_GLOBAL, IDatt, name=attname)
        status = nf90_copy_att(ncoutid, NF90_GLOBAL, trim(attname), ncidout, NF90_GLOBAL)
        if (status /= NF90_NOERR) call handle_err(status,errstring="copying attribute "//trim(attname))
     end do
  end if
  if ( do_restart ) then
#ifdef DEBUG
     write(*,*) "Creating file:",trim(bfm_restart)//".nc"," containing",ntime_res,"time frames"
     write(*,*) "Start copying global attributes ..."
#endif
     do IDatt=1,nGlobalAtts_res
        status=nf90_inq_attname(ncresid, NF90_GLOBAL, IDatt, name=attname)
        status = nf90_copy_att(ncresid, NF90_GLOBAL, trim(attname), ncidres, NF90_GLOBAL)
        if (status /= NF90_NOERR) call handle_err(status,errstring="copying attribute "//trim(attname))
     end do
  end if


  ! Tracks of the variables that have to be stored
  if ( do_output ) then 
     allocate(bfmvarid_out(nVars_out))
     n_bfmvar_out=0
     do IDvar=1,nVars_out
        status=nf90_inquire_variable(ncoutid, IDvar, ndims=ndims, name=varname)
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        do n = 1 , NSAVE
           if ( trim(var_save(n)) == trim(varname) ) then
              n_bfmvar_out = n_bfmvar_out + 1
              bfmvarid_out(n_bfmvar_out)= IDvar
              write(*,*) "Assigned output ",trim(varname), " with ID:", IDvar
           endif
        enddo
     enddo
     write(*,*) "Total Output Variables ", n_bfmvar_out
     write(*,*)
     if (n_bfmvar_out == 0) then
        write(*,*) "Selected output variables do not match the content of the input files.", n_bfmvar_out
        stop
     endif
     ! Assign dimensions and attributes to variables 
     do n = 1 , n_bfmvar_out
        IDvar = bfmvarid_out(n)
        ! inquire variable
        status=nf90_inquire_variable(ncoutid, IDvar, ndims=ndims, name=varname)
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        status=nf90_inquire_variable(ncoutid, IDvar, dimids=dimids(1:ndims))
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        status = nf90_inquire_dimension(ncoutid, dimids(ndims), name = DimName)
        write(*,*) "Define variable: ",trim(varname)," with ID: ",IDvar
#ifdef DEBUG
        write(*,*) "from file ",trim(fname_out)
        write(*,*) "last dimension name ",trim(DimName)
#endif
        if (DimName /= "time" .OR. ndims == 1) cycle ! enter only with time-varying variables
        ! check the dimension of the variable
        status = nf90_inquire_dimension(ncoutid, dimids(1), name = DimName)
        if (DimName == "oceanpoint") then
           ! 3D array
           status = nf90_def_var(ncidout, trim(varname), NF90_REAL, (/ IDx, IDy, IDz, IDtime /), IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
        else
           ! 2D array
           status = nf90_def_var(ncidout, trim(varname), NF90_REAL, (/ IDx, IDy, IDtime /), IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
        end if
        ! copy attributes
        status = nf90_copy_att(ncoutid, IDvar, "units", ncidout, IDtarget)
        if (status /= NF90_NOERR) call handle_err(status)
        status = nf90_copy_att(ncoutid, IDvar, "long_name", ncidout, IDtarget)
        if (status /= NF90_NOERR) call handle_err(status)
        ! Add fill value
        status = nf90_put_att(ncidout, IDtarget, "_FillValue", NF90_FILL_REAL)
        if (status /= NF90_NOERR) call handle_err(status)
        ! Add reference coordinate names (needed with CDO)
        status = nf90_put_att(ncidout, IDtarget, "coordinates", "lon lat")
        if (status /= NF90_NOERR) call handle_err(status)
     end do ! n_bfmvar_out 
  end if
  if ( do_restart ) then
     n_bfmvar_res = 0
     status = nf90_inq_dimid(ncresid, "d3vars", IDvars_res)
     if (status == NF90_NOERR) then
        status = nf90_inquire_dimension(ncresid, IDvars_res, len=n_bfmvar_res)
        if (status /= NF90_NOERR) call handle_err(status)
        write(*,*) "Restart Vars: ",n_bfmvar_res

        ! for pH variable
        n_bfmvar_res = n_bfmvar_res + 1

        allocate(res_names(n_bfmvar_res))
        allocate(res_units(n_bfmvar_res))
        allocate(res_long(n_bfmvar_res))

        status = nf90_inq_varid(ncresid, "D3STATE", IDstate_res)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
        status=nf90_inquire_variable(ncresid, IDstate_res, ndims=ndims_state_res)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
        allocate(dimsid_res(ndims_state_res))

        status=nf90_inquire_variable(ncresid, IDstate_res, dimids=dimsid_res)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
        status = nf90_inquire_dimension(ncresid, dimsid_res(2), name=dimname3d_res)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var in "//fname_res)
        deallocate(dimsid_res)

        status = nf90_inq_varid(ncresid, "D3STATE_NAME", IDname_res)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_NAME var in "//fname_res)
        step_start_arr2 = (/ 1, 1 /)
        step_count_arr2 = (/ LEN(res_names), n_bfmvar_res-1 /)
        status = nf90_get_var(ncresid, IDname_res, res_names, start=step_start_arr2, count=step_count_arr2)
        if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_NAME var in "//fname_res)

        status = nf90_inq_varid(ncresid, "D3STATE_UNITS", IDunits_res)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_UNITS var in "//fname_res)
        status = nf90_get_var(ncresid, IDunits_res, res_units, start=step_start_arr2, count=step_count_arr2)
        if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_UNITS var in "//fname_res)

        status = nf90_inq_varid(ncresid, "D3STATE_LONG", IDlong_res)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_LONG var in "//fname_res)
        status = nf90_get_var(ncresid, IDlong_res, res_long, start=step_start_arr2, count=step_count_arr2)
        if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_LONG var in "//fname_res)

        allocate(bfmvarid_res(n_bfmvar_res))

        !Define pH variable in last position in output file
        res_names(n_bfmvar_res) = "pH"
        res_units(n_bfmvar_res) = "pH"
        res_long(n_bfmvar_res)  = "acidity or basicity of an aqueous solution"

        do n = 1 , n_bfmvar_res
#ifdef DEBUG
           write(*,*) "Define variable: ",trim(res_names(n))
           write(*,*) "from file ",trim(fname_res)
           write(*,*) "last dimension name ",trim(dimname3d_res)
#endif
           ! check the dimension of the variable
           if (dimname3d_res == "oceanpoint") then
              ! 3D array
              status = nf90_def_var(ncidres, trim(res_names(n)), NF90_REAL, (/ IDx, IDy, IDz, IDtime /), IDtarget)
              if (status /= NF90_NOERR) call handle_err(status)
              bfmvarid_res(n)= IDtarget
           else
              write(*,*) "Variable suppose to be 3D oceanpoint."
              stop
           end if
           ! copy attributes
           status = nf90_put_att(ncidres, IDtarget, "units", trim(res_units(n)))
           if (status /= NF90_NOERR) call handle_err(status)
           status = nf90_put_att(ncidres, IDtarget, "long_name", trim(res_long(n)))
           if (status /= NF90_NOERR) call handle_err(status)
           ! Add fill value
           status = nf90_put_att(ncidres, IDtarget, "_FillValue", NF90_FILL_REAL)
           if (status /= NF90_NOERR) call handle_err(status)
           ! Add reference coordinate names (needed with CDO)
           status = nf90_put_att(ncidres, IDtarget, "coordinates", "lon lat")
           if (status /= NF90_NOERR) call handle_err(status)
           ! inquire variable
           write(*,*) "Variable saved: ",trim(res_names(n))," : ",n
        end do ! n_bfmvar_res

        deallocate(res_names)
        deallocate(res_units)
        deallocate(res_long)
     end if ! 3D restart var
  end if


  ! copy time variable and attributes
  if ( do_output ) then
     status = nf90_inq_varid(ncoutid, "time", IDtimetmp)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring time var in "//fname_out)
     allocate(time_out(ntime_out))
     status = nf90_get_var(ncoutid, IDtimetmp, time_out, start = (/ 1 /), count = (/ ntime_out /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading time values from"//fname_out)
     status = nf90_def_var(ncidout, "time", NF90_REAL, (/ IDtime /), IDvartime_out)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_copy_att(ncoutid, IDtimetmp, "units", ncidout, IDvartime_out)
     if (status /= NF90_NOERR) call handle_err(status,errstring="copying time units")     
  end if
  if ( do_restart ) then
     ! copy time variable and attributes
     status = nf90_inq_varid(ncresid, "time", IDtimetmp)
     allocate(time_res(ntime_res))
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring time var in "//fname_res)
     status = nf90_get_var(ncresid, IDtimetmp, time_res, start = (/ 1 /), count = (/ ntime_res /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading time values from"//fname_res)
     status = nf90_def_var(ncidres, "time", NF90_REAL, (/ IDtime /), IDvartime_res)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_copy_att(ncresid, IDtimetmp, "units", ncidres, IDvartime_res)
     if (status /= NF90_NOERR) call handle_err(status,errstring="copying time units")     
  end if


  ! exit definition mode
  if ( do_output ) then
     status = nf90_enddef(ncidout)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while exiting definition mode in output")
  end if
  if ( do_restart ) then
     status = nf90_enddef(ncidres)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while exiting definition mode in restart")
  end if


  ! write time values
  if ( do_output ) then
     status = nf90_put_var(ncidout, IDvartime_out, time_out, start = (/ 1 /), count = (/ ntime_out /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: time output")
     deallocate(time_out)
  end if
  if ( do_restart ) then
     status = nf90_put_var(ncidres, IDvartime_res, time_res, start = (/ 1 /), count = (/ ntime_res /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: time restart")
     deallocate(time_res)
  end if


  !write fill value
  if ( do_output ) then
     allocate(tmpfillvar3d(jpiglo,jpjglo,jpkglo,ntime_out),tmpfillvar2d(jpiglo,jpjglo,ntime_out))
     tmpfillvar3d = NF90_FILL_REAL
     tmpfillvar2d = NF90_FILL_REAL
     do n = 1 , n_bfmvar_out
        IDvar = bfmvarid_out(n)
        ! inquire ID in the output file
        call handle_err( nf90_inquire_variable(ncoutid, IDvar, ndims=ndims, name=varname) )
        call handle_err( nf90_inq_varid(ncidout, varname, IDtarget) )
        call handle_err( nf90_inquire_variable(ncidout, IDtarget, ndims=ndims) )
        ! write empty variable values
        SELECT CASE (ndims)
        case (4)
           ! 3D array
           call handle_err( nf90_put_var(ncidout, IDtarget, tmpfillvar3d), &
                errstring="variable:"//trim(varname))
        case (3)
           ! 2D array
           call handle_err( nf90_put_var(ncidout, IDtarget, tmpfillvar2d), &
                errstring="variable:"//trim(varname))
        case default
           write(*,'(A,I4)') "invalid dimension size: ",ndims
        end select
     end do
     deallocate(tmpfillvar3d,tmpfillvar2d)
  endif
  if ( do_restart ) then
     allocate(tmpfillvar3d(jpiglo,jpjglo,jpkglo,ntime_res),tmpfillvar2d(jpiglo,jpjglo,ntime_res))
     tmpfillvar3d = NF90_FILL_REAL
     tmpfillvar2d = NF90_FILL_REAL
     do n = 1 , n_bfmvar_res
        IDtarget = bfmvarid_res(n)
        ! inquire ID in the output file
        call handle_err( nf90_inquire_variable(ncidres, IDtarget, ndims=ndims, name=varname) )
        ! write empty variable values
        SELECT CASE (ndims)
        case (4)
           ! 3D array
           call handle_err( nf90_put_var(ncidres, IDtarget, tmpfillvar3d), &
                errstring="variable:"//trim(varname))
        case (3)
           ! 2D array
           call handle_err( nf90_put_var(ncidres, IDtarget, tmpfillvar2d), &
                errstring="variable:"//trim(varname))
        case default
           write(*,'(A,I4,A,A,A,I4)') "invalid dimension size: ",ndims, " for var: '",trim(varname), "'  ID: ", IDtarget
        end select
     end do
     deallocate(tmpfillvar3d,tmpfillvar2d)
  end if

  !close input/output files
  if ( do_output ) then
     status = nf90_close(ncoutid)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while closing BFM input file")
     status = nf90_close(ncidout)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while closing BFM output file")
  end if
  if ( do_restart ) then
     status = nf90_close(ncresid)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while closing BFM restart input file")
     status = nf90_close(ncidres)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while closing BFM restart output file")
  end if

#ifdef DEBUG
  write(*,*) "Output file created"
#endif

  return

end subroutine create_outputfile
