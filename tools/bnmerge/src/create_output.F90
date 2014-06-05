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
module create_output

  use netcdf
  use mod_bnmerge, ONLY : handle_err, RLEN
  use mod_bnmerge, ONLY : nc_compres,nc_shuffle,nc_deflate,nc_defllev
  implicit none

#ifdef PARAL
  include 'mpif.h'
  integer :: nthread=0
#endif

  ! !PUBLIC MEMBER FUNCTIONS:
  public create_output_init

contains

  subroutine create_output_init
    use netcdf
    use mod_bnmerge, ONLY: chunk_fname, bfm_restart, out_fname, do_restart, do_output, n_bfmvar_out, n_bfmvar_res

    implicit none
    integer :: ncidout, ncidres, ncid_chunk_out, ncid_chunk_res
    integer :: id_array_out(4), id_array_res(4)
    integer :: IDvartime_out, IDvartime_res
    integer :: ntime_out, ntime_res

    ! Create file and attribute and dimensions
    if ( do_output  ) call define_output(out_fname, chunk_fname, ncidout, ncid_chunk_out, id_array_out, IDvartime_out, ntime_out)
    if ( do_restart ) call define_output(bfm_restart, bfm_restart, ncidres, ncid_chunk_res, id_array_res, IDvartime_res, ntime_res)

    write(*,*) "Writing variables..."
    if( do_output  ) call define_variables_out(ncid_chunk_out, ncidout, id_array_out, n_bfmvar_out )
    if( do_restart ) call define_variables_res(ncid_chunk_res, ncidres, id_array_res, n_bfmvar_res )

    write(*,*) "Exit definition mode"
    ! exit definition mode
    if ( do_output  ) call handle_err(nf90_enddef(ncidout),errstring="while exiting definition mode in output")
    if ( do_restart ) call handle_err(nf90_enddef(ncidres),errstring="while exiting definition mode in restart")

    write(*,*) "Filling values"
    ! write time values and fill values
    if( do_output  ) call fill_values(ncid_chunk_out, ncidout, IDvartime_out, ntime_out, n_bfmvar_out)
    if( do_restart ) call fill_values(ncid_chunk_res, ncidres, IDvartime_res, ntime_res, n_bfmvar_res)

    write(*,*) "Closing files"
    !close input/output files
    if ( do_output ) then
       call handle_err(nf90_close(ncid_chunk_out),errstring="while closing BFM input file")
       call handle_err(nf90_close(ncidout),errstring="while closing BFM output file")
       write(*,*) "Output merged file created"
    end if
    if ( do_restart ) then
       call handle_err(nf90_close(ncid_chunk_res),errstring="while closing BFM restart input file")
       call handle_err(nf90_close(ncidres),errstring="while closing BFM restart output file")
       write(*,*) "Restart merged file created"
    end if

  end subroutine create_output_init


  subroutine define_output(fname, chunkname, ncid_out, ncid_chunk, id_array, id_vartime, ntime)
    use netcdf
    use mod_bnmerge, ONLY : jpiglo, jpjglo, jpkglo, inp_dir, out_dir, ln_mask

    implicit none
    character(len=100), intent(in)  :: fname, chunkname
    integer,            intent(out) :: ncid_out, ncid_chunk, id_array(4), id_vartime, ntime

    integer                        :: status
    integer                        :: nGlobalAtts 
    integer                        :: IDx, IDy, IDz, IDtime, IDatt, IDunlimdim, IDtimetmp, IDocepnt, IDsrfpnt, IDbtnpnt
    character(len = NF90_MAX_NAME) :: globattname


#ifdef PARAL
    status = nf90_create_par(path=trim(out_dir)//"/"//trim(fname)//".nc", &
         cmode=IOR(NF90_NETCDF4, NF90_MPIIO), &
         ncid=ncid_out, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
#else
    status = nf90_create(trim(out_dir)//"/"//trim(fname)//".nc", NF90_NETCDF4, ncid_out)
#endif   
    if(status /= NF90_NOERR) call handle_err(status,errstring="A file named "//trim(fname)//".nc already exists!" )

    ! Define the dimensions
    status = nf90_def_dim(ncid_out, "x", jpiglo, IDx)
    if (status /= NF90_NOERR) call handle_err(status,errstring="reading X in out")
    status = nf90_def_dim(ncid_out, "y", jpjglo, IDy)
    if (status /= NF90_NOERR) call handle_err(status,errstring="reading Y in out")
    status = nf90_def_dim(ncid_out, "depth", jpkglo, IDz)
    if (status /= NF90_NOERR) call handle_err(status,errstring="reading depth in out")
    status = nf90_def_dim(ncid_out, "time", NF90_UNLIMITED, IDtime)
    if(status /= NF90_NOERR) call handle_err(status,errstring="reading time in out")

    ! read variables from domain 0000 and copy attributes
#ifdef PARAL
    status = nf90_open_par(path = trim(inp_dir)//"/"//trim(chunkname)//"_0000.nc", &
         cmode=IOR(NF90_NOWRITE, NF90_MPIIO), &
         ncid = ncid_chunk, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
#else
    status = nf90_open(path = trim(inp_dir)//"/"//trim(chunkname)//"_0000.nc", mode = NF90_NOWRITE, ncid = ncid_chunk)
#endif
    if (status /= NF90_NOERR) call handle_err(status,errstring="opening named "//trim(chunkname)//"_0000.nc!" )
    status = nf90_inquire(ncid_chunk, nAttributes=nGlobalAtts, unlimitedDimId=IDunlimdim)
    if (status /= NF90_NOERR) call handle_err(status,errstring="reading info")
    status = nf90_inquire_dimension(ncid_chunk, IDunlimdim, len = ntime)
    if (status /= NF90_NOERR) call handle_err(status,errstring="reading time")

    ! define geographic variables
    if ( ln_mask ) then
       call handle_err (nf90_def_var(ncid_out, "mask", NF90_REAL, (/ IDx, IDy, IDz /), IDatt))
       if ( nc_compres )  &
          call handle_err( nf90_def_var_deflate(ncid_out, IDatt, nc_shuffle, nc_deflate, nc_defllev) )
       call handle_err (nf90_put_att(ncid_out, IDatt, "coordinates", "lon lat"))
    end if
    call handle_err (nf90_def_var(ncid_out, "lat", NF90_REAL, (/ IDx, IDy /), IDatt))
    call handle_err( nf90_put_att(ncid_out, IDatt, "_FillValue", NF90_FILL_REAL) )
    call handle_err (nf90_put_att(ncid_out, IDatt, "units", "degrees_north"))
    call handle_err (nf90_def_var(ncid_out, "lon", NF90_REAL, (/ IDx, IDy /), IDatt))
    call handle_err( nf90_put_att(ncid_out, IDatt, "_FillValue", NF90_FILL_REAL) )
    call handle_err (nf90_put_att(ncid_out, IDatt, "units", "degrees_east"))
    call handle_err (nf90_def_var(ncid_out, "depth", NF90_REAL, (/ IDz /), IDatt))
    call handle_err( nf90_put_att(ncid_out, IDatt, "_FillValue", NF90_FILL_REAL) )
    call handle_err (nf90_put_att(ncid_out, IDatt, "long_name", "depth_below_sea"))
    call handle_err (nf90_put_att(ncid_out, IDatt, "units", "m"))
    call handle_err (nf90_put_att(ncid_out, IDatt, "positive", "down"))
    call handle_err (nf90_put_att(ncid_out, IDatt, "axis", "Z"))

    ! copy global attributes
#ifdef DEBUG
    write(*,*) "Creating file:",trim(fname)//".nc"," containing",ntime,"time frames"
    write(*,*) "Start copying global attributes ..."
#endif
    do IDatt=1,nGlobalAtts
       status=nf90_inq_attname(ncid_chunk, NF90_GLOBAL, IDatt, name=globattname)
       status = nf90_copy_att(ncid_chunk, NF90_GLOBAL, trim(globattname), ncid_out, NF90_GLOBAL)
       if (status /= NF90_NOERR) call handle_err(status,errstring="copying attribute "//trim(globattname))
    end do

    ! define time variable and attributes
#ifdef DEBUG
    write(*,*) "define time variable and attributes ..."
#endif
    status = nf90_inq_varid(ncid_chunk, "time", IDtimetmp)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring time var in "//chunkname)
    status = nf90_def_var(ncid_out, "time", NF90_REAL, (/ IDtime /), id_vartime)
    if (status /= NF90_NOERR) call handle_err(status)
    status = nf90_copy_att(ncid_chunk, IDtimetmp, "units", ncid_out, id_vartime)
    if (status /= NF90_NOERR) call handle_err(status,errstring="copying time units")

    id_array = (/ IDx, IDy, IDz, IDtime /)

  end subroutine define_output


  subroutine fill_values(ncid_chunk, ncid_out, id_vartime, ntime, nvars)
    use netcdf
    use mod_bnmerge, ONLY : jpiglo, jpjglo, jpkglo

    implicit none

    integer, intent(in) :: ncid_chunk, ncid_out, id_vartime, ntime, nvars

    real, allocatable :: time(:)
    integer           :: IDtimetmp, status

    integer                        :: IDtarget, ndims
    character(len = NF90_MAX_NAME) :: varname
    real(RLEN), allocatable        :: tmpfillvar3d(:,:,:,:), tmpfillvar2d(:,:,:)

    allocate(time(ntime))
    status = nf90_inq_varid(ncid_chunk, "time", IDtimetmp)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring time var")
    status = nf90_get_var(ncid_chunk, IDtimetmp, time, start = (/ 1 /), count = (/ ntime /))
    if (status /= NF90_NOERR) call handle_err(status,errstring="reading time values from")
    status = nf90_put_var(ncid_out, id_vartime, time, start = (/ 1 /), count = (/ ntime /))
    if (status /= NF90_NOERR) call handle_err(status,errstring="variable: time output")
    deallocate(time)

    !write fill value
    allocate(tmpfillvar3d(jpiglo,jpjglo,jpkglo,ntime),tmpfillvar2d(jpiglo,jpjglo,ntime))
    tmpfillvar3d = NF90_FILL_DOUBLE
    tmpfillvar2d = NF90_FILL_DOUBLE
    do IDtarget = 1 , nvars
       ! inquire ID in the output file
       call handle_err( nf90_inquire_variable(ncid_out, IDtarget, ndims=ndims, name=varname), &
               errstring="inquiring target in out variable" )
       ! write empty variable values
       select case (ndims)
       case (4)
          ! 3D array
          call handle_err( nf90_put_var(ncid_out, IDtarget, tmpfillvar3d), &
               errstring="variable:"//trim(varname))
       case (3)
          ! 2D array
          call handle_err( nf90_put_var(ncid_out, IDtarget, tmpfillvar2d), &
               errstring="variable:"//trim(varname))
          ! case default
          !   write(*,'(A,I4,A,A,A,I4)') "invalid dimension size: ",ndims, " for var: '",trim(varname), "'  ID: ", IDtarget
       end select
    end do
    deallocate(tmpfillvar3d,tmpfillvar2d)
  end subroutine fill_values


  subroutine define_variables_out(ncid_chunk, ncidout, id_array_out, n_var_total )
    use netcdf
    use mod_bnmerge, ONLY: replace_char
    use mod_bnmerge, ONLY: NSAVE, TYPE_OCE, TYPE_BTN, TYPE_SRF, var_save, bfmvarid_out, bfmvartype_out, bfmvartarget_out

    implicit none
    integer, intent(in)  :: ncid_chunk, ncidout, id_array_out(4)
    integer, intent(out) :: n_var_total

    integer                        :: n, nvars, IDvar, IDtarget, status, ndims, dimids_out(4)
    character(len = NF90_MAX_NAME) :: varname, dimname_out, string


    ! Tracks of the variables that have to be stored
    status = nf90_inquire(ncid_chunk, nVariables=nvars)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring numnber of variables in OUTPUT ")
    allocate(bfmvarid_out(nvars))
    n_var_total=0
    do IDvar=1,nvars
       status=nf90_inquire_variable(ncid_chunk, IDvar, name=varname)
       if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
       do n = 1 , NSAVE
          string = var_save(n)
          call replace_char(str=string, tar='()', rep='_')
          if ( ( trim(var_save(n)) == trim(varname) ) .OR. ( trim(string) == trim(varname) ) ) then
             n_var_total = n_var_total + 1
             bfmvarid_out(n_var_total)   = IDvar
          endif
       enddo
    enddo
    write(*,*) "Total Output Variables  ", n_var_total
    if (n_var_total == 0) then
       write(*,*) "Selected output variables do not match the content of the input files.", n_var_total
       stop
    endif

    ! Assign dimensions and attributes to variables 
    allocate(bfmvartype_out(n_var_total))
    allocate(bfmvartarget_out(n_var_total))
    do n = 1 , n_var_total
       IDvar = bfmvarid_out(n)
       ! inquire variable
       status=nf90_inquire_variable(ncid_chunk, IDvar, ndims=ndims, name=varname)
       if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
       status=nf90_inquire_variable(ncid_chunk, IDvar, dimids=dimids_out(1:ndims))
       if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
       status = nf90_inquire_dimension(ncid_chunk, dimids_out(ndims), name = dimname_out)
       if (dimname_out /= "time" .OR. ndims == 1) cycle ! enter only with time-varying variables
       ! check the dimension of the variable
       status = nf90_inquire_dimension(ncid_chunk, dimids_out(1), name = dimname_out)
#ifdef DEBUG
       write(*,*) "Define variable: ",trim(varname)
       write(*,*) "from file chunk in output file"
       write(*,*) "last dimension name ",trim(dimname_out)
#endif
       select case(dimname_out)
       case ("oceanpoint") ! 3D array
          bfmvartype_out(n) = TYPE_OCE
          status = nf90_def_var(ncidout, trim(varname), NF90_DOUBLE, &
               id_array_out, IDtarget)
       case("surfacepoint") ! 2D surface
          bfmvartype_out(n) = TYPE_SRF
          status = nf90_def_var(ncidout, trim(varname), NF90_DOUBLE, &
               (/ id_array_out(1), id_array_out(2), id_array_out(4) /), IDtarget)
       case("bottompoint") ! 2D bottom
          bfmvartype_out(n) = TYPE_BTN
          status = nf90_def_var(ncidout, trim(varname), NF90_DOUBLE, &
               (/ id_array_out(1), id_array_out(2), id_array_out(4) /), IDtarget)
       case default
          write(*,*) "Dimension name not recognized: "//trim(dimname_out)
          stop
       end select
       if ( nc_compres )  &
          call handle_err( nf90_def_var_deflate(ncidout, IDtarget, nc_shuffle, nc_deflate, nc_defllev) )
       bfmvartarget_out(n) = IDtarget

       ! copy attributes
       status = nf90_copy_att(ncid_chunk, IDvar, "units", ncidout, IDtarget)
       if (status /= NF90_NOERR) call handle_err(status)
       status = nf90_copy_att(ncid_chunk, IDvar, "long_name", ncidout, IDtarget)
       if (status /= NF90_NOERR) call handle_err(status)
       ! Add fill value
       status = nf90_put_att(ncidout, IDtarget, "_FillValue", NF90_FILL_DOUBLE)
       if (status /= NF90_NOERR) call handle_err(status)
       ! Add reference coordinate names (needed with CDO)
       status = nf90_put_att(ncidout, IDtarget, "coordinates", "lon lat")
       if (status /= NF90_NOERR) call handle_err(status)
    end do ! n_var_total 

  end subroutine define_variables_out


  subroutine define_variables_res(ncid_chunk, ncidres, id_array_res, n_var_total)
    use netcdf
    use mod_bnmerge, ONLY: TYPE_OCE, TYPE_BTN, TYPE_SRF, TYPE_PH, bfmvarid_res, bfmvartype_res, bfmvartarget_res

    implicit none
    integer, intent(in)  :: ncid_chunk, ncidres, id_array_res(4)
    integer, intent(out) :: n_var_total

    integer                        :: n, IDvar, IDtarget, IDph, status
    integer                        :: nbfmvars3d_res, nbfmvars2dben_res, nbfmvars2dice_res
    character(len=64), allocatable :: res_names(:),res_units(:),res_long(:)
    integer                        :: step_start_arr2(2), step_count_arr2(2)
    integer                        :: dimids_res(2)
    character(len = NF90_MAX_NAME) :: dimname_res

    n_var_total = 0

    !Get the total number of variables
    nbfmvars3d_res    = 0
    nbfmvars2dben_res = 0
    nbfmvars2dice_res = 0
    ! for 3d variables
    status = nf90_inq_dimid(ncid_chunk, "d3vars", IDvar)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring d3vars")
    status = nf90_inquire_dimension(ncid_chunk, IDvar, len=nbfmvars3d_res)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring number of 3D vars")
    ! for pH variable
    nbfmvars3d_res = nbfmvars3d_res + 1
    ! for 2d benthic variables
    status = nf90_inq_dimid(ncid_chunk, "d2vars_ben", IDvar)
    if (status == NF90_NOERR) then
       status = nf90_inquire_dimension(ncid_chunk, IDvar, len=nbfmvars2dben_res)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring number of 2D ben vars")
    end if
    ! for 2d ice variables
    status = nf90_inq_dimid(ncid_chunk, "d2vars_ice", IDvar)
    if (status == NF90_NOERR) then
       status = nf90_inquire_dimension(ncid_chunk, IDvar, len=nbfmvars2dice_res)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring number of 2D ice vars")
    end if

    !allocate global arrays
    allocate(bfmvarid_res(nbfmvars3d_res+nbfmvars2dben_res+nbfmvars2dice_res))
    allocate(bfmvartype_res(nbfmvars3d_res+nbfmvars2dben_res+nbfmvars2dice_res))
    allocate(bfmvartarget_res(nbfmvars3d_res+nbfmvars2dben_res+nbfmvars2dice_res))
   
    ! define 3d variables (always present)
    allocate(res_names(nbfmvars3d_res))
    allocate(res_units(nbfmvars3d_res))
    allocate(res_long(nbfmvars3d_res))
    step_start_arr2 = (/ 1, 1 /)
    step_count_arr2 = (/ LEN(res_names), nbfmvars3d_res-1 /)

    status = nf90_inq_varid(ncid_chunk, "D3STATE", IDvar)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE var")
    status=nf90_inquire_variable(ncid_chunk, IDvar, dimids=dimids_res)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE dims")
    status = nf90_inquire_dimension(ncid_chunk, dimids_res(2), name=dimname_res)
    if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE dims")

    status = nf90_inq_varid(ncid_chunk, "D3STATE_NAME", IDvar)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_NAME var")
    status = nf90_get_var(ncid_chunk, IDvar, res_names, start=step_start_arr2, count=step_count_arr2)
    if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_NAME var")

    status = nf90_inq_varid(ncid_chunk, "D3STATE_UNITS", IDvar)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_UNITS var")
    status = nf90_get_var(ncid_chunk, IDvar, res_units, start=step_start_arr2, count=step_count_arr2)
    if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_UNITS var")

    status = nf90_inq_varid(ncid_chunk, "D3STATE_LONG", IDVar)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_LONG var")
    status = nf90_get_var(ncid_chunk, IDvar, res_long, start=step_start_arr2, count=step_count_arr2)
    if (status /= NF90_NOERR) call handle_err(status,errstring="getting D3STATE_LONG var")

    !Define pH variable in last position
    status = nf90_inq_varid(ncid_chunk, "pH", IDph)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring pH var")
    res_names(nbfmvars3d_res) = "pH"
    res_units(nbfmvars3d_res) = "pH"
    res_long(nbfmvars3d_res)  = "acidity or basicity of an aqueous solution"

    do n = 1 , nbfmvars3d_res
#ifdef DEBUG
       write(*,*) "Define variable "//trim(res_names(n))//": ",n
       write(*,*) "from file chunk in RESTART"
       write(*,*) "last dimension name ",trim(dimname_res)
#endif
       ! check the dimension of the variable
       if (dimname_res == "oceanpoint") then ! 3D array
          status = nf90_def_var(ncidres, trim(res_names(n)), NF90_DOUBLE, id_array_res, IDtarget)
          if (status /= NF90_NOERR) call handle_err(status)
          if ( nc_compres )  &
            call handle_err( nf90_def_var_deflate(ncidres, IDtarget, nc_shuffle, nc_deflate, nc_defllev) ) 
       else
          write(*,*) "Variable suppose to be 3D oceanpoint."
          stop
       end if
       ! add attributes
       status = nf90_put_att(ncidres, IDtarget, "units", trim(res_units(n)))
       if (status /= NF90_NOERR) call handle_err(status)
       status = nf90_put_att(ncidres, IDtarget, "long_name", trim(res_long(n)))
       if (status /= NF90_NOERR) call handle_err(status)
       ! Add fill value
       status = nf90_put_att(ncidres, IDtarget, "_FillValue", NF90_FILL_DOUBLE)
       if (status /= NF90_NOERR) call handle_err(status)
       ! Add reference coordinate names (needed with CDO)
       status = nf90_put_att(ncidres, IDtarget, "coordinates", "lon lat")
       if (status /= NF90_NOERR) call handle_err(status)
       !fill the global array for merging
       n_var_total = n_var_total + 1
       if( res_names(n) == "pH" ) then
          bfmvarid_res(n_var_total)   = IDph 
          bfmvartype_res(n_var_total) = TYPE_PH
       else
          bfmvarid_res(n_var_total)   = n
          bfmvartype_res(n_var_total) = TYPE_OCE
       endif
       bfmvartarget_res(n_var_total) = IDtarget
    end do ! nbfmvars3d_res
    write(*,*) "Total Restart 3D Variables ", nbfmvars3d_res

    deallocate(res_names)
    deallocate(res_units)
    deallocate(res_long)


    ! define 2d benthic variables (optional)
    if (nbfmvars2dben_res > 0) then
       allocate(res_names(nbfmvars2dben_res))
       allocate(res_units(nbfmvars2dben_res))
       allocate(res_long(nbfmvars2dben_res))
       step_start_arr2 = (/ 1, 1 /)
       step_count_arr2 = (/ LEN(res_names), nbfmvars2dben_res /)

       status = nf90_inq_varid(ncid_chunk, "D2STATE_BEN", IDvar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_BEN var")
       status=nf90_inquire_variable(ncid_chunk, IDvar, dimids=dimids_res)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_BEN dims")
       status = nf90_inquire_dimension(ncid_chunk, dimids_res(2), name=dimname_res)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_BEN dims")

       status = nf90_inq_varid(ncid_chunk, "D2STATE_BEN_NAME", IDvar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_BEN_NAME var")
       status = nf90_get_var(ncid_chunk, IDvar, res_names, start=step_start_arr2, count=step_count_arr2)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_BEN_NAME var")

       status = nf90_inq_varid(ncid_chunk, "D2STATE_BEN_UNITS", IDvar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_BEN_UNITS var")
       status = nf90_get_var(ncid_chunk, IDvar, res_units, start=step_start_arr2, count=step_count_arr2)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_BEN_UNITS var")

       status = nf90_inq_varid(ncid_chunk, "D2STATE_BEN_LONG", IDVar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_BEN_LONG var")
       status = nf90_get_var(ncid_chunk, IDvar, res_long, start=step_start_arr2, count=step_count_arr2)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_BEN_LONG var")

       do n = 1 , nbfmvars2dben_res
#ifdef DEBUG
          write(*,*) "Define variable: ",trim(res_names(n))
          write(*,*) "from file chunk in RESTART"
          write(*,*) "last dimension name ",trim(dimname_res)
#endif
          ! check the dimension of the variable
          if (dimname_res == "bottompoint") then ! 2D array
             status = nf90_def_var(ncidres, trim(res_names(n)), NF90_DOUBLE, &
                  (/ id_array_res(1), id_array_res(2), id_array_res(4) /), IDtarget)
             if (status /= NF90_NOERR) call handle_err(status)
          if ( nc_compres )  &
            call handle_err( nf90_def_var_deflate(ncidres, IDtarget, nc_shuffle, nc_deflate, nc_defllev) )
          else
             write(*,*) "Variable suppose to be 2D bottompoint."
             stop
          end if
          ! add attributes
          status = nf90_put_att(ncidres, IDtarget, "units", trim(res_units(n)))
          if (status /= NF90_NOERR) call handle_err(status)
          status = nf90_put_att(ncidres, IDtarget, "long_name", trim(res_long(n)))
          if (status /= NF90_NOERR) call handle_err(status)
          ! Add fill value
          status = nf90_put_att(ncidres, IDtarget, "_FillValue", NF90_FILL_DOUBLE)
          if (status /= NF90_NOERR) call handle_err(status)
          ! Add reference coordinate names (needed with CDO)
          status = nf90_put_att(ncidres, IDtarget, "coordinates", "lon lat")
          if (status /= NF90_NOERR) call handle_err(status)
          !fill the global array for merging
          n_var_total = n_var_total + 1
          bfmvarid_res(n_var_total)   = n
          bfmvartype_res(n_var_total) = TYPE_BTN
          bfmvartarget_res(n_var_total) = IDtarget
       end do ! nbfmvars2dben_res
       write(*,*) "Total Restart Benthic Variables ", nbfmvars2dben_res

       deallocate(res_names)
       deallocate(res_units)
       deallocate(res_long)
    end if ! 2D restart ben var

    ! define 2d ice variables (optional)
    if (nbfmvars2dice_res > 0) then
       allocate(res_names(nbfmvars2dice_res))
       allocate(res_units(nbfmvars2dice_res))
       allocate(res_long(nbfmvars2dice_res))
       step_start_arr2 = (/ 1, 1 /)
       step_count_arr2 = (/ LEN(res_names), nbfmvars2dice_res /)

       status = nf90_inq_varid(ncid_chunk, "D2STATE_ICE", IDvar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_ICE var")
       status=nf90_inquire_variable(ncid_chunk, IDvar, dimids=dimids_res)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_ICE dims")
       status = nf90_inquire_dimension(ncid_chunk, dimids_res(2), name=dimname_res)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_ICE dims")

       status = nf90_inq_varid(ncid_chunk, "D2STATE_ICE_NAME", IDvar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_ICE_NAME var")
       status = nf90_get_var(ncid_chunk, IDvar, res_names, start=step_start_arr2, count=step_count_arr2)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_ICE_NAME var")

       status = nf90_inq_varid(ncid_chunk, "D2STATE_ICE_UNITS", IDvar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_ICE_UNITS var")
       status = nf90_get_var(ncid_chunk, IDvar, res_units, start=step_start_arr2, count=step_count_arr2)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_ICE_UNITS var")

       status = nf90_inq_varid(ncid_chunk, "D2STATE_ICE_LONG", IDVar)
       if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D2STATE_ICE_LONG var")
       status = nf90_get_var(ncid_chunk, IDvar, res_long, start=step_start_arr2, count=step_count_arr2)
       if (status /= NF90_NOERR) call handle_err(status,errstring="getting D2STATE_ICE_LONG var")

       do n = 1 , nbfmvars2dice_res
#ifdef DEBUG
          write(*,*) "Define variable: ",trim(res_names(n))
          write(*,*) "from file chunk in RESTART"
          write(*,*) "last dimension name ",trim(dimname_res)
#endif
          ! check the dimension of the variable
          if (dimname_res == "surfacepoint") then ! 2D array
             status = nf90_def_var(ncidres, trim(res_names(n)), NF90_DOUBLE, &
                  (/ id_array_res(1), id_array_res(2), id_array_res(4) /), IDtarget)
             if (status /= NF90_NOERR) call handle_err(status)
          if ( nc_compres )  &
            call handle_err( nf90_def_var_deflate(ncidres, IDtarget, nc_shuffle, nc_deflate, nc_defllev) )
          else
             write(*,*) "Variable suppose to be 2D surfacepoint."
             stop
          end if
          ! add attributes
          status = nf90_put_att(ncidres, IDtarget, "units", trim(res_units(n)))
          if (status /= NF90_NOERR) call handle_err(status)
          status = nf90_put_att(ncidres, IDtarget, "long_name", trim(res_long(n)))
          if (status /= NF90_NOERR) call handle_err(status)
          ! Add fill value
          status = nf90_put_att(ncidres, IDtarget, "_FillValue", NF90_FILL_DOUBLE)
          if (status /= NF90_NOERR) call handle_err(status)
          ! Add reference coordinate names (needed with CDO)
          status = nf90_put_att(ncidres, IDtarget, "coordinates", "lon lat")
          if (status /= NF90_NOERR) call handle_err(status)
          !fill the global array for merging
          n_var_total = n_var_total + 1
          bfmvarid_res(n_var_total)   = n
          bfmvartype_res(n_var_total) = TYPE_SRF
          bfmvartarget_res(n_var_total) = IDtarget
       end do ! nbfmvars2dice_res
       write(*,*) "Total Restart Ice Variables ", nbfmvars2dice_res

       deallocate(res_names)
       deallocate(res_units)
       deallocate(res_long)
    end if ! 2D restart ice var

  end subroutine define_variables_res

end module create_output
