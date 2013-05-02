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

  use netcdf
  use mod_bnmerge
  implicit none
  integer           :: ncid, ncbfmid, status
  integer           :: p,n,d,t
  character(LEN=172) :: fname
  integer           :: ndims, nVars, nGlobalAtts
  integer           :: IDtime,IDtimetmp,IDunlimdim,IDvartime
  integer           :: IDx, IDy, IDz, IDvar, IDtarget, IDatt
  integer           :: ntime
  real, allocatable, dimension(:) :: time
  integer           :: npoints
  integer :: vartype,dimids(4),dimlen(4)
  character(len = NF90_MAX_NAME) :: DimName,varname,attname
  

     status = nf90_create(trim(out_dir)//"/"//trim(chunk_fname)//".nc", NF90_NOCLOBBER, ncid)
     if(status /= NF90_NOERR) call handle_err(status,errstring="A file named "//trim(chunk_fname)//".nc already exists!" )
     ! Define the dimensions
     status = nf90_def_dim(ncid, "time", NF90_UNLIMITED, IDtime)
     if(status /= NF90_NOERR) call handle_err(status)
     status = nf90_def_dim(ncid, "x", jpiglo, IDx)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_def_dim(ncid, "y", jpjglo, IDy)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_def_dim(ncid, "depth", jpk, IDz)
     if (status /= NF90_NOERR) call handle_err(status)
     ! read variables from domain 0000 and copy attributes
     fname = trim(inp_dir)//"/"//trim(chunk_fname)//"_0000.nc" 
     status = nf90_open(path = fname, mode = NF90_WRITE, ncid = ncbfmid)
     if (status /= NF90_NOERR) call handle_err(status)     
     status = nf90_inquire(ncbfmid, nDims, nVars, nGlobalAtts, IDunlimdim)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(ncbfmid, IDunlimdim, len = ntime)
     if (status /= NF90_NOERR) call handle_err(status)
     ! define geographic variables
     if (ln_mask) then
       call handle_err (nf90_def_var(ncid, "mask", NF90_REAL, (/ IDx, IDy, IDz /), IDtarget))
       call handle_err (nf90_put_att(ncid, IDtarget, "coordinates", "lon lat"))
     endif
     call handle_err (nf90_def_var(ncid, "lat", NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err (nf90_put_att(ncid, IDtarget, "units", "degrees_north"))
     call handle_err (nf90_def_var(ncid, "lon", NF90_REAL, (/ IDx, IDy /), IDtarget))
     call handle_err (nf90_put_att(ncid, IDtarget, "units", "degrees_east"))
     call handle_err (nf90_def_var(ncid, "depth", NF90_REAL, (/ IDz /), IDtarget))
     call handle_err (nf90_put_att(ncid, IDtarget, "long_name", "depth_below_sea"))
     call handle_err (nf90_put_att(ncid, IDtarget, "units", "m"))
     call handle_err (nf90_put_att(ncid, IDtarget, "positive", "down"))
     call handle_err (nf90_put_att(ncid, IDtarget, "axis", "Z"))
     ! copy global attributes
#ifdef DEBUG
        write(*,*) "Creating file:",trim(chunk_fname)//".nc"," containing",ntime,"time frames"
        write(*,*) "Start copying global attributes ..."
#endif
     do IDatt=1,nGlobalAtts
        status=nf90_inq_attname(ncbfmid, NF90_GLOBAL, IDatt, name=attname)
        status = nf90_copy_att(ncbfmid, NF90_GLOBAL, trim(attname), ncid, NF90_GLOBAL)
        if (status /= NF90_NOERR) call handle_err(status,errstring="copying attribute "//trim(attname))
     end do
     ! Tracks of the variables that have to be stored
     allocate(bfmvarid(nVars))
     n_bfmvar=0
     do IDvar=1,nVars
       status=nf90_inquire_variable(ncbfmid, IDvar, ndims=ndims, name=varname)
       if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
       do n = 1 , NSAVE
           if ( trim(var_save(n)) == trim(varname) ) then
              n_bfmvar = n_bfmvar + 1
              bfmvarid(n_bfmvar)= IDvar
              write(*,*) "Assigned output ",trim(varname), " with ID:", IDvar
           endif
       enddo 
     enddo
     write(*,*) "Total Output Variables ", n_bfmvar
     write(*,*)
     if (n_bfmvar == 0) then
        write(*,*) "Selected output variables do not match the content of the input files.", n_bfmvar
        stop
     endif

     ! Assign dimensions and attributes to variables 
     do n = 1 , n_bfmvar
        IDvar = bfmvarid(n)
        ! inquire variable
        status=nf90_inquire_variable(ncbfmid, IDvar, ndims=ndims, name=varname)
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        status=nf90_inquire_variable(ncbfmid, IDvar, dimids=dimids(1:ndims))
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        status = nf90_inquire_dimension(ncbfmid, dimids(ndims), name = DimName)
        write(*,*) "Define variable: ",trim(varname)," with ID: ",IDvar
#ifdef DEBUG
        write(*,*) "from file ",trim(fname)
        write(*,*) "last dimension name ",trim(DimName)
#endif
        if (DimName /= "time" .OR. ndims == 1) cycle ! enter only with time-varying variables
           ! check the dimension of the variable
           status = nf90_inquire_dimension(ncbfmid, dimids(1), name = DimName)
           if (DimName == "oceanpoint") then
           ! 3D array
              status = nf90_def_var(ncid, trim(varname), NF90_REAL, (/ IDx, IDy, IDz, IDtime /), IDtarget)
              if (status /= NF90_NOERR) call handle_err(status)
           else
           ! 2D array
              status = nf90_def_var(ncid, trim(varname), NF90_REAL, (/ IDx, IDy, IDtime /), IDtarget)
              if (status /= NF90_NOERR) call handle_err(status)
           end if
           ! copy attributes
           status = nf90_copy_att(ncbfmid, IDvar, "units", ncid, IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
           status = nf90_copy_att(ncbfmid, IDvar, "long_name", ncid, IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
           ! Add fill value
           status = nf90_put_att(ncid, IDtarget, "_FillValue", NF90_FILL_REAL)
           if (status /= NF90_NOERR) call handle_err(status)
           ! Add reference coordinate names (needed with CDO)
           status = nf90_put_att(ncid, IDtarget, "coordinates", "lon lat")
           if (status /= NF90_NOERR) call handle_err(status)
     end do ! n_bfmvar 

     ! copy time variable and attributes
     status = nf90_inq_varid(ncbfmid, "time", IDtimetmp)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring time var in "//fname)
     allocate(time(ntime))
     status = nf90_get_var(ncbfmid, IDtimetmp, time, start = (/ 1 /), count = (/ ntime /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading time values from"//fname)
     status = nf90_def_var(ncid, "time", NF90_REAL, (/ IDtime /), IDvartime)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_copy_att(ncbfmid, IDtimetmp, "units", ncid, IDvartime)
     if (status /= NF90_NOERR) call handle_err(status,errstring="copying time units")     
     ! close the 0000 BFM netcdf file
     status = nf90_close(ncbfmid)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while closing BFM input file")
     ! exit definition mode
     status = nf90_enddef(ncid)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while exiting definition mode")
     ! write time values
     status = nf90_put_var(ncid, IDvartime, time, start = (/ 1 /), count = (/ ntime /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: time")
     status = nf90_close(ncid)
     if (status /= NF90_NOERR) call handle_err(status)
#ifdef DEBUG
        write(*,*) "Output file created with ",ntime,"time frames"
#endif

     return

end subroutine create_outputfile
