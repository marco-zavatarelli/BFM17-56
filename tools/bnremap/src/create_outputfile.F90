! BFM_NEMO-REMAP bnremap
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

subroutine create_outputfile

  use netcdf
  use mod_bnremap
  implicit none
  integer           :: status
  integer           :: ndims, nVars, nGlobalAtts
  integer           :: IDtime,IDtimetmp,IDunlimdim,IDvartime
  integer           :: IDx, IDy, IDz, IDvar, IDtarget, IDatt
  real, allocatable, dimension(:) :: time
  integer :: dimids(4)
  character(len = NF90_MAX_NAME) :: DimName,varname,attname
  

     status = nf90_create(trim(out_fname), NF90_NOCLOBBER, IDncout)
     if(status /= NF90_NOERR) call handle_err(status,errstring="opening file")
        write(*,*) "Start definition ..."
     ! Define the dimensions
     status = nf90_def_dim(IDncout, "time", NF90_UNLIMITED, IDtime)
     if(status /= NF90_NOERR) call handle_err(status,errstring="defining dim time")
     status = nf90_def_dim(IDncout, "x", jpi, IDx)
     if (status /= NF90_NOERR) call handle_err(status,errstring="defining dim x")
     status = nf90_def_dim(IDncout, "y", jpj, IDy)
     if (status /= NF90_NOERR) call handle_err(status,errstring="defining dim y")
     status = nf90_def_dim(IDncout, "z", jpk, IDz)
     if (status /= NF90_NOERR) call handle_err(status,errstring="defining dim z")
     ! define geographic variables
        write(*,*) "Start defgining variables ..."
              status = nf90_def_var(IDncout, "mask", NF90_REAL, (/ IDx, IDy, IDz /), IDtarget)
              status = nf90_def_var(IDncout, "lat", NF90_REAL, (/ IDx, IDy /), IDtarget)
              status = nf90_def_var(IDncout, "lon", NF90_REAL, (/ IDx, IDy /), IDtarget)
     ! copy global attributes
#ifdef DEBUG
        write(*,*) "Start copying global attributes ..."
#endif
     ! read variables from input file and copy attributes
     status = nf90_inquire(IDncin, nDims, nVars, nGlobalAtts, IDunlimdim)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(IDncin, IDunlimdim, len = ntime)
     if (status /= NF90_NOERR) call handle_err(status)
     do IDatt=1,nGlobalAtts
        status=nf90_inq_attname(IDncin, NF90_GLOBAL, IDatt, name=attname)
        status = nf90_copy_att(IDncin, NF90_GLOBAL, trim(attname), IDncout, NF90_GLOBAL)
        if (status /= NF90_NOERR) call handle_err(status,errstring="copying attribute "//trim(attname))
     end do
     ! keep tracks of the variables that are stored
     allocate(bfmvarid(nVars))
     n_bfmvar = 0
     do IDvar=1,nVars
        ! inquire variable
        status=nf90_inquire_variable(IDncin, IDvar, ndims=ndims, name=varname)
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        status=nf90_inquire_variable(IDncin, IDvar, dimids=dimids(1:ndims))
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        status = nf90_inquire_dimension(IDncin, dimids(ndims), name = DimName)
        if (old_version) then
           status = scan(varname,".")
           if (status > 0) varname = varname(1:2)//varname(4:4)
        end if
#ifdef DEBUG
        write(*,*) "Reading variable ",trim(varname)
        write(*,*) "last dimension name ",trim(DimName)
#endif
        if (DimName /= "time" .OR. ndims == 1) cycle ! enter only with time-varying variables
           ! keep tracks of the variables that are stored
           n_bfmvar = n_bfmvar + 1
           bfmvarid(n_bfmvar) = IDvar
           ! check if it's a 2D or 3D variable
           status = nf90_inquire_dimension(IDncin, dimids(1), name = DimName)
           if (DimName == ocepointname) then
              status = nf90_def_var(IDncout, trim(varname), NF90_REAL, (/ IDx, IDy, IDz, IDtime /), IDtarget)
              if (status /= NF90_NOERR) call handle_err(status)
           else
              status = nf90_def_var(IDncout, trim(varname), NF90_REAL, (/ IDx, IDy, IDtime /), IDtarget)
              if (status /= NF90_NOERR) call handle_err(status)
           end if
           ! copy attributes
           status = nf90_copy_att(IDncin, IDvar, "units", IDncout, IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
           status = nf90_copy_att(IDncin, IDvar, "long_name", IDncout, IDtarget)
           if (status /= NF90_NOERR) call handle_err(status)
     end do ! nvars
     ! copy time variable and attributes
     status = nf90_inq_varid(IDncin, "time", IDtimetmp)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring time var in "//in_fname)
     allocate(timE(ntime))
     status = nf90_get_var(IDncin, IDtimetmp, time, start = (/ 1 /), count = (/ ntime /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="reading time values from"//in_fname)
     status = nf90_def_var(IDncout, "time", NF90_REAL, (/ IDtime /), IDvartime)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_copy_att(IDncin, IDtimetmp, "units", IDncout, IDvartime)
     if (status /= NF90_NOERR) call handle_err(status,errstring="copying time units")     
     ! exit definition mode
     status = nf90_enddef(IDncout)
     if (status /= NF90_NOERR) call handle_err(status,errstring="while exiting definition mode")
     ! write time values
     status = nf90_put_var(IDncout, IDvartime, time, start = (/ 1 /), count = (/ ntime /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: time")
#ifdef DEBUG
        write(*,*) "Output file created!"
#endif

     return

end subroutine create_outputfile
