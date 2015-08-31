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
subroutine remap_vars
  use netcdf
  use mod_bnremap
  implicit none
  integer           :: i,j,k,d,status
  integer           :: ndims, nVars, nGlobalAtts, IDunlimdim
  integer           :: IDvar,IDocepnt
  integer           :: lenoce
  integer           :: noce
  real, allocatable, dimension(:,:) :: chunk
  real, allocatable, dimension(:,:,:,:) :: bfmvar3d
  real, allocatable, dimension(:,:,:)   :: bfmvar2d
  integer :: dimids(4),dimlen(4)
  character(len = NF90_MAX_NAME) :: varname


     status = nf90_inquire(IDncin, nDims, nVars, nGlobalAtts, IDunlimdim)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(IDncin, IDunlimdim, len = ntime)
     ! inquire 1D dimension
     status = nf90_inq_dimid(IDncin, ocepointname, IDocepnt)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(IDncin, IDocepnt, len = lenoce)

     ! allocate local 3D and 2D variables
     allocate(bfmvar3d(jpi,jpj,jpk,ntime))
     allocate(bfmvar2d(jpi,jpj,ntime))

     ! loop over the variables
     do d=1,n_bfmvar
        status=nf90_inquire_variable(IDncin, bfmvarid(d), ndims=ndims, name=varname)
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        if (old_version) then
           status = scan(varname,".")
           if (status > 0) varname = varname(1:2)//varname(4:4)
        end if

#ifdef DEBUG
        write(*,*) "Writing variable: "//trim(varname)
#endif
        status=nf90_inquire_variable(IDncin, bfmvarid(d), dimids=dimids(1:ndims))
        status = nf90_inquire_dimension(IDncin, dimids(1), len = dimlen(1))
        allocate(chunk(dimlen(1),ntime))
        ! read all time stamp of chunk data                                                    
        status = nf90_get_var(IDncin, bfmvarid(d), chunk(:,:), start = (/ 1, 1 /),     &            
                            count = (/ dimlen(1), ntime /))                                 
        if (status /= NF90_NOERR) call handle_err(status,errstring="variable: "//trim(varname))
        ! inquire ID in the output file
        status = nf90_inq_varid(IDncout, varname, IDvar)

        if (dimlen(1) == lenoce) then ! 3D variable
           bfmvar3d=NF90_FILL_REAL
           ! loop sequence is mandatory
           noce = 1
           do k = 1,jpk
              do j=1,jpj
                 do i=1,jpi
                    if (mask(i,j,k) > 0.0_RLEN) then
                       bfmvar3d(i,j,k,:) = chunk(noce,:)
                       noce = noce+1
                    end if
                 end do
              end do
           end do
           status = nf90_put_var(IDncout, IDvar, bfmvar3d, &
                    start = (/ 1, 1, 1, 1 /), count = (/ jpi, jpj, jpk, ntime /))
           if (status /= NF90_NOERR) &
              call handle_err(status,errstring="variable: "//trim(varname))
        else
           bfmvar2d=NF90_FILL_REAL
           ! loop sequence is mandatory
           noce = 1
           do j=1,jpj
              do i=1,jpi
                 if (mask(i,j,k) > 0.0_RLEN) then
                    bfmvar2d(i,j,:) = chunk(noce,:)
                    noce = noce+1
                 end if
              end do
           end do
           status = nf90_put_var(IDncout, IDvar, bfmvar2d, &
                    start = (/ 1, 1, 1 /), count = (/ jpi, jpj, ntime /))
           if (status /= NF90_NOERR) &
              call handle_err(status,errstring="variable: "//trim(varname))
        end if
        deallocate(chunk)
     end do ! variables

     ! close the netcdf file
     status = nf90_close(IDncin)
     if (status /= NF90_NOERR) call handle_err(status)
     deallocate(bfmvar3d)
     deallocate(bfmvar2d)

  ! write global grid specifications
     status = nf90_inq_varid(IDncout, "mask", IDvar)
     status = nf90_put_var(IDncout, IDvar, mask, start = (/ 1, 1, 1 /),     &
                            count = (/ jpi, jpj, jpk /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     status = nf90_inq_varid(IDncout, "lat", IDvar)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring variable: lat")
     status = nf90_put_var(IDncout, IDvar, lat, start = (/ 1, 1 /),     &
                            count = (/ jpi, jpj /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lat")
     status = nf90_inq_varid(IDncout, "lon", IDvar)
     status = nf90_put_var(IDncout, IDvar, lon, start = (/ 1, 1 /),     &
                            count = (/ jpi, jpj /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lon")

     status = nf90_close(IDncout)
     if (status /= NF90_NOERR) call handle_err(status)

  ! clean-up memory
  deallocate(mask)
  deallocate(lon)
  deallocate(lat)

end subroutine remap_vars
