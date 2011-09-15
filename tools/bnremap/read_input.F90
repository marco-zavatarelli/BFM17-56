#define DEBUG
! BFM_NEMO-MERGE bnremap
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

subroutine read_input

  use mod_bnremap
  implicit none
  integer,parameter  :: namlst=10
  integer            :: status,n,ndims,IDvar
  integer :: dimids(4),dimlen(4)
  real,allocatable,dimension(:,:,:,:) :: tmpmask

  namelist /bnremap_nml/ in_fname,out_fname, &
                         meshmask,mesh_flag,old_version


  ! Reading directory names and file name specification
  open(namlst,file='bnremap.nml',action='read',status='old',err=99)
  read(namlst,nml=bnremap_nml,err=98)
  close(namlst)

  ! open the input file
  status = nf90_open(path = trim(in_fname), mode = NF90_SHARE, ncid = IDncin)
  if (status /= NF90_NOERR) call handle_err(status)

  ! if requested, read mask from meshmask file
  if ( mesh_flag) then
     status = nf90_open(path = meshmask, mode = NF90_NOWRITE, ncid = IDmaskfile)
     if (status /= NF90_NOERR) call handle_err(status)
#ifdef DEBUG
    write(*,*) "Opening meshmask file:",trim(meshmask)
#endif
     maskname = "tmask"
     lonname = "nav_lon"
     latname = "nav_lat"
  else
     IDmaskfile = IDncin
     maskname = "mask"
     lonname  = "lon"
     latname  = "lat"
  end if

  if (old_version) then
     ocepointname = "ncomp"
  else
     ocepointname = "oceanpoint"
  end if

     ! read mask and lat-lon data
     status = nf90_inq_varid(IDmaskfile, trim(maskname), IDvar)
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     status=nf90_inquire_variable(IDmaskfile, IDvar, ndims=ndims)
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     status=nf90_inquire_variable(IDmaskfile, IDvar, dimids=dimids(1:ndims))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     do n=1,ndims
        status = nf90_inquire_dimension(IDmaskfile, dimids(n), len = dimlen(n))
     end do
     jpi = dimlen(1)
     jpj = dimlen(2)
     jpk = dimlen(3)
#ifdef DEBUG
     write(*,*) "dimensions from input file:",jpi,jpj,jpk
#endif
     allocate(mask(jpi,jpj,jpk))
     if (mesh_flag) then
        allocate(tmpmask(jpi,jpj,jpk,1))
        status = nf90_get_var(IDmaskfile, IDvar, tmpmask, start = (/ 1, 1, 1, 1 /),     &
                            count = (/ jpi, jpj, jpk, 1 /))
        if (status /= NF90_NOERR) call handle_err(status,errstring="reading variable: mask")
        mask(:,:,:) = tmpmask(:,:,:,1)
        deallocate(tmpmask)
     else
        status = nf90_get_var(IDmaskfile, IDvar, mask, start = (/ 1, 1, 1 /),     &
                            count = (/ jpi, jpj, jpk /))
     end if
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     write(*,*) "mask OK"
     allocate(lon(jpi,jpj))
     status = nf90_inq_varid(IDmaskfile, trim(lonname), IDvar)
     status = nf90_get_var(IDmaskfile, IDvar, lon, start = (/ 1, 1 /),     &
                            count = (/ jpi, jpj /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lon")
     write(*,*) "lon OK"
     allocate(lat(jpi,jpj))
     status = nf90_inq_varid(IDmaskfile, trim(latname), IDvar)
     status = nf90_get_var(IDmaskfile, IDvar, lat, start = (/ 1, 1 /),     &
                            count = (/ jpi, jpj /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lat")
     write(*,*) "lat OK"

  if ( mesh_flag) then
     status = nf90_close(IDmaskfile)
     if (status /= NF90_NOERR) call handle_err(status,errstring="closing meshmask file")
  end if

      return

99    stop 'I could not open bnremap.nml'
98    stop 'I could not read bnremap_nml'

end subroutine read_input
