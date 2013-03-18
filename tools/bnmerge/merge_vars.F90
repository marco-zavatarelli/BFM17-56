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
subroutine merge_vars
  use netcdf
  use mod_bnmerge
  implicit none
  integer           :: ncbfmid, status
  integer           :: p,i,j,k,d,n
  character(LEN=4), allocatable, dimension(:)  :: procname
  character(LEN=172) :: fname
  integer           :: ndims, nVars, nGlobalAtts, IDunlimdim, IDmask
  integer           :: IDx, IDy, IDz, IDocepnt, IDsrfpnt, IDvar, ncid
  integer           :: lenoce, lensrf, ntime
  integer           :: noce, nsrf 
  integer           :: nimpp, njmpp, nlci, nlcj
  real, allocatable, dimension(:,:) :: lat, lon
  real, allocatable, dimension(:,:,:) :: mask
  real, allocatable, dimension(:) :: oceanpoint, depth
  real, allocatable, dimension(:,:) :: chunk
  real, allocatable, dimension(:,:,:,:) :: bfmvar3d
  real, allocatable, dimension(:,:,:)   :: bfmvar2d
  integer :: vartype,dimids(4),dimlen(4),zflag
  character(len = NF90_MAX_NAME) :: DimName,varname,attname

  zflag = 0

  ! Allocate global masks
  allocate(maskglo(jpiglo,jpjglo,jpk))
  allocate(latglo(jpiglo,jpjglo))
  allocate(longlo(jpiglo,jpjglo))
  allocate(depth(jpk))

  ! Initialisations
  maskglo=NF90_FILL_REAL
  latglo = 0
  longlo = 0
  depth = 0
  ocepoints = 0
  srfpoints = 0

  allocate(procname(jpnij))
  do p=1,jpnij
     ! build the file name for each process (start from 0)
     write(procname(p),'(I4.4)') p-1
  end do

  ! open the output file
     fname = trim(out_dir)//"/"//trim(chunk_fname)//".nc"
     status = nf90_open(path = fname, mode = NF90_WRITE, ncid = ncid)
     if (status /= NF90_NOERR) call handle_err(status)
#ifdef DEBUG
     write(*,*) "Output file: "   
     write(*,*) trim(fname)
#endif

  do p=1,jpnij
     ! build the file name for each process (start from 0)
     fname = trim(inp_dir)//"/"//trim(chunk_fname)//"_"//procname(p)//".nc"
     status = nf90_open(path = fname, mode = NF90_SHARE, ncid = ncbfmid)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire(ncbfmid, nDims, nVars, nGlobalAtts, IDunlimdim)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(ncbfmid, IDunlimdim, len = ntime)
#ifdef DEBUG
     write(*,*)
     write(*,*) "Processing file : "
     write(*,*) trim(fname)
     write(*,*) "===================="
     write(*,*) "Domain:",p-1
     write(*,*) "No of dimensions = ",nDims
     write(*,*) "No of variables  = ",nVars
     write(*,*) "No of time step  = ",ntime
#endif
     ! inquire dimensions
     status = nf90_inq_dimid(ncbfmid, "x", IDx)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inq_dimid(ncbfmid, "y", IDy)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inq_dimid(ncbfmid, "z", IDz)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inq_dimid(ncbfmid, "oceanpoint", IDocepnt)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(ncbfmid, IDocepnt, len = lenoce)
     ! check if we have erroneously included land domains
     ! and skip the domain
     if (lenoce==1) cycle
     status = nf90_inq_dimid(ncbfmid, "surfacepoint", IDsrfpnt)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(ncbfmid, IDsrfpnt, len = lensrf)

     ! read oceanpoint data
     status = nf90_inq_varid(ncbfmid, "oceanpoint", IDvar)
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: oceanpoint")
     allocate(oceanpoint(lenoce))
     status = nf90_get_var(ncbfmid, IDvar, oceanpoint, start = (/ 1 /),     &
                            count = (/ lenoce /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: oceanpoint")

     ! read vertical depths only once
     if (zflag == 0) then 
     call handle_err(nf90_inq_varid(ncbfmid, "z", IDvar),errstring="Error inquiring depth values")
     call handle_err(nf90_get_var(ncbfmid, IDvar,depth),errstring="Error in getting depth values")
     zflag = 1
     endif
     ! read mask and lat-lon data
     status = nf90_inq_varid(ncbfmid, "mask", IDvar)
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     status=nf90_inquire_variable(ncbfmid, IDvar, xtype=vartype, ndims=ndims)
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     status=nf90_inquire_variable(ncbfmid, IDvar, dimids=dimids(1:ndims))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     do n=1,ndims
        status = nf90_inquire_dimension(ncbfmid, dimids(n), len = dimlen(n))
     end do
     allocate(mask(dimlen(1),dimlen(2),dimlen(3)))
     status = nf90_get_var(ncbfmid, IDvar, mask, start = (/ 1, 1, 1 /),     &
                            count = (/ dimlen(1),dimlen(2),dimlen(3) /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     allocate(lon(dimlen(1),dimlen(2)))
     status = nf90_inq_varid(ncbfmid, "lon", IDvar)
     status = nf90_get_var(ncbfmid, IDvar, lon, start = (/ 1, 1 /),     &
                            count = (/ dimlen(1),dimlen(2) /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lon")
     allocate(lat(dimlen(1),dimlen(2)))
     status = nf90_inq_varid(ncbfmid, "lat", IDvar)
     status = nf90_get_var(ncbfmid, IDvar, lat, start = (/ 1, 1 /),     &
                            count = (/ dimlen(1),dimlen(2) /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lat")
#ifdef DEBUG
     write(*,*) "Mask type", vartype,NF90_FLOAT
     write(*,*) "size mask:",size(mask,1),size(mask,2),size(mask,3)
     write(*,*) "size lat:",size(lat,1),size(lat,2)
     write(*,*) "size lon:",size(lon,1),size(lon,2)
#endif

     ! get the coordinates of the sub-domain
     nimpp = nimppt(p)
     njmpp = njmppt(p)
     nlci  = nlcit(p)
     nlcj  = nlcjt(p)
     ! insert the subdomain in the global domain
     ! note that mask is double!
     maskglo(nimpp:nimpp+nlci-1,njmpp:njmpp+nlcj-1,:) = mask(:,:,:)
     latglo(nimpp:nimpp+nlci-1,njmpp:njmpp+nlcj-1) = lat(:,:)
     longlo(nimpp:nimpp+nlci-1,njmpp:njmpp+nlcj-1) = lon(:,:)
     deallocate(lon)
     deallocate(lat)

     ! allocate local 3D and 2D variables
     jpi = nlci
     jpj = nlcj
     allocate(bfmvar3d(jpi,jpj,jpk,ntime))
     allocate(bfmvar2d(jpi,jpj,ntime))
#ifdef DEBUG
     write(*,*) "domain specifications nimpp,nlci,njmpp,nlcj,jpi,jpj,jpk,ntime:", &
                nimpp,nlci,njmpp,nlcj,jpi,jpj,jpk,ntime
     write(*,*) "size maskglo:",size(maskglo,1),size(maskglo,2),size(maskglo,3)
     write(*,*) "size latglo:",size(latglo,1),size(latglo,2)
#endif

     ! loop over the variables
     do d=1,n_bfmvar
        status=nf90_inquire_variable(ncbfmid, bfmvarid(d), xtype=vartype, ndims=ndims, name=varname)
        if (status /= NF90_NOERR) call handle_err(status,errstring="Inquire variable: "//trim(varname))
        status=nf90_inquire_variable(ncbfmid, bfmvarid(d), dimids=dimids(1:ndims))
        status = nf90_inquire_dimension(ncbfmid, dimids(1), len = dimlen(1))
        allocate(chunk(dimlen(1),ntime))
        ! read all time stamp of chunk data                                                    
        status = nf90_get_var(ncbfmid, bfmvarid(d), chunk(:,:), start = (/ 1, 1 /),     &            
                            count = (/ dimlen(1), ntime /))                                 
        if (status /= NF90_NOERR) call handle_err(status,errstring="Get variable: "//trim(varname))
        ! inquire ID in the output file
        status = nf90_inq_varid(ncid, varname, IDvar)
#ifdef DEBUG
        write(*,*)
        write(*,*) "Writing variable: "//trim(varname)
        write(*,*) "BFM Dimension:",dimlen(1)
#endif

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
           status = nf90_put_var(ncid, IDvar, bfmvar3d, &
                    start = (/ nimpp, njmpp, 1, 1 /), count = (/ jpi, jpj, jpk, ntime /))
#ifdef DEBUG
        write(*,*) "3D var, Dimensions: ",jpi, jpj, jpk, ntime
#endif
           if (status /= NF90_NOERR) &
              call handle_err(status,errstring="Put 3D domain: "//procname(p)//" variable: "//trim(varname))
        else
           bfmvar2d=NF90_FILL_REAL
           ! loop sequence is mandatory
           noce = 1
           do j=1,jpj
              do i=1,jpi
                 if (mask(i,j,1) > 0.0_RLEN) then
                    bfmvar2d(i,j,:) = chunk(noce,:)
                    noce = noce+1
                 end if
              end do
           end do
           status = nf90_put_var(ncid, IDvar, bfmvar2d, &
                    start = (/ nimpp, njmpp, 1 /), count = (/ jpi, jpj, ntime /))
#ifdef DEBUG
        write(*,*) "2D var, Dimensions: ",jpi, jpj, ntime
#endif
           if (status /= NF90_NOERR) &
              call handle_err(status,errstring="Put 2D domain: "//procname(p)//" variable: "//trim(varname))
        end if
        deallocate(chunk)
     end do ! variables

     ! close the netcdf file
     call handle_err(nf90_close(ncbfmid))
     deallocate(mask)
     deallocate(oceanpoint)
     deallocate(bfmvar3d)
     deallocate(bfmvar2d)

#ifdef DEBUG
     status = nf90_inquire(ncid, nDims, nVars, nGlobalAtts, IDunlimdim)
     if (status /= NF90_NOERR) call handle_err(status)
     status = nf90_inquire_dimension(ncid, IDunlimdim, len = ntime)
     write(*,*) "Current number of time frames:",ntime
#endif
  end do ! processes

  ! write global grid specifications 
  ! Oce-land points mask
  if (ln_mask) then
     status = nf90_inq_varid(ncid, "mask", IDmask)
     status = nf90_put_var(ncid, IDmask, real(maskglo,4), start = (/ 1, 1, 1 /),     &
                            count = (/ jpiglo, jpjglo, jpk /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="Writing: mask")
  endif
  ! Latitude
     status = nf90_inq_varid(ncid, "lat", IDmask)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring variable: lat")
     status = nf90_put_var(ncid, IDmask, real(latglo,4), start = (/ 1, 1 /),     &
                            count = (/ jpiglo, jpjglo /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="Writing: lat")
  ! Longitude
     status = nf90_inq_varid(ncid, "lon", IDmask)
     status = nf90_put_var(ncid, IDmask, real(longlo,4), start = (/ 1, 1 /),     &
                            count = (/ jpiglo, jpjglo /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="Writing: lon")
  ! Depth levels
  call handle_err(nf90_inq_varid(ncid, "depth", IDmask))
  call handle_err(nf90_put_var(ncid, IDmask, real(depth,4)),errstring="Writing: depth")   
  ! close file
     status = nf90_close(ncid)
     if (status /= NF90_NOERR) call handle_err(status)

  ! clean-up memory
  deallocate(maskglo)
  deallocate(longlo)
  deallocate(latglo)
  deallocate(depth)
end subroutine merge_vars
