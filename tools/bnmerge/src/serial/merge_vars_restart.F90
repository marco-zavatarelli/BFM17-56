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
subroutine merge_vars_restart

  use netcdf
  use mod_bnmerge
  implicit none

  integer           :: ncresid, status
  integer           :: p,i,j,k,d,n
  character(LEN=4), allocatable, dimension(:)  :: procname
  character(LEN=172) :: fname
  integer           :: ndims, nVars, nGlobalAtts, IDmask
  integer           :: IDx, IDy, IDz, IDocepnt, IDsrfpnt, IDbtnpnt, IDvar, ncid
  integer           :: lenoce, lensrf, lenbtn
  integer           :: noce, nsrf 
  integer           :: nimpp, njmpp
  real, allocatable, dimension(:,:) :: lat, lon
  real, allocatable, dimension(:,:,:) :: mask
  real, allocatable, dimension(:) :: bottompoint, depth
  real, allocatable, dimension(:) :: chunk
  real, allocatable, dimension(:,:,:,:) :: bfmvar3d
  real, allocatable, dimension(:,:,:)   :: bfmvar2d
  integer :: vartype,dimids(4),dimlen(4),zflag
  character(len = NF90_MAX_NAME) :: dimname,varname,attname

  integer :: iniI, iniJ, finI, finJ, resI, resJ
  integer :: Istart, Icount, Jstart, Jcount
  integer :: step_start_arr3(3), step_count_arr3(3), step_start_arr4(4), step_count_arr4(4)

  real :: chunktmp
  real :: bottompointtmp

  integer                        :: ID3dvars, ID3dname, IDph
  character(len=64), allocatable :: res_names(:)

#ifdef DEBUG
  write(*,*) "Starting merging..."   
#endif

  zflag = 0

  ! Allocate global masks
  allocate(maskglo(jpiglo,jpjglo,jpk))
  allocate(latglo(jpiglo,jpjglo))
  allocate(longlo(jpiglo,jpjglo))
  allocate(depth(jpk))

  ! Initialisations
  maskglo=NF90_FILL_REAL
  latglo = NF90_FILL_REAL
  longlo = NF90_FILL_REAL
  depth = NF90_FILL_REAL
  ocepoints = 0
  srfpoints = 0

  allocate(procname(jpnij))
  do p=1,jpnij
     ! build the file name for each process (start from 0)
     write(procname(p),'(I4.4)') p-1
  end do

  ! open the output file
  fname = trim(out_dir)//"/"//trim(bfm_restart)//".nc"
  status = nf90_open(path = fname, mode = NF90_WRITE, ncid = ncid)
  if (status /= NF90_NOERR) call handle_err(status,errstring="opening named "//trim(fname)//".nc!" )
#ifdef DEBUG
  write(*,*) "Output file: "   
  write(*,*) trim(fname)
#endif

  do p=1,jpnij
     ! build the file name for each process (start from 0)
     fname = trim(inp_dir)//"/"//trim(bfm_restart)//"_"//procname(p)//".nc"
     status = nf90_open(path = fname, mode = NF90_SHARE, ncid = ncresid)
     if (status /= NF90_NOERR) call handle_err(status,errstring="opening named "//trim(fname)//".nc!" )
     status = nf90_inquire(ncresid, nDims, nVars, nGlobalAtts)
     if (status /= NF90_NOERR) call handle_err(status)
#ifdef DEBUG
     write(*,*)
     write(*,*) "Processing file : "
     write(*,*) trim(fname)
     write(*,*) "===================="
     write(*,*) "Domain:",p-1
     write(*,*) "No of dimensions = ",nDims
     write(*,*) "No of variables  = ",nVars
#endif
     ! inquire dimensions
     status = nf90_inq_dimid(ncresid, "x", IDx)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring x")
     status = nf90_inq_dimid(ncresid, "y", IDy)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring y")
     status = nf90_inq_dimid(ncresid, "z", IDz)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring z")
     status = nf90_inq_dimid(ncresid, "oceanpoint", IDocepnt)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring oceanpoint")
     status = nf90_inquire_dimension(ncresid, IDocepnt, len = lenoce)
     ! check if we have erroneously included land domains
     ! and skip the domain
     if (lenoce==1) cycle
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring dim oceanpoint")
     status = nf90_inq_dimid(ncresid, "surfacepoint", IDsrfpnt)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring surfacepoint")
     status = nf90_inquire_dimension(ncresid, IDsrfpnt, len = lensrf)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring dim surfacepoint")
     status = nf90_inq_dimid(ncresid, "bottompoint", IDbtnpnt)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring bottompoint")
     status = nf90_inquire_dimension(ncresid, IDbtnpnt, len = lenbtn)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring dim bottompoint")

     allocate(bottompoint(lenbtn))
     status = nf90_inq_varid(ncresid, "bottompoint", IDbtnpnt)
     if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring var bottompoint")
     status = nf90_get_var(ncresid, IDbtnpnt, bottompoint, start = (/ 1 /),     &
          count = (/ lenbtn /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: bottompoint")

     ! read vertical depths only once
     if (zflag == 0) then 
        call handle_err(nf90_inq_varid(ncresid, "z", IDvar),errstring="Error inquiring depth values")
        call handle_err(nf90_get_var(ncresid, IDvar,depth),errstring="Error in getting depth values")
        zflag = 1
     endif

     ! read mask and lat-lon data
     status = nf90_inq_varid(ncresid, "mask", IDvar)
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     status=nf90_inquire_variable(ncresid, IDvar, xtype=vartype, ndims=ndims)
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     status=nf90_inquire_variable(ncresid, IDvar, dimids=dimids(1:ndims))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     do n=1,ndims
        status = nf90_inquire_dimension(ncresid, dimids(n), len = dimlen(n))
     end do
     allocate(mask(dimlen(1),dimlen(2),dimlen(3)))
     status = nf90_get_var(ncresid, IDvar, mask, start = (/ 1, 1, 1 /),     &
          count = (/ dimlen(1),dimlen(2),dimlen(3) /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
     allocate(lon(dimlen(1),dimlen(2)))
     status = nf90_inq_varid(ncresid, "lon", IDvar)
     status = nf90_get_var(ncresid, IDvar, lon, start = (/ 1, 1 /),     &
          count = (/ dimlen(1),dimlen(2) /))
     if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lon")
     allocate(lat(dimlen(1),dimlen(2)))
     status = nf90_inq_varid(ncresid, "lat", IDvar)
     status = nf90_get_var(ncresid, IDvar, lat, start = (/ 1, 1 /),     &
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
     jpi  = nlcit(p)
     jpj  = nlcjt(p)
     !!dont write frame in chunk which overlap the other neighbours chunks
     iniI = 2 ; iniJ = 2
     finI = 1 ; finJ = 1
     if ( nimpp == 1 ) iniI = 1
     if ( njmpp == 1 ) iniJ = 1
     if( (nimpp + jpi) > jpiglo ) finI = 0
     if( (njmpp + jpj) > jpjglo ) finJ = 0 
     Istart = nimpp+iniI-1      ; Jstart = njmpp+iniJ-1
     Icount = jpi-finI-(iniI-1) ; Jcount = jpj-finJ-(iniJ-1)

     ! insert the subdomain in the global domain
     ! note that mask is double!
     maskglo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1,:) = mask(iniI:jpi-finI,iniJ:jpj-finJ,:) 
     latglo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1)    = lat(iniI:jpi-finI,iniJ:jpj-finJ)
     longlo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1)    = lon(iniI:jpi-finI,iniJ:jpj-finJ)
     deallocate(lon)
     deallocate(lat)

     ! allocate local 3D and 2D variables
     allocate(bfmvar3d(jpi,jpj,jpk,1))
     allocate(bfmvar2d(jpi,jpj,1))
#ifdef DEBUG
     write(*,*) "domain specifications nimpp,jpi,njmpp,jpj,jpi,jpj,jpk:", &
          nimpp,jpi,njmpp,jpj,jpi,jpj,jpk
     write(*,*) "size maskglo:",size(maskglo,1),size(maskglo,2),size(maskglo,3)
     write(*,*) "size latglo:",size(latglo,1),size(latglo,2)
#endif

     ! read variable data
     status = nf90_inq_varid(ncresid, "D3STATE", ID3dvars)
     if (status == NF90_NOERR) then
        status=nf90_inquire_variable(ncresid, ID3dvars, ndims=ndims)
        if (status /= NF90_NOERR) call handle_err(status,errstring="Inquire variable: D3STATE")
        status=nf90_inquire_variable(ncresid, ID3dvars, dimids=dimids(1:ndims))
        status = nf90_inquire_dimension(ncresid, dimids(2), name=dimname, len = dimlen(2))
        allocate(chunk(dimlen(2)))
        allocate(res_names(n_bfmvar3d))
        status = nf90_inq_varid(ncresid, "D3STATE_NAME", ID3dname)
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_NAME var in "//fname)
        status = nf90_get_var(ncresid, ID3dname, res_names, start=(/ 1, 1 /), count=(/ LEN(res_names), n_bfmvar3d-1 /))
        if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring D3STATE_NAME var in "//fname)
        res_names(n_bfmvar3d) = "pH"

        WRITE(*,*) "DIMLEN: ", dimlen(2)
        ! loop over the variables
        do d=1,n_bfmvar3d
           varname = res_names(d)
           ! read all time stamp of chunk data                                                    
           if( varname .eq. "pH" ) then
              status = nf90_inq_varid(ncresid, "pH", IDph)
              if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring pH var in "//fname)
              status = nf90_get_var(ncresid, IDph, chunk(:), start = (/ 1 /), count = (/ dimlen(2) /))
              if (status /= NF90_NOERR) call handle_err(status,errstring="Get variable ID: "//trim(varname))
           else
              status = nf90_get_var(ncresid, ID3dvars, chunk(:), start = (/ d, 1 /), count = (/ 1, dimlen(2) /))
              if (status /= NF90_NOERR) call handle_err(status,errstring="Get variable ID: "//trim(varname))
           end if
           ! inquire ID in the output file
           status = nf90_inq_varid(ncid, varname, IDvar)
#ifdef DEBUG
           write(*,*)
           write(*,*) "Writing variable: "//trim(varname)
           write(*,*) "BFM Dimension:",dimlen(2)
#endif

           iniI = 2
           iniJ = 2
           finI = 1
           finJ = 1
           if ( nimpp == 1 ) then
              iniI = 1
           end if
           if ( njmpp == 1 ) then
              iniJ = 1
           end if
           if( (nimpp + jpi) > jpiglo ) then
              finI = 0
           endif
           if( (njmpp + jpj) > jpjglo ) then
              finJ = 0
           endif

           select case (dimname)

           case("oceanpoint")    ! 3D variable
              bfmvar3d=NF90_FILL_REAL
              ! loop sequence is mandatory
              noce = 1
              do k = 1,jpk
                 do j=1,jpj
                    do i=1,jpi
                       if (mask(i,j,k) > 0.0_RLEN) then
                          bfmvar3d(i,j,k,1) = chunk(noce)
                          noce = noce+1
                       end if
                    end do
                 end do
              end do
              step_start_arr4 = (/ Istart, Jstart, 1, 1 /)
              step_count_arr4 = (/ Icount, Jcount, jpk, 1 /)
              status = nf90_put_var(ncid, IDvar, bfmvar3d(iniI:jpi-finI,iniJ:jpj-finJ,:,1), &
                   start = step_start_arr4, count = step_count_arr4)
#ifdef DEBUG
              write(*,*) "3D var, Dimensions: ",jpi, jpj, jpk
#endif
              if (status /= NF90_NOERR) &
                   call handle_err(status,errstring="Put 3D domain: "//procname(p)//" variable: "//trim(varname))

           case ("surfacepoint")  ! 2D variable
              bfmvar2d=NF90_FILL_REAL
              ! loop sequence is mandatory
              noce = 1
              do j=1,jpj
                 do i=1,jpi
                    if (mask(i,j,1) > 0.0_RLEN) then
                       bfmvar2d(i,j,1) = chunk(noce)
                       noce = noce+1
                    end if
                 end do
              end do
              step_start_arr3 = (/ Istart, Jstart, 1 /)
              step_count_arr3= (/ Icount, Jcount, 1 /)
              status = nf90_put_var(ncid, IDvar, bfmvar2d(iniI:jpi-finI,iniJ:jpj-finJ,1), &
                   start = step_start_arr3, count = step_count_arr3)
#ifdef DEBUG
              write(*,*) "2D Surface var, Dimensions: ",jpi, jpj
#endif
              if (status /= NF90_NOERR) &
                   call handle_err(status,errstring="Put 2D domain: "//procname(p)//" variable: "//trim(varname))

           case ("bottompoint")   ! 2D variable from the bottom

              ! reorder bottom data
              do i=1 , lenbtn-1
                 do j=i+1 , lenbtn
                    if( bottompoint(i) .gt. bottompoint(j)  ) then
                       chunktmp    = chunk(i)
                       chunk(i)    = chunk(j)
                       chunk(j)    = chunktmp
                       bottompointtmp = bottompoint(i)
                       bottompoint(i) = bottompoint(j)
                       bottompoint(j) = bottompointtmp
                    endif
                 enddo
              enddo

              bfmvar2d=NF90_FILL_REAL
              ! loop sequence is mandatory
              noce = 1
              do j=1,jpj
                 do i=1,jpi
                    if (mask(i,j,1) > 0.0_RLEN) then
                       bfmvar2d(i,j,1) = chunk(noce)
                       noce = noce+1
                    end if
                 end do
              end do
              step_start_arr3 = (/ Istart, Jstart, 1 /)
              step_count_arr3= (/ Icount, Jcount, 1 /)
              status = nf90_put_var(ncid, IDvar, bfmvar2d(iniI:jpi-finI,iniJ:jpj-finJ, 1), &
                   start = step_start_arr3, count = step_count_arr3)
#ifdef DEBUG
              write(*,*) "2D Bottom var, Dimensions: ",jpi, jpj
#endif
              if (status /= NF90_NOERR) &
                   call handle_err(status,errstring="Put 2D domain: "//procname(p)//" variable: "//trim(varname))
           end select

        end do ! variables

        deallocate(chunk)
        deallocate(res_names)

     end if

     ! close the netcdf file
     call handle_err(nf90_close(ncresid))
     deallocate(mask)
     deallocate(bottompoint)
     deallocate(bfmvar3d)
     deallocate(bfmvar2d)

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

end subroutine merge_vars_restart
