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
module merge_vars
  use netcdf
  use mod_bnmerge, ONLY : handle_err, RLEN, ZERO, &
       jpiglo, jpjglo, jpkglo, latglo, longlo, maskglo

  implicit none
  real, allocatable, dimension(:) :: depth

  ! !PUBLIC MEMBER FUNCTIONS:
  public merge_vars_init

contains

  subroutine merge_vars_init
    use mod_bnmerge, ONLY : n_bfmvar_out, bfmvarid_out, bfmvartype_out, bfmvartarget_out, &
         n_bfmvar_res, bfmvarid_res, bfmvartype_res, bfmvartarget_res, &
         dimslen_out, dimslen_res, &
         do_output, do_restart, &
         out_dir, chunk_fname, bfm_restart, &
         jpnij

    use mod_bnmerge, ONLY : nimppt, njmppt, ibonit, ibonjt, nlcit , nlcjt, nldit, nldjt, nleit, nlejt


    implicit none
    integer            :: procnum
    integer            :: ncid_out, ncid_res
    integer            :: zflag,gflag

#ifdef DEBUG
    write(*,*) "Starting merging..."   
#endif

    ! intialize flags to not repeat operations
    zflag = 0 ! do only once
    gflag = 0 ! do only in output or restart

    ! Allocate global masks
    allocate(maskglo(jpiglo,jpjglo,jpkglo))
    allocate(latglo(jpiglo,jpjglo))
    allocate(longlo(jpiglo,jpjglo))
    allocate(depth(jpkglo))

    ! Initialisations
    maskglo=NF90_FILL_REAL
    latglo = NF90_FILL_REAL
    longlo = NF90_FILL_REAL
    depth = NF90_FILL_REAL

    ! open the output file
    if ( do_output ) then
       write(*,*) "Merging output init..."
       call handle_err(nf90_open(path = trim(out_dir)//"/"//trim(chunk_fname)//".nc", &
            mode = NF90_WRITE, ncid = ncid_out), &
            errstring="opening output file ")
    endif
    if( do_restart ) then
       write(*,*) "Merging restart init..."
       call handle_err(nf90_open(path = trim(out_dir)//"/"//trim(bfm_restart)//".nc", &
            mode = NF90_WRITE, ncid = ncid_res), &
            errstring="opening restart file ")
    end if

    !main execution
    if ( do_output ) then
       write(*,*) "Merging output main..."
       do procnum=1,jpnij
          call merge_vars_proc(.FALSE., chunk_fname, procnum, ncid_out, &
               n_bfmvar_out, bfmvarid_out, bfmvartype_out, bfmvartarget_out, dimslen_out, &
               zflag, gflag)
       end do ! processes
       gflag = 1
    end if
    if ( do_restart ) then
       write(*,*) "Merging restart main..."
       do procnum=1,jpnij
          call merge_vars_proc(.TRUE. , bfm_restart, procnum, ncid_res, &
               n_bfmvar_res, bfmvarid_res, bfmvartype_res, bfmvartarget_res, dimslen_res, &
               zflag, gflag)
       end do ! processes
    end if


    if ( do_output ) then
       call merge_vars_globals(ncid_out)
       write(*,*) "Finished merge output"
    endif
    if ( do_restart ) then
       call merge_vars_globals(ncid_res)
       write(*,*) "Finished merge restart"
    end if

    ! clean-up memory
    deallocate(maskglo)
    deallocate(longlo)
    deallocate(latglo)
    deallocate(depth)
    if( allocated(bfmvarid_out) ) deallocate(bfmvarid_out)
    if( allocated(bfmvarid_res) ) deallocate(bfmvarid_res)
    if( allocated(bfmvartype_out) ) deallocate(bfmvartype_out)
    if( allocated(bfmvartype_res) ) deallocate(bfmvartype_res)
    if( allocated(bfmvartarget_out) ) deallocate(bfmvartarget_out)
    if( allocated(bfmvartarget_res) ) deallocate(bfmvartarget_res)

    deallocate (nimppt)
    deallocate (njmppt)
    deallocate (ibonit)
    deallocate (ibonjt)
    deallocate (nlcit)
    deallocate (nlcjt)
    deallocate (nldit)
    deallocate (nldjt)
    deallocate (nleit)
    deallocate (nlejt)

  end subroutine merge_vars_init



  subroutine merge_vars_proc(is_restart, fname_in, procnum, ncid, n_vars, bfmvarid, bfmvartypes, bfmvartargets, dimslen, &
       zflag, gflag)
    use mod_bnmerge, ONLY : inp_dir, nimppt, njmppt, nlcit , nlcjt, TYPE_OCE, TYPE_BTN, TYPE_SRF, TYPE_PH

    implicit none
    integer, intent(inout)           :: zflag
    integer, intent(in)              :: procnum, ncid, n_vars, gflag
    integer,dimension(:), intent(in) :: bfmvarid, bfmvartypes, bfmvartargets, dimslen
    character(LEN=100), intent(in)   :: fname_in
    logical, intent(in)              :: is_restart

    integer             :: jpi, jpj ! number of grid points along i and j (proc)
    integer             :: i,j,k,d,n
    integer             :: ncinid, status
    character(LEN=172)  :: fname
    integer             :: ndims, nVars, nGlobalAtts, IDunlimdim
    integer             :: IDbtnpnt, IDvar

    integer           :: ntime
    integer           :: noce 
    integer           :: nimpp, njmpp

    real, allocatable, dimension(:,:)     :: lat, lon
    real, allocatable, dimension(:,:,:)   :: mask
    integer, allocatable, dimension(:)       :: bottompoint
    real(RLEN), allocatable, dimension(:,:)     :: chunk
    real(RLEN), allocatable, dimension(:,:,:,:) :: bfmvar3d
    real(RLEN), allocatable, dimension(:,:,:)   :: bfmvar2d

    integer :: vartype,dimids(4),dimlen(4)

    character(len = NF90_MAX_NAME) :: dimname
#ifdef DEBUG
    character(len = NF90_MAX_NAME) :: varname
#endif

    integer :: iniI, iniJ, finI, finJ
    integer :: Istart, Icount, Jstart, Jcount
    integer :: step_start_arr2(2), step_count_arr2(2), &
         step_start_arr3(3), step_count_arr3(3), &
         step_start_arr4(4), step_count_arr4(4)

    real(RLEN), allocatable, dimension(:) :: chunktmp
    integer                               :: bottompointtmp

    integer                        :: ID3dvars, ID2dbenvars, ID2dicevars

    character(LEN=4) :: procname

    ! build the file name for each process (start from 0)
    write(procname,'(I4.4)') procnum-1

    ! build the file name for each process (start from 0)
    fname = trim(inp_dir)//"/"//trim(fname_in)//"_"//procname//".nc"
    status = nf90_open(path = fname, mode = NF90_SHARE, ncid = ncinid)
    if (status /= NF90_NOERR) call handle_err(status)
    status = nf90_inquire(ncinid, nDims, nVars, nGlobalAtts, IDunlimdim)
    if (status /= NF90_NOERR) call handle_err(status)
    status = nf90_inquire_dimension(ncinid, IDunlimdim, len = ntime)
#ifdef DEBUG
    write(*,*)
    write(*,*) "Processing file : "
    write(*,*) trim(fname)
    write(*,*) "===================="
    write(*,*) "Domain:",procnum-1
    write(*,*) "No of dimensions = ",nDims
    write(*,*) "No of variables  = ",nVars
    write(*,*) "No of time step  = ",ntime
#endif

    ! read bottompoint data
    allocate(bottompoint(dimslen(TYPE_BTN)))
    status = nf90_inq_varid(ncinid, "bottompoint", IDbtnpnt)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring var bottompoint")
    status = nf90_get_var(ncinid, IDbtnpnt, bottompoint, start = (/ 1 /),     &
         count = (/ dimslen(TYPE_BTN) /))
    if (status /= NF90_NOERR) call handle_err(status,errstring="variable: bottompoint")

    ! read mask and lat-lon data
    status = nf90_inq_varid(ncinid, "mask", IDvar)
    if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
    status=nf90_inquire_variable(ncinid, IDvar, xtype=vartype, ndims=ndims)
    if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
    status=nf90_inquire_variable(ncinid, IDvar, dimids=dimids(1:ndims))
    if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")
    do n=1,ndims
       status = nf90_inquire_dimension(ncinid, dimids(n), len = dimlen(n))
    end do
    allocate(mask(dimlen(1),dimlen(2),dimlen(3)))
    status = nf90_get_var(ncinid, IDvar, mask, start = (/ 1, 1, 1 /),     &
         count = (/ dimlen(1),dimlen(2),dimlen(3) /))
    if (status /= NF90_NOERR) call handle_err(status,errstring="variable: mask")

    ! get the coordinates of the sub-domain
    nimpp = nimppt(procnum)
    njmpp = njmppt(procnum)
    jpi  = nlcit(procnum)
    jpj  = nlcjt(procnum)

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
    if ( gflag == 0 ) then
       allocate(lon(dimlen(1),dimlen(2)))
       status = nf90_inq_varid(ncinid, "lon", IDvar)
       status = nf90_get_var(ncinid, IDvar, lon, start = (/ 1, 1 /),     &
            count = (/ dimlen(1),dimlen(2) /))
       if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lon")
       allocate(lat(dimlen(1),dimlen(2)))
       status = nf90_inq_varid(ncinid, "lat", IDvar)
       status = nf90_get_var(ncinid, IDvar, lat, start = (/ 1, 1 /),     &
            count = (/ dimlen(1),dimlen(2) /))
       if (status /= NF90_NOERR) call handle_err(status,errstring="variable: lat")

       maskglo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1,:) = mask(iniI:jpi-finI,iniJ:jpj-finJ,:) 
       latglo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1)    = lat(iniI:jpi-finI,iniJ:jpj-finJ)
       longlo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1)    = lon(iniI:jpi-finI,iniJ:jpj-finJ)

       ! read vertical depths only once
       if (zflag == 0) then 
          call handle_err(nf90_inq_varid(ncinid, "z", IDvar),errstring="Error inquiring depth values")
          call handle_err(nf90_get_var(ncinid, IDvar,depth),errstring="Error in getting depth values")
          zflag = 1
       endif

#ifdef DEBUG
       write(*,*) "Mask type: ", vartype,NF90_FLOAT
       write(*,*) "size mask: ",size(mask,1),size(mask,2),size(mask,3)
       write(*,*) "size lat:  ",size(lat,1),size(lat,2)
       write(*,*) "size lon:  ",size(lon,1),size(lon,2)
       write(*,*) "domain specifications nimpp,jpi,njmpp,jpj,jpi,jpj,jpkglo,ntime:", &
            nimpp,jpi,njmpp,jpj,jpi,jpj,jpkglo,ntime
       write(*,*) "size maskglo: ",size(maskglo,1),size(maskglo,2),size(maskglo,3)
       write(*,*) "size latglo:  ",size(latglo,1),size(latglo,2)
       write(*,*) "size depth:   ",size(depth)
#endif

       deallocate(lon)
       deallocate(lat)
    endif

    ! allocate local 3D and 2D variables
    allocate(bfmvar3d(jpi,jpj,jpkglo,ntime))
    allocate(bfmvar2d(jpi,jpj,ntime))

    if ( is_restart ) then
       ! read variable IDs
       status = nf90_inq_varid(ncinid, "D3STATE", ID3dvars)
       status = nf90_inq_varid(ncinid, "D2STATE_BEN", ID2dbenvars)
       status = nf90_inq_varid(ncinid, "D2STATE_ICE", ID2dicevars)
    end if
    
    allocate(chunk(MAX(dimslen(TYPE_OCE), dimslen(TYPE_BTN), dimslen(TYPE_SRF)),ntime))

    ! loop over the variables
    do d=1,n_vars
       select case(bfmvartypes(d))
       case(TYPE_OCE)
          dimname="oceanpoint"
          if( is_restart ) then
             IDvar = ID3dvars
             step_start_arr2 = (/ bfmvarid(d), 1                 /)
             step_count_arr2 = (/ 1          , dimslen(TYPE_OCE) /)
          else           
             IDvar = bfmvarid(d)
             step_start_arr2 = (/ 1                , 1     /)
             step_count_arr2 = (/ dimslen(TYPE_OCE), ntime /)
          endif
       case(TYPE_BTN)
          dimname="bottompoint"
          if( is_restart ) then
             IDvar = ID2dbenvars
             step_start_arr2 = (/ bfmvarid(d), 1                 /)
             step_count_arr2 = (/ 1          , dimslen(TYPE_BTN) /)
          else           
             IDvar = bfmvarid(d)
             step_start_arr2 = (/ 1                , 1     /)
             step_count_arr2 = (/ dimslen(TYPE_BTN), ntime /)
          endif
       case(TYPE_SRF)
          dimname="surfacepoint"
          if( is_restart ) then
             IDvar = ID2dicevars
             step_start_arr2 = (/ bfmvarid(d), 1                 /)
             step_count_arr2 = (/ 1          , dimslen(TYPE_SRF) /)
          else           
             IDvar = bfmvarid(d)
             step_start_arr2 = (/ 1                , 1     /)
             step_count_arr2 = (/ dimslen(TYPE_SRF), ntime /)
          endif
       case(TYPE_PH)
          dimname="oceanpoint"
          IDvar = bfmvarid(d)
          step_start_arr2 = (/ 1                , 1     /)
          step_count_arr2 = (/ dimslen(TYPE_OCE), ntime /)
       end select
       ! read all time stamp of chunk data
       status = nf90_get_var(ncinid, IDvar, chunk, start = step_start_arr2, count = step_count_arr2 )
       if (status /= NF90_NOERR) call handle_err(status,errstring="Get variable in "//fname_in)
#ifdef DEBUG
       call handle_err(nf90_inquire_variable(ncid, bfmvartargets(d), name=varname),errstring="Getting name")
       write(*,'(A,i3)') "Writing variable : "//trim(varname)//", Proc: "//procname//", ID: ",bfmvartargets(d)
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
          bfmvar3d=NF90_FILL_DOUBLE
          ! loop sequence is mandatory
          noce = 0
          do k = 1,jpkglo
             do j=1,jpj
                do i=1,jpi
                   if (mask(i,j,k) > ZERO) then
                      noce = noce+1
                      bfmvar3d(i,j,k,:) = chunk(noce,:)
                   end if
                end do
             end do
          end do
          step_start_arr4 = (/ Istart, Jstart, 1, 1 /)
          step_count_arr4 = (/ Icount, Jcount, jpkglo, ntime /)
          status = nf90_put_var(ncid, bfmvartargets(d), bfmvar3d(iniI:jpi-finI,iniJ:jpj-finJ,:,:), &
               start = step_start_arr4, count = step_count_arr4)
          if (status /= NF90_NOERR) call handle_err(status,errstring="Put 3D domain: "//procname)
#ifdef DEBUG
          write(*,'(A,i3,i3,i3,i6)') "3D var "//trim(varname)//", Proc: "//procname//", Dims: ",jpi, jpj, jpkglo, ntime
#endif
       case ("surfacepoint")  ! 2D variable
          bfmvar2d=NF90_FILL_DOUBLE
          ! loop sequence is mandatory
          noce = 0
          do j=1,jpj
             do i=1,jpi
                if (mask(i,j,1) > ZERO) then
                   noce = noce+1
                   bfmvar2d(i,j,:) = chunk(noce,:)
                end if
             end do
          end do
          step_start_arr3 = (/ Istart, Jstart, 1 /)
          step_count_arr3= (/ Icount, Jcount, ntime /)
          status = nf90_put_var(ncid, bfmvartargets(d), bfmvar2d(iniI:jpi-finI,iniJ:jpj-finJ,:), &
               start = step_start_arr3, count = step_count_arr3)
          if (status /= NF90_NOERR) call handle_err(status,errstring="Put 2D domain: "//procname)
#ifdef DEBUG
          write(*,'(A,i3,i3,i3,i6)') "2D Surface var "//varname//", Proc: "//procname//", Dims: ",jpi, jpj, ntime
#endif
       case ("bottompoint")   ! 2D variable from the bottom

          ! allocate temporal array for bottompoint
          allocate(chunktmp(ntime))

          ! reorder bottom data
          do i=1 , dimslen(TYPE_BTN)-1
             do j=i+1 , dimslen(TYPE_BTN)
                if( bottompoint(i) .gt. bottompoint(j)  ) then
                   chunktmp(:)   = chunk(i,:)
                   chunk(i,:)    = chunk(j,:)
                   chunk(j,:)    = chunktmp(:)
                   bottompointtmp = bottompoint(i)
                   bottompoint(i) = bottompoint(j)
                   bottompoint(j) = bottompointtmp
                endif
             enddo
          enddo

          deallocate(chunktmp)

          bfmvar2d=NF90_FILL_DOUBLE
          ! loop sequence is mandatory
          noce = 0
          do j=1,jpj
             do i=1,jpi
                if (mask(i,j,1) > ZERO) then
                   noce = noce+1
                   bfmvar2d(i,j,:) = chunk(noce,:)
                end if
             end do
          end do
          step_start_arr3 = (/ Istart, Jstart, 1 /)
          step_count_arr3= (/ Icount, Jcount, ntime /)
          status = nf90_put_var(ncid, bfmvartargets(d), bfmvar2d(iniI:jpi-finI,iniJ:jpj-finJ,:), &
               start = step_start_arr3, count = step_count_arr3)
          if (status /= NF90_NOERR) call handle_err(status,errstring="Put 2D domain: "//procname)
#ifdef DEBUG
          write(*,'(A,i3,i3,i3,i6)') "2D Bottom var "//varname//", Proc: "//procname//", Dims: ",jpi, jpj, ntime
#endif
       end select

#ifdef DEBUG
       write(*,'(A,i6,A,i3)') "Wrote: ",noce," points for variable: "//trim(varname)//" ID: ", bfmvartargets(d)
#endif

    end do ! variables

    ! close the netcdf file
    call handle_err(nf90_close(ncinid))
    deallocate(mask)
    deallocate(bottompoint)
    deallocate(bfmvar3d)
    deallocate(bfmvar2d)
    deallocate(chunk)

#ifdef DEBUG
    status = nf90_inquire(ncid, nDims, nVars, nGlobalAtts, IDunlimdim)
    if (status /= NF90_NOERR) call handle_err(status)
    status = nf90_inquire_dimension(ncid, IDunlimdim, len = ntime)
    write(*,*) "Current number of time frames:",ntime
#endif
  end subroutine merge_vars_proc




  subroutine merge_vars_globals(ncid)
    use mod_bnmerge, ONLY : ln_mask
    
    implicit none
    integer, intent(in) :: ncid

    integer             :: IDvar
    integer             :: status
    
    ! write global grid specifications
    ! Oce-land points mask
    if (ln_mask) then
       status = nf90_inq_varid(ncid, "mask", IDvar)
       status = nf90_put_var(ncid, IDvar, real(maskglo,4), start = (/ 1, 1, 1 /),     &
            count = (/ jpiglo, jpjglo, jpkglo /))
       if (status /= NF90_NOERR) call handle_err(status,errstring="Output writing: mask")
    endif
    ! Latitude
    status = nf90_inq_varid(ncid, "lat", IDvar)
    if (status /= NF90_NOERR) call handle_err(status,errstring="inquiring variable in output: lat")
    status = nf90_put_var(ncid, IDvar, real(latglo,4), start = (/ 1, 1 /),     &
         count = (/ jpiglo, jpjglo /))
    if (status /= NF90_NOERR) call handle_err(status,errstring="Output writing: lat")
    ! Longitude
    status = nf90_inq_varid(ncid, "lon", IDvar)
    status = nf90_put_var(ncid, IDvar, real(longlo,4), start = (/ 1, 1 /),     &
         count = (/ jpiglo, jpjglo /))
    if (status /= NF90_NOERR) call handle_err(status,errstring="Output writing: lon")
    ! Depth levels
    call handle_err(nf90_inq_varid(ncid, "depth", IDvar))
    call handle_err(nf90_put_var(ncid, IDvar, real(depth,4)),errstring="Output writing: depth")   
    ! close file
    status = nf90_close(ncid)
    if (status /= NF90_NOERR) call handle_err(status,errstring="Closing otuput")
  end subroutine merge_vars_globals

end module merge_vars
