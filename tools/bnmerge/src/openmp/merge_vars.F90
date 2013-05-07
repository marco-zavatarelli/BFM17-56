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
  use pnetcdf
  use mod_bnmerge
  use mpi

  !$ USE omp_lib           ! Note OpenMP sentinel

  implicit none
  character(LEN=172) :: fname
  integer ncid

  integer, parameter :: FID=1, NDIMS=FID+1, NVARS=NDIMS+1, NGATTS=NVARS+1
  integer, parameter :: IDUNDIM=NGATTS+1, IDMASK=IDUNDIM+1, IDOCE=IDMASK+1, IDLON=IDOCE+1, IDLAT=IDLON+1 
  integer, parameter :: NTIMEK=1, LENOCEK=NTIMEK+1, LENSRFK=LENOCEK+1
  integer, parameter :: NPARS=IDLAT, NPARSK=LENSRFK

  integer :: p, n, id, nthreads=0, nthread=0, d
  integer :: nimpp, njmpp
  real(4), allocatable, dimension(:,:) :: lat, lon
  integer(kind=MPI_OFFSET_KIND) :: step_start(1), step_count(1), &
       step_start_arr3(3), step_count_arr3(3), &
       step_start_arr2(2), step_count_arr2(2), &
       step_start_arr4(4), step_count_arr4(4)

  integer, allocatable, dimension(:,:)                       :: A_infofile
  integer(kind=MPI_OFFSET_KIND), allocatable, dimension(:,:) :: A_infofile_kind
  character(LEN=4), allocatable, dimension(:)                :: A_procname
  character(LEN=172), allocatable, dimension(:)              :: A_fname
  integer(kind=MPI_OFFSET_KIND), allocatable, dimension(:,:) :: A_dimlen
  integer, allocatable, dimension(:,:)                       :: A_dimids
  real(4), dimension(:), pointer                             :: depth !only read once the vertical depth

  real(4), allocatable, dimension(:,:,:,:) :: bfmvar3d
  real(4), allocatable, dimension(:,:,:)   :: bfmvar2d

  character(len = NF_MAX_NAME) :: varname
  integer :: vartype, varndims, vardimids(4), varnatts
  integer(kind=MPI_OFFSET_KIND) :: vardimlen(4)

  character(len = 200) :: tmp_str
  integer noce, vardID, i, j, k, vardIDmask, vardIDlat, vardIDlon, vardIDdep

  integer iniI, iniJ, finI, finJ, resI, resJ
  integer(kind=MPI_OFFSET_KIND) :: Istart, Icount, Jstart, Jcount

  type ppChunk
     real(4), allocatable, dimension(:,:) :: p
  end type ppChunk
  type(ppChunk), dimension(:,:), allocatable :: A_chunk
  type(ppChunk), dimension(:),   allocatable :: A_lat, A_lon
  type ppOce
     real(4), dimension(:), allocatable :: p
  end type ppOce
  type(ppOce), dimension(:), allocatable :: A_oceanpoint
  type ppMask
     real(4), dimension(:,:,:), allocatable :: p
  end type ppMask
  type(ppMask), dimension(:), allocatable :: A_mask


  !--------------------------------------------------------------------------------
  !0. OMP setup

  !$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
  !$OMP MASTER
  !$      nthreads = omp_get_num_threads()
  !$      WRITE(*,*) 'Running OMP with ',nthreads,' thread(s).'
  !$      CALL OMP_SET_NESTED(.FALSE.)      ! disables nested parallelism
  !$OMP END MASTER
  !$OMP END PARALLEL

  !--------------------------------------------------------------------------------

  ! Allocate global masks
  allocate(maskglo(jpiglo,jpjglo,jpk))
  allocate(latglo(jpiglo,jpjglo))
  allocate(longlo(jpiglo,jpjglo))
  allocate(depth(jpk))

  ! Initialisations
  maskglo = NF_FILL_REAL
  latglo  = NF_FILL_REAL
  longlo  = NF_FILL_REAL
  depth   = NF_FILL_REAL


  ! Allocate array of pointers
  allocate(A_oceanpoint(jpnij))
  allocate(A_mask(jpnij))
  allocate(A_lat(jpnij))
  allocate(A_lon(jpnij))
  allocate(A_chunk(jpnij,n_bfmvar))

  ! Allocate arrays
  allocate(A_infofile(jpnij,NPARS))
  allocate(A_infofile_kind(jpnij,NPARSK))
  allocate(A_procname(jpnij))
  allocate(A_fname(jpnij))
  allocate(A_dimlen(jpnij,4))
  allocate(A_dimids(jpnij,4))



  !open the output file
  fname = trim(out_dir)//"/"//trim(chunk_fname)//".nc"
  call handle_err( nfmpi_open(MPI_COMM_WORLD, path = fname, mode = NF_WRITE, mpi_info=MPI_INFO_NULL, ncid = ncid) )
#ifdef DEBUG
  write(*,*) "Output file: ", trim(fname)
#endif 

  !open the first chunk to read vertical depth once
  write(A_procname(1),'(I4.4)') 0
  A_fname(1) = trim(inp_dir)//"/"//trim(chunk_fname)//"_"//A_procname(1)//".nc"
  call handle_err( nfmpi_open(MPI_COMM_WORLD, path = A_fname(1), mode = NF_NOWRITE, mpi_info=MPI_INFO_NULL, ncid = A_infofile(1,FID)))
  call handle_err( nfmpi_inq_varid(A_infofile(1,FID), "z", id),errstring="Error inquiring depth values")
  call handle_err( nfmpi_begin_indep_data( A_infofile(1,FID)) )
  call handle_err( nfmpi_get_var_real(A_infofile(1,FID), id, depth),errstring="Error in getting depth values")
  call handle_err( nfmpi_end_indep_data( A_infofile(1,FID)) )
  call handle_err( nfmpi_close(A_infofile(1,FID)) )

  do p=1,jpnij
     !open chunk file
     write(A_procname(p),'(I4.4)') p-1
     A_fname(p) = trim(inp_dir)//"/"//trim(chunk_fname)//"_"//A_procname(p)//".nc"
     call handle_err( nfmpi_open(MPI_COMM_WORLD, path = A_fname(p), mode = NF_NOWRITE, mpi_info=MPI_INFO_NULL, ncid = A_infofile(p,FID)))
  end do

  ! loop over all the files getting info
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP PRIVATE ( id, nthread, nimpp, njmpp, lat, lon )                      &
  !$OMP PRIVATE ( step_start, step_count, step_start_arr2, step_count_arr2, step_start_arr3, step_count_arr3, step_start_arr4, step_count_arr4 ) &
  !$OMP SHARED  ( inp_dir, chunk_fname, A_fname )     &
  !$OMP SHARED  ( jpnij, nimppt, njmppt, nlcit, nlcjt, maskglo, latglo, longlo )     &
  !$OMP SHARED  ( A_infofile, A_infofile_kind, A_oceanpoint, A_chunk, A_procname, A_mask, A_lon, A_lat, A_dimlen, A_dimids ) &
  !$OMP PRIVATE ( d, jpi,jpj, bfmvar2d, bfmvar3d, varname, vartype, varndims, vardimids, vardimlen, varnatts, tmp_str ) &
  !$OMP SHARED  ( jpiglo, jpjglo, jpk, n_bfmvar, bfmvarid, ncid ) &
  !$OMP PRIVATE ( vardID, noce, i, j, k ) &
  !$OMP PRIVATE ( iniI, iniJ, finI, finJ, resI, resJ, Istart, Jstart, Icount, Jcount )
  do p=1,jpnij
#ifdef DEBUG
     !$ nthread = omp_get_thread_num()
     !$ WRITE(*,*) 'Running OMP thread: ', nthread, ' ID: ', p
#endif

     ! get attributes
     call handle_err(nfmpi_inq(A_infofile(p,FID), A_infofile(p,NDIMS), A_infofile(p,NVARS), &
          A_infofile(p,NGATTS), A_infofile(p,IDUNDIM)))

     ! get dimensions
     call handle_err(nfmpi_inq_dimlen(A_infofile(p,FID), A_infofile(p,IDUNDIM), len = A_infofile_kind(p,NTIMEK)))
     call handle_err(nfmpi_inq_dimid(A_infofile(p,FID), "x", id))
     call handle_err(nfmpi_inq_dimid(A_infofile(p,FID), "y", id))
     call handle_err(nfmpi_inq_dimid(A_infofile(p,FID), "z", id))
     call handle_err(nfmpi_inq_dimid(A_infofile(p,FID), "oceanpoint", A_infofile(p,IDOCE)))
     call handle_err(nfmpi_inq_dimlen(A_infofile(p,FID), A_infofile(p,IDOCE), len = A_infofile_kind(p,LENOCEK)))

     ! check if we have erroneously included land domains and skip the domain
     if ( A_infofile_kind(p,LENOCEK)==1 ) call handle_err( NF_EBADID, errstring="land domains"  )
     call handle_err(nfmpi_inq_dimid(A_infofile(p,FID), "surfacepoint", id))
     call handle_err(nfmpi_inq_dimlen(A_infofile(p,FID), id, len = A_infofile_kind(p,LENSRFK)))

     ! read mask and lat-lon data
     call handle_err(nfmpi_inq_varid(A_infofile(p,FID), "mask", A_infofile(p,IDMASK)),errstring="variable: mask id")
     call handle_err(nfmpi_inq_varndims(A_infofile(p,FID), A_infofile(p,IDMASK), ndims=A_infofile(p,NDIMS)), &
          errstring="variable: mask dims")
     call handle_err(nfmpi_inq_vardimid(A_infofile(p,FID), A_infofile(p,IDMASK), dimids=A_dimids(p,1:A_infofile(p,NDIMS))), &
          errstring="variable: mask dims ID")
     do n=1,A_infofile(p,NDIMS)
        call handle_err(nfmpi_inq_dimlen(A_infofile(p,FID), A_dimids(p,n), len = A_dimlen(p,n)))
     end do
     call handle_err(nfmpi_inq_varid(A_infofile(p,FID), "lon", A_infofile(p,IDLON)))
     call handle_err(nfmpi_inq_varid(A_infofile(p,FID), "lat", A_infofile(p,IDLAT)))

     ! read oceanpoint
     allocate( A_oceanpoint(p)%p(A_infofile_kind(p,LENOCEK)) )
     step_start = (/ 1 /)
     step_count = (/ A_infofile_kind(p,LENOCEK) /)
     call handle_err(nfmpi_get_vara_real_all(A_infofile(p,FID), A_infofile(p,IDOCE), step_start, step_count, A_oceanpoint(p)%p ), &
          errstring="variable: oceanpoint in chunks "//trim(A_procname(p)))

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

     ! read mask
     allocate( A_mask(p)%p(A_dimlen(p,1),A_dimlen(p,2),A_dimlen(p,3)) )
     step_start_arr3 = (/ 1, 1, 1 /)
     step_count_arr3 = (/ A_dimlen(p,1),A_dimlen(p,2),A_dimlen(p,3) /)
     call handle_err(nfmpi_get_vara_real_all(A_infofile(p,FID), A_infofile(p,IDMASK), step_start_arr3, step_count_arr3, A_mask(p)%p), &
          errstring="variable: mask get")
     maskglo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1,:) = A_mask(p)%p(iniI:jpi-finI,iniJ:jpj-finJ,:)
     ! read lon
     allocate(A_lon(p)%p(A_dimlen(p,1),A_dimlen(p,2)))
     step_start_arr2 = (/ 1, 1 /)
     step_count_arr2 = (/ A_dimlen(p,1),A_dimlen(p,2) /)
     call handle_err(nfmpi_get_vara_real_all(A_infofile(p,FID), A_infofile(p,IDLON), step_start_arr2, step_count_arr2, A_lon(p)%p ), & 
          errstring="variable: lon")
     longlo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1) = A_lon(p)%p(iniI:jpi-finI,iniJ:jpj-finJ)

     ! read lat
     allocate(A_lat(p)%p(A_dimlen(p,1),A_dimlen(p,2)))
     step_start_arr2 = (/ 1, 1 /)
     step_count_arr2 = (/ A_dimlen(p,1), A_dimlen(p,2) /)
     call handle_err(nfmpi_get_vara_real_all(A_infofile(p,FID), A_infofile(p,IDLAT), step_start_arr2, step_count_arr2, A_lat(p)%p ), &
          errstring="variable: lat")
     latglo(Istart:Istart+Icount-1,Jstart:Jstart+Jcount-1) = A_lat(p)%p(iniI:jpi-finI,iniJ:jpj-finJ)

     allocate( bfmvar3d(jpi,jpj,jpk,A_infofile_kind(p,NTIMEK)) )
     allocate( bfmvar2d(jpi,jpj,A_infofile_kind(p,NTIMEK)) )

     ! loop over the variables
     do d=1,n_bfmvar
        ! get the info from the variable
        call handle_err( nfmpi_inq_var(A_infofile(p,FID), bfmvarid(d), varname, vartype, varndims, vardimids, varnatts), &
             errstring="Inquire variable: "//trim(varname) )
        call handle_err( nfmpi_inq_dimlen(A_infofile(p,FID), vardimids(1), vardimlen(1)) )
        call handle_err( nfmpi_inq_dimlen(A_infofile(p,FID), vardimids(2), vardimlen(2)) )
        if ( vardimlen(2) /= A_infofile_kind(p,NTIMEK) ) then
           write(tmp_str,*)  "Var "//trim(varname)//" not equal times: ", vardimlen(2), "/=", A_infofile_kind(p,NTIMEK)
           call handle_err(123,tmp_str)
        end if

        step_start_arr2 = (/ 1, 1 /)
        step_count_arr2 = (/ vardimlen(1), vardimlen(2) /)
        allocate(A_chunk(p,d)%p( vardimlen(1), vardimlen(2) ))

#ifdef DEBUG
        !$ nthread = omp_get_thread_num()
        write(*, '(A,I3,A,I3,A,I2,A,I2,A,I5,A,I5 )') " Thr: ", nthread, " Fid: ", p, " var: "//trim(varname)//" type: ", vartype, " ndims: ", varndims , " DIMLEN: ", vardimlen(1), " - ", vardimlen(2)
#endif

        ! read all time stamp of chunk data
        !$OMP CRITICAL
        call handle_err( nfmpi_get_vara_real_all(A_infofile(p,FID), bfmvarid(d), step_start_arr2, step_count_arr2, A_chunk(p,d)%p ), &
             errstring="Get variable: "//trim(varname) )
        !$OMP END CRITICAL

        ! inquire ID in the output file
        call handle_err( nfmpi_inq_varid(ncid, varname, vardID) )

        !write the variable content to memory
        if ( vardimlen(1) == A_infofile_kind(p,LENOCEK) ) then ! 3D variable
           bfmvar3d=NF_FILL_REAL
           noce = 1
           do k=1,jpk
              do j=1,jpj
                 do i=1,jpi
                    if (A_mask(p)%p(i,j,k) > 0.0_RLEN) then
                       bfmvar3d(i,j,k,:) = A_chunk(p,d)%p(noce,:)
                       noce = noce+1
                    end if
                 end do
              end do
           end do

           step_start_arr4 = (/ Istart, Jstart, 1, 1 /)
           step_count_arr4= (/ Icount, Jcount, jpk, A_infofile_kind(p,NTIMEK) /)
           !$OMP CRITICAL
           call handle_err(  &
                nfmpi_put_vara_real_all(ncid, vardID, step_start_arr4, step_count_arr4, bfmvar3d(iniI:jpi-finI,iniJ:jpj-finJ,:,:)), &
                errstring="Put 3D domain: "//A_procname(p)//" variable: "//trim(varname) )
           !$OMP END CRITICAL                      
        else
           bfmvar2d=NF_FILL_REAL
           noce = 1
           do j=1,jpj
              do i=1,jpi
                 if (A_mask(p)%p(i,j,1) > 0.0_RLEN) then
                    bfmvar2d(i,j,:) = A_chunk(p,d)%p(noce,:)
                    noce = noce+1
                 end if
              end do
           end do

           step_start_arr3 = (/ Istart, Jstart, 1 /)
           step_count_arr3= (/ Icount, Jcount, A_infofile_kind(p,NTIMEK) /)
           !$OMP CRITICAL
           call handle_err( & 
                nfmpi_put_vara_real_all(ncid, vardID, step_start_arr3, step_count_arr3, bfmvar2d(iniI:jpi-finI,iniJ:jpj-finJ,:)), &
                errstring="Put 2D domain: "//A_procname(p)//" variable: "//trim(varname))
           !$OMP END CRITICAL
        end if

        deallocate(A_chunk(p,d)%p)
     end do

     deallocate(bfmvar3d)
     deallocate(bfmvar2d)


#ifdef DEBUG
     !$ nthread = omp_get_thread_num()
     !$ WRITE(*,*) 'END OMP thread: ', nthread, ' ID: ', p
#endif

  end do
  !$OMP END PARALLEL DO


  ! write global grid specifications 
  ! Oce-land points mask
  if (ln_mask) then
     call handle_err( nfmpi_inq_varid(ncid, "mask", vardIDmask), errstring="inquiring variable: mask" )
     step_start_arr3 = (/ 1, 1, 1 /)
     step_count_arr3 = (/ jpiglo, jpjglo, jpk /)
     call handle_err( nfmpi_put_vara_real_all(ncid, vardIDmask, step_start_arr3, step_count_arr3 , real(maskglo,4)), &
          errstring="Writing: mask")
  endif

  ! Latitude
  call handle_err( nfmpi_inq_varid(ncid, "lat", vardIDlat) ,errstring="inquiring variable: lat")
  step_start_arr2 = (/ 1, 1 /)
  step_count_arr2 = (/ jpiglo, jpjglo /)
  call handle_err( nfmpi_put_vara_real_all(ncid, vardIDlat, step_start_arr2, step_count_arr2 , real(latglo,4)), &
       errstring="Writing: lat")

  ! Longitude
  call handle_err( nfmpi_inq_varid(ncid, "lon", vardIDlon) ,errstring="inquiring variable: lon")
  step_start_arr2 = (/ 1, 1 /)
  step_count_arr2 = (/ jpiglo, jpjglo /)
  call handle_err( nfmpi_put_vara_real_all(ncid, vardIDlon, step_start_arr2, step_count_arr2 , real(longlo,4)), &
       errstring="Writing: lon")

  ! Depth levels
  call handle_err( nfmpi_inq_varid(ncid, "depth", vardIDdep))
  call handle_err( nfmpi_begin_indep_data(ncid) )
  call handle_err( nfmpi_put_var_real(ncid, vardIDdep, real(depth,4)), &
       errstring="Writing: depth")
  call handle_err( nfmpi_end_indep_data(ncid) )



#ifdef DEBUG
  WRITE(*,*) 'Closing Files...'
#endif

  ! close netcdf chunk files
  do p=1,jpnij
     call handle_err(nfmpi_close(A_infofile(p,FID)))
  end do

  !close the output file
  call handle_err(nfmpi_close(ncid))

#ifdef DEBUG
  WRITE(*,*) 'Deallocating in Merge Vars...'
#endif

  deallocate(maskglo)
  deallocate(longlo)
  deallocate(latglo)
  deallocate(depth)
  do p=1,jpnij
     if ( allocated(A_oceanpoint(p)%p) ) deallocate(A_oceanpoint(p)%p)
     if ( allocated(A_mask(p)%p) ) deallocate(A_mask(p)%p)
     if ( allocated(A_lat(p)%p) ) deallocate(A_lat(p)%p)
     if ( allocated(A_lon(p)%p) ) deallocate(A_lon(p)%p)
  end do
  deallocate(A_oceanpoint)
  deallocate(A_mask)
  deallocate(A_lat)
  deallocate(A_lon)
  deallocate(A_chunk)
  deallocate(A_infofile)
  deallocate(A_infofile_kind)
  deallocate(A_procname)
  deallocate(A_fname)
  deallocate(A_dimlen)
  deallocate(A_dimids)

end subroutine merge_vars
