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

subroutine read_input

  use mod_bnmerge, ONLY : jpkglo,jpiglo,jpjglo, jpnij, &
       nimppt, njmppt, &
       nlcit , nlcjt, &
       chunk_fname, bfm_restart, out_fname, inp_dir, out_dir, ln_mask, var_save, cf_nml_bnmerge, &
       nc_compres, nc_shuffle, nc_deflate, nc_defllev

  implicit none

  integer :: jn
  integer :: jpi,jpj,jpk
  integer,parameter :: inum=199
  integer           :: iost,n
  character(LEN=120) :: dummyline,layout
  integer,parameter    :: namlst=10,unit=11
  integer, allocatable, dimension(:) :: nldit, nldjt, nleit, nlejt
  namelist /bnmerge_nml/ chunk_fname,bfm_restart,out_fname,inp_dir,out_dir,layout,ln_mask,var_save, &
                         nc_compres,nc_shuffle,nc_defllev
  ! Reading directory names and file name specification
  open(namlst,file=trim(cf_nml_bnmerge),action='read',status='old',err=99)
  read(namlst,nml=bnmerge_nml,err=98)
  close(namlst)
  ! Set use of NetCDF compression for output file
  if ( nc_defllev > 0 ) nc_deflate=1
  if ( nc_deflate == 1 ) then
    write (*,*) 'bnmerge NetCDF data compression with this setup:.'
    write (*,*) " - Shuffling (0=off/1=on) :  ",nc_shuffle
    if ( nc_deflate == 1 ) then
      write (*,*) " - Data deflation active with compression level (0-9) : ",nc_defllev
    else
      write (*,*) " - Data deflation is not active."
    endif
    write (*,*) ''
  else
    write (*,*) 'bnmerge is NOT using NetCDF data compression.'
    write (*,*) ''
  endif

  ! read processor layout from layout.dat file 
  open ( UNIT=inum, FILE=layout, FORM='FORMATTED', ACCESS='SEQUENTIAL',   &
       STATUS='UNKNOWN', ERR=100, IOSTAT=iost)
#ifdef DEBUG
  if ( iost == 0 ) then
     write (*,*) '     file   = ', trim(layout),' opened correctly'
     write (*,*) '     unit   = ', inum
     write (*,*)
  end if
#endif
  read (inum,'(6i8)',iostat=iost) jpnij,jpi,jpj,jpk,jpiglo,jpjglo

  ! get rid of the problem with different forms of the layout.dat file
  if ( iost /= 0 ) read (inum,'(6i8)',iostat=iost) jpnij,jpi,jpj,jpk,jpiglo,jpjglo
  read (inum,'(a)') dummyline
  allocate (nimppt(jpnij))
  allocate (njmppt(jpnij))
  allocate (nldit(jpnij))
  allocate (nldjt(jpnij))
  allocate (nleit(jpnij))
  allocate (nlejt(jpnij))
  allocate (nlcit(jpnij))
  allocate (nlcjt(jpnij))
  do n = 1, jpnij
     read (inum,'(9i5)') jn, nlcit(jn), nlcjt(jn), &
          nldit(jn), nldjt(jn), &
          nleit(jn), nlejt(jn), &
          nimppt(jn), njmppt(jn)
  end do
  close(inum)
  jpkglo = jpk

  deallocate(nldit, nldjt, nleit, nlejt)

  write (*,*) ' === ',trim(layout),' === '
  write (*,'(a)') '   jpnij     jpi     jpj     jpk  jpiglo  jpjglo'
  write (*,'(6i8)') jpnij,jpi,jpj,jpk,jpiglo,jpjglo
  write(*,*)
  return

100 continue
  if ( iost /= 0 ) then
     write (*,*)
     write (*,*) ' ===>>>> : problem opening file: ', trim(layout)
     write (*,*) ' =======   ===  '
     write (*,*) '           unit   = ', inum
     write (*,*) '           iostat = ', iost
     write (*,*) '           verify the file '
     write (*,*)
     stop 'Stop in read_input'
  endif
99 write (*,*) 'I could not open the namelist: ', trim(cf_nml_bnmerge); write (*,*)
  stop
98 write (*,*) 'I could not read the namelist: ', trim(cf_nml_bnmerge); write (*,*)
  stop

end subroutine read_input
