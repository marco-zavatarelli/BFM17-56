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

  use mod_bnmerge
  integer,parameter :: inum=199
  integer           :: iost,n
  character(LEN=120) :: dummyline,layout
  integer,parameter    :: namlst=10,unit=11
  namelist /bnmerge_nml/ inp_dir,out_dir,chunk_fname,layout,ln_mask,var_save

      var_save="NotVar"
      ! Reading directory names and file name specification
      open(namlst,file=trim(cf_nml_bnmerge),action='read',status='old',err=99)
      read(namlst,nml=bnmerge_nml,err=98)
      close(namlst)
      ! read processor layout from layout.dat file 
      open ( UNIT=inum, FILE=layout, FORM='FORMATTED', ACCESS='SEQUENTIAL',   &
               STATUS='UNKNOWN', ERR=100, IOSTAT=iost)
      if ( iost == 0 ) then
               write (*,*) '     file   = ', trim(layout),' opened correctly'
               write (*,*) '     unit   = ', inum
               write (*,*)
      end if
      read (inum,'(6i8)',iostat=iost) jpnij,jpi,jpj,jpk,jpiglo,jpjglo
      ! get rid of the problem with different forms of the layout.dat file
      if ( iost /= 0 ) read (inum,'(6i8)',iostat=iost) jpnij,jpi,jpj,jpk,jpiglo,jpjglo
      read (inum,'(a)') dummyline
      allocate (nimppt(jpnij))
      allocate (njmppt(jpnij))
      allocate (ibonit(jpnij))
      allocate (ibonjt(jpnij))
      allocate (nlcit(jpnij))
      allocate (nlcjt(jpnij))
      allocate (nldit(jpnij))
      allocate (nldjt(jpnij))
      allocate (nleit(jpnij))
      allocate (nlejt(jpnij))
      do n = 1, jpnij
         read (inum,'(9i5)') jn, nlcit(jn), nlcjt(jn), &
                             nldit(jn), nldjt(jn), &
                             nleit(jn), nlejt(jn), &
                             nimppt(jn), njmppt(jn)
      end do
      close(inum)
      write (*,*) ' === ',trim(layout),' === '
      write (*,'(a)') '   jpnij     jpi     jpj     jpk  jpiglo  jpjglo'
      write (*,'(6i8)') jpnij,jpi,jpj,jpk,jpiglo,jpjglo
      write(*,*)
      return

100   continue
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
99    write (*,*) 'I could not open the namelist: ', trim(cf_nml_bnmerge); write (*,*)
      stop
98    write (*,*) 'I could not read the namelist: ', trim(cf_nml_bnmerge); write (*,*)
      stop

end subroutine read_input
