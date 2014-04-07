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

program bnmerge

  use mod_bnmerge, ONLY: GET_ARGUMENTS, chunk_fname, bfm_restart, out_fname, do_restart, do_output, tick, tock
  use create_output, ONLY: create_output_init
  use merge_vars, ONLY: merge_vars_init

  implicit none
#ifdef PARAL
  include 'mpif.h'
  character(len=MPI_MAX_PROCESSOR_NAME) name 
  integer err, rank, nprocs, namelen
#endif
  real(8)    :: calctime=0.
  integer :: calc=0

#ifdef PARAL
  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, err)
  call MPI_GET_PROCESSOR_NAME(name, namelen, err)
  write(*,*) 'Process: ', trim(name), ' Rank:', rank, ' Nprocs:', nprocs
#endif

  call tick(calc)

  call GET_ARGUMENTS
  call read_input

  write(*,*) "Starting bnmerge..."
  write(*,*)

  if( chunk_fname .NE. "" ) then
     do_output=.TRUE.
  else
     do_output=.FALSE.
  end if

  if( bfm_restart .NE. "" ) then
     do_restart=.TRUE.
  else
     do_restart=.FALSE.
  end if
  
  if( TRIM(out_fname) .EQ. "" .and. do_output ) then
     out_fname = chunk_fname
     write(*,*) "out_fname not provided. Use chunk_fname instead for data output name."
     write(*,*)
  end if

  if (do_output)  write(*,*) "Output data file is: ", TRIM(out_fname),'.nc'
  if (do_restart) write(*,*) "Output restart file is: ", TRIM(bfm_restart),'.nc'
  write(*,*)

#ifdef DEBUG
    write(*,*) "Merge Output  files? ", do_output
    write(*,*) "Merge Restart files? ", do_restart
#endif

  call create_output_init
  call merge_vars_init

  write(*,*) 
  write(*,*) 'End bnmerge'
  write(*,*) 

!#ifdef PARAL
!  CALL MPI_Finalize(err)
!#endif

  calctime = tock(calc)
  print *,'Timing summary'
  print *,'Calc: ', calctime

end program bnmerge
