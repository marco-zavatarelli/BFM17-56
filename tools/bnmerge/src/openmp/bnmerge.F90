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


subroutine tick(t)
  integer, intent(OUT) :: t
  call system_clock(t)
end subroutine tick

! returns time in seconds from now to time described by t
real function tock(t)
  integer, intent(in) :: t
  integer :: now, clock_rate

  call system_clock(now,clock_rate)

  tock = real(now - t)/real(clock_rate)
end function tock


program bnmerge
  use mod_bnmerge, ONLY: GET_ARGUMENTS
  use mpi

  integer err, rank, nprocs, namelen
  character(len = MPI_MAX_PROCESSOR_NAME) :: name
  real :: calctime=0

  call tick(calc)

  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, err)
  call MPI_GET_PROCESSOR_NAME(name, namelen, err)

  write(*,*) 'Process: ', trim(name), ' Rank:', rank, ' Nprocs:', nprocs

  call GET_ARGUMENTS
  call read_input

  call create_outputfile
  call merge_vars

  WRITE(*,*) 'End bnmerge'

  CALL MPI_Finalize(err)

  calctime = tock(calc)
  print *,'Timing summary'
  print *,'Calc: ', calctime


end program bnmerge
