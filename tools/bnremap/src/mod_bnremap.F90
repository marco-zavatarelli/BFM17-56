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

module mod_bnremap

  use netcdf
  implicit none
  integer,    parameter, public :: RLEN=selected_real_kind(12,307)
  real(RLEN), parameter, public :: ZERO=0.0_RLEN
  real(RLEN), parameter, public :: ONE=1.0_RLEN

   integer, public ::               &  
      ntime     ,                   &  !: number of timesteps
      jpi       ,                   &  !: number of grid points along i 
      jpj       ,                   &  !: number of grid points along j 
      jpk                              !: number of grid points along k

   ! masks
   real, public, allocatable, dimension(:,:)   :: lat, lon 
   real, public, allocatable, dimension(:,:,:) :: mask

   ! no. of ocean points in the whole domain (volume and surface)
   integer,public           :: ocepoints, srfpoints

   ! NetCDF IDs  of variables to be remapped
   integer,public                                  :: n_bfmvar
   integer,allocatable,dimension(:),public         :: bfmvarid

   character(LEN=170),public            :: in_fname, out_fname, meshmask
   character(LEN=3),public              :: newBFMname
   character(LEN=4),public              :: oldBFMname
   character(LEN=NF90_MAX_NAME),public  :: maskname,latname,lonname,ocepointname
   logical,public                       :: old_version = .FALSE.
   logical,public                       :: mesh_flag = .FALSE.
   integer,public                       :: IDmaskfile,IDncin,IDncout
   
   public 

   contains 

   subroutine handle_err(iret,errstring)
   use netcdf
   implicit none
     integer,intent(in)  :: iret
     character(len=*),optional,intent(in) :: errstring
     if (iret .ne. NF90_NOERR) then
       write(*,*) "====== NetCDF Error ======"
       if (present(errstring)) write(*,*) errstring
       write(*,*) NF90_STRERROR(iret)
       stop "stop in function handle_err"
     endif
   end subroutine handle_err
   
end module mod_bnremap
