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

module mod_bnmerge

  integer,    parameter, public :: RLEN=selected_real_kind(12,307)
  real(RLEN), parameter, public :: ZERO=0.0_RLEN
  real(RLEN), parameter, public :: ONE=1.0_RLEN

   integer, public ::               &  
      jpi       ,                   &  !: number of grid points along i (proc)
      jpj       ,                   &  !: number of grid points along j (proc)
      jpk       ,                   &  !: number of grid points along k (proc)
      jpiglo    ,                   &  !: number of grid points along i (global domain)
      jpjglo    ,                   &  !: number of grid points along j (global domain)
      jpni      ,                   &  !: number of processors following i 
      jpnj      ,                   &  !: number of processors following j
      jpnij                            !: nb of local domain = nb of processors 
      !                                !  ( <= jpni x jpnj )

   integer, public, allocatable, dimension(:) ::   &  !:
      nimppt, njmppt,  &  !: i-, j-indexes for each processor
      ibonit, ibonjt,  &  !: i-, j- processor neighbour existence
      nlcit , nlcjt,   &  !: dimensions of every subdomain
      nldit , nldjt,   &  !: first, last indoor index for each i-domain
      nleit , nlejt       !: first, last indoor index for each j-domain

   ! masks
   real, public, allocatable, dimension(:,:)                 :: latglo, longlo ! FLOAT
   real(RLEN), public, allocatable, target, dimension(:,:,:) :: maskglo

   ! no. of ocean points in the whole domain (volume and surface)
   integer,public           :: ocepoints, srfpoints

   ! NetCDF IDs  of variables to be merged
   integer,public                                  :: n_bfmvar
   integer,allocatable,dimension(:),public         :: bfmvarid

   ! namelist variables
   character(LEN=100) :: inp_dir, out_dir, out_fname
   integer,parameter  :: NSAVE=120      ! Maximum no variables which can be saved
   character(len=64),dimension(NSAVE):: var_save
   logical :: ln_grid=.FALSE.

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
   
end module mod_bnmerge
