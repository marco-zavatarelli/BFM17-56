! BFM_NEMO chlsat
!    Copyright (C) 2014 Marcello Vichi (marcello.vichi@bo.ingv.it)
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

module mod_chlsat

   integer,    parameter, public :: RLEN=selected_real_kind(12,307)
   real(RLEN), parameter, public :: ZERO=0.0_RLEN
   real(RLEN), parameter, public :: ONE=1.0_RLEN

   integer, public :: jpi,jpj,jpk,ntime

   ! masks
   real, public, allocatable, dimension(:)             :: depth ! FLOAT
   real, public, allocatable, dimension(:,:)           :: lat, lon ! FLOAT
   real(RLEN), public, allocatable, dimension(:,:,:)   :: mask,e3t
   real(RLEN), public, allocatable, dimension(:,:,:,:) :: chla,eps
   real(RLEN), public, allocatable, dimension(:,:,:)   :: chlsat_od,chlsat_01,chlsat_001
   real(RLEN), public, allocatable, dimension(:,:,:)   :: ezd_od,ezd_01,ezd_001
   real(RLEN), public, allocatable, dimension(:,:,:,:) :: gpp,rsp
   real(RLEN), public, allocatable, dimension(:,:,:)   :: gpp_01,gpp_001
   real(RLEN), public, allocatable, dimension(:,:,:)   :: npp_01,npp_001

   ! namelist variables
   character(LEN=400) :: cf_nml='chlsat.nml'     ! namelist name
   character(LEN=100) :: inp_dir, out_dir, out_fname, chla_fname,eps_fname
   character(LEN=100) :: mask_fname
   character(LEN=20)  :: chla_name,eps_name
   real(RLEN)         :: tolerance=ZERO
   logical            :: compute_eps=.FALSE.,compute_chlsat=.TRUE.,compute_intpp=.FALSE.
   real(RLEN)         :: p_eps0=0.0435_RLEN, p_epsChla=0.03
   character(LEN=100) :: gpp_fname,rsp_fname
   character(LEN=20)  :: gpp_name,rsp_name

   public 
   contains 
!
!   ------------------------------------------------------------------------------    
!    Compute photic zone depth and average chlorophyll
!    given a certain optical threshold tau (e.g. tau = -log(0.01))
!   ------------------------------------------------------------------------------
!
   subroutine chlcalc(tau,ezd,chlsat)
   use netcdf
   implicit none
   real(RLEN),                         intent(in)  :: tau
   real(RLEN),dimension(jpi,jpj,ntime),intent(out) :: ezd,chlsat
   integer                                         :: i,j,k,t
   real(RLEN)                                      :: store,dep
   real(RLEN)                                      :: expterm,optlen,norm

   ezd(:,:,:) = ZERO
   chlsat(:,:,:) = ZERO
   do t = 1,ntime
      do j = 1,jpj
         do i = 1,jpi
          store = ZERO
          dep = ZERO
          norm = ZERO
          if (mask(i,j,1) > ZERO) then
            do k = 1,jpk-1
               optlen = eps(i,j,k,t)*e3t(i,j,k)
               store = store + optlen
               dep = dep + e3t(i,j,k)
               if ((store .lt. tau) .and. (mask(i,j,k)>ZERO)) then
                  expterm = exp(-optlen) 
                  chlsat(i,j,t) = chlsat(i,j,t)+chla(i,j,k,t)*expterm

                  norm = norm+expterm
               else
                  ! special case if absorption is high in the first layer
                  if (k==1) then
                     ezd(i,j,t) = e3t(i,j,1)
                     chlsat(i,j,t) = chla(i,j,1,t)
                  else
                     ezd(i,j,t) = dep - e3t(i,j,k)*0.5_RLEN
                     chlsat(i,j,t) = chlsat(i,j,t)/norm
                  end if
                  exit
               end if
            end do
          else
            ezd(i,j,t) = NF90_FILL_REAL
            chlsat(i,j,t) = NF90_FILL_REAL
          end if
         end do
      end do
   end do
   end subroutine chlcalc

!
!   ------------------------------------------------------------------------------    
!    Additional parameterizations
!   ------------------------------------------------------------------------------
!

!
!   ------------------------------------------------------------------------------    
!    Integrated gross and net primary production
!    given a certain optical threshold tau (e.g. tau = -log(0.01))
!   ------------------------------------------------------------------------------
!
   subroutine intppcalc(tau,intgpp,intnpp)
   use netcdf
   implicit none
   real(RLEN),                         intent(in)  :: tau
   real(RLEN),dimension(jpi,jpj,ntime),intent(out) :: intgpp,intnpp
   integer                                         :: i,j,k,t
   real(RLEN)                                      :: store
   real(RLEN)                                      :: optlen

   intgpp(:,:,:) = ZERO
   intnpp(:,:,:) = ZERO
   do t = 1,ntime
      do j = 1,jpj
         do i = 1,jpi
          store = ZERO
          if (mask(i,j,1) > ZERO) then
            do k = 1,jpk-1
               optlen = eps(i,j,k,t)*e3t(i,j,k)
               store = store + optlen
               if ((store .lt. tau) .and. (mask(i,j,k)>ZERO)) then
                  intgpp(i,j,t) = intgpp(i,j,t)+gpp(i,j,k,t)*e3t(i,j,k)
                  intnpp(i,j,t) = intnpp(i,j,t)+(gpp(i,j,k,t)-rsp(i,j,k,t))*e3t(i,j,k)
               else
                  ! special case if absorption is high in the first layer
                  if (k==1) then
                     intgpp(i,j,t) = gpp(i,j,1,t)*e3t(i,j,1)
                     intnpp(i,j,t) = (gpp(i,j,1,t)-rsp(i,j,1,t))*e3t(i,j,1)
                  end if
                  if (intgpp(i,j,t)>NF90_FILL_REAL) then
                     intgpp(i,j,t)=NF90_FILL_REAL
                     write(*,*) 'WARNING: incompatible mask for gpp at point',i,j,k,t,gpp(i,j,1,t)
                  end if
                  if (intnpp(i,j,t)>NF90_FILL_REAL) then
                     intnpp(i,j,t)=NF90_FILL_REAL
                     write(*,*) 'WARNING: incompatible mask for rsp at point',i,j,k,t,rsp(i,j,1,t)
                  end if
                  exit
               end if
            end do
          else
            intgpp(i,j,t) = NF90_FILL_REAL
            intnpp(i,j,t) = NF90_FILL_REAL
          end if
         end do
      end do
   end do
   end subroutine intppcalc

!
!   ------------------------------------------------------------------------------    
!    Handle errors of NetCDF operations
!   ------------------------------------------------------------------------------
!
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
!
!   ------------------------------------------------------------------------------    
!    Retrieve the namelist file name and location
!   ------------------------------------------------------------------------------
!
   subroutine get_arguments()
    
   INTEGER            :: iargc, jarg
   CHARACTER(len=400) :: cr 
   CHARACTER(LEN=2), DIMENSION(1), PARAMETER :: &
        clist_opt = (/ '-f' /)     

   ! --- get arguments
   jarg = 0

   ! set the default name
   if (jarg == iargc()) then
          PRINT *, ' Input namelist not specified !!! '
          PRINT *, ' Setting the name to the default (chlsat.nml). ' 
          PRINT *, ' -----'
          PRINT *, ' To specify a different filename use the option -f, e.g. '
          PRINT *, '       chlsat -f namelist.nml'
          PRINT *, ''
          cr=trim(cf_nml)
   endif

   DO WHILE ( jarg < iargc() )
      !!
      jarg = jarg + 1
      CALL getarg(jarg,cr)
      !!
      SELECT CASE (trim(cr))
 
      CASE('-f')
        
         IF ( jarg + 1 > iargc() ) THEN
            PRINT *, 'ERROR: Missing namelist name!'; goto 123; 
         ELSE
            !!
            jarg = jarg + 1
            CALL getarg(jarg,cr)
            !!
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: ', trim(cr), ' is definitively not the name of the namelist!'
               goto 123;
            END IF
            !!
            !!
         END IF

      CASE DEFAULT
         PRINT *, ' Unrecognized input option !!! '
         PRINT *, ' Setting the name to the default (chlsat.nml). ' 
         PRINT *, ' -----'
         PRINT *, ' To specify a different filename use the option -f, e.g. '
         PRINT *, '       chlsat -f namelist.nml'
         PRINT *, '' 
         cr=trim(cf_nml)
      END SELECT       
 
   END DO

    !
    ! set the namelist name
    cf_nml=trim(cr)

    PRINT *, ''
    PRINT *, 'Namelist is: ', trim(cf_nml)   
    PRINT *, ''

   return

123 PRINT*, 'If you see this something is wrong with the specified namelist! ' 
   STOP
   end subroutine get_arguments


end module mod_chlsat
