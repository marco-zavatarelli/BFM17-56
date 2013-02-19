#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: trc_ini_bfm.F90
!
! !INTERFACE:
   subroutine trc_ini_bfm()
!
! !DESCRIPTION:
!  Initialise the BFM in NEMO
!  Main communication of array dimensions between
!  BFM and NEMO
!  Initialization of variables and netcdf output
!
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,          &
                  NO_BOXES_X, NO_BOXES_Y, NO_BOXES_Z,  &
                  NO_BOXES_XY, NO_D3_BOX_DIAGNOSS,     &
                  NO_STATES,Depth,D3STATE,PELRIVER
#ifdef INCLUDE_BEN
   use mem, only: NO_D2_BOX_STATES, NO_D2_BOX_DIAGNOSS, &
                  D2STATE
#endif
   use mem, only: Volume, Area, Area2d
   use mem, only: ppO2o,ppN1p,ppN3n,ppN4n,ppN5s
#ifdef INCLUDE_PELCO2
   use mem, only: ppO3c,ppO3h
#endif
   use global_mem, only: RLEN,ZERO,LOGUNIT,NML_OPEN,NML_READ, &
                         error_msg_prn,ONE
   use constants,  only: SEC_PER_DAY
   use api_bfm
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   use netcdf_bfm, only: init_netcdf_rst_bfm,read_rst_bfm
   use time,       only: bfmtime, julian_day 
   ! NEMO modules
   USE trcnam_trp, only: ln_trczdf_exp,ln_trcadv_cen2,ln_trcadv_tvd
   use trc
   use oce_trc
   use iom_def,    only:jpdom_data
   use iom
   use sbc_oce, only: ln_rnf
   use trc_oce, only: etot3
   use trcdta
   use trcbc
   use dom_oce, only: nyear, nmonth, nday

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   !! * Substitutions
#include "domzgr_substitute.h90"

   integer    :: i,j,k,ll,m
   integer    :: status
   integer,parameter    :: namlst=10,unit=11
   integer,allocatable  :: ocepoint(:),surfpoint(:),botpoint(:)
   logical,allocatable  :: mask1d(:)
   integer              :: nc_id ! logical unit for data initialization
   character(len=20)    :: start_time
   real(RLEN)           :: julianday
   REAL(RLEN)           :: ztraf, zmin, zmax, zmean, zdrift
!EOP
!-----------------------------------------------------------------------
!BOC
   !-------------------------------------------------------
   ! Initial time
   !-------------------------------------------------------
   write(start_time,'(I4.4,a,I2.2,a,I2.2)') nyear,'-',nmonth,'-',nday
   call julian_day(nyear,nmonth,nday,0,0,julianday)
   bfmtime%date0    = start_time
   bfmtime%time0    = julianday
   bfmtime%timeEnd  = julianday + ((nitend-nit000)* rdt) / SEC_PER_DAY
   bfmtime%step0    = nit000 - 1
   bfmtime%timestep = rdt
   bfmtime%stepnow  = nit000 - 1
   bfmtime%stepEnd  = nitend
   !-------------------------------------------------------
   ! Force Euler timestepping for the BFM
   !-------------------------------------------------------
   ln_top_euler = .TRUE.

   !-------------------------------------------------------
   ! Initial allocations
   !-------------------------------------------------------
   allocate(SEAmask(jpi,jpj,jpk))  ! allocate  masks for
   allocate(BOTmask(jpi,jpj,jpk))  ! array packing
   allocate(SRFmask(jpi,jpj,jpk))
   SEAmask = .FALSE.
   BOTmask = .FALSE.
   SRFmask = .FALSE.

   allocate(ZEROS(jpi,jpj,jpk))   ! allocate ancillary pack mask
   ZEROS = ZERO
   allocate(rtmp3Da(jpi,jpj,jpk)) ! allocate temporary 3D array
   allocate(rtmp3Db(jpi,jpj,jpk)) ! allocate temporary 3D array
   rtmp3Da = ZERO; rtmp3Db = ZERO

   where (tmask > ZERO)         ! assign logical Land-Sea mask
     SEAmask = .TRUE.
   elsewhere
     SEAmask = .FALSE.
   end where

   !-------------------------------------------------------
   ! Prepares the spatial information and the masks
   ! for the bottom and surface grid points.
   ! 3D boolean arrays with .T. at the first
   ! and last ocean layers only
   !-------------------------------------------------------
   do j = 1,jpj
      do i = 1,jpi
         rtmp3Da(i,j,:) = fse3t(i,j,:)
         rtmp3Db(i,j,:) = e1t(i,j)*e2t(i,j)
         if (SEAmask(i,j,1)) then
            SRFmask(i,j,1) = .TRUE.
            do k = 1,jpk
               if (.not.SEAmask(i,j,k)) then
                  BOTmask(i,j,k-1) =  .TRUE.
                  exit
               end if
            end do
         end if
      end do
   end do

   !---------------------------------------------
   ! Set the dimensions
   !---------------------------------------------
   NO_BOXES_X  = jpi
   NO_BOXES_Y  = jpj
   NO_BOXES_Z  = jpk
   NO_BOXES    = count(SEAmask)
   NO_BOXES_XY = count(SRFmask)
   NO_STATES   = NO_D3_BOX_STATES * NO_BOXES
#ifdef INCLUDE_BEN
   NO_STATES = NOSTATES + NO_BOXES_XY*NO_D2_BOX_STATES
#endif

   !-------------------------------------------------------
   ! Allocate and build the indices of ocean points in
   ! the 3D nemo arrays
   !-------------------------------------------------------
#ifndef USEPACK
   allocate (iwet(NO_BOXES))
   allocate (jwet(NO_BOXES))
   allocate (kwet(NO_BOXES))
   ll=0
   do k = 1, jpk
      do j = 1, jpj
         do i = 1, jpi
            if (SEAmask(i,j,k)) then
               ll=ll+1
               iwet(ll)=i
               jwet(ll)=j
               kwet(ll)=k
            endif
         enddo
      enddo
    enddo
#endif

   !-------------------------------------------------------
   ! Compressed coordinates for netcdf output
   !-------------------------------------------------------
   lon_len = NO_BOXES_X
   lat_len = NO_BOXES_Y
   allocate(ocepoint(NO_BOXES))
   allocate(mask1d(1:NO_BOXES_X*NO_BOXES_Y*NO_BOXES_Z))
   mask1d = reshape(SEAmask,(/NO_BOXES_X*NO_BOXES_Y*NO_BOXES_Z/))
   ocepoint = find(mask1d,NO_BOXES)
   allocate(surfpoint(NO_BOXES_XY))
   mask1d = reshape(SRFmask,(/NO_BOXES_X*NO_BOXES_Y*NO_BOXES_Z/))
   surfpoint = find(mask1d,NO_BOXES_XY)
   allocate(botpoint(NO_BOXES_XY))
   mask1d = reshape(BOTmask,(/NO_BOXES_X*NO_BOXES_Y*NO_BOXES_Z/))
   botpoint = find(mask1d,NO_BOXES_XY)
   deallocate(mask1d)

   !-------------------------------------------------------
   ! Prepares the array containing the indices of the
   ! elements in pelagic BFM 1D arrays that have a benthic layer
   !-------------------------------------------------------
   allocate(BOTindices(NO_BOXES_XY)); BOTindices = 0
   allocate(btmp1D(NO_BOXES))
   btmp1D = pack(BOTmask,SEAmask)
   BOTindices = find(btmp1D,NO_BOXES_XY)

   !-------------------------------------------------------
   ! Prepares the array containing the indices of the
   ! elements in pelagic BFM 1D arrays that have a surface layer
   !-------------------------------------------------------
   allocate(SRFindices(NO_BOXES_XY)); SRFindices = 0
   btmp1D = pack(SRFmask,SEAmask)
   SRFindices = find(btmp1D,NO_BOXES_XY)
   deallocate(btmp1D)

   !-------------------------------------------------------
   ! Prepares the array containing the total amount per var
   !-------------------------------------------------------
   allocate(D3STATE_tot(NO_D3_BOX_STATES))
#ifdef INCLUDE_BEN
   allocate(D2STATE_tot(NO_D2_BOX_STATES))
#endif

   !---------------------------------------------
   ! Assign the rank of the process
   ! (meaningful only with key_mpp)
   !---------------------------------------------
   parallel_rank = narea-1

   !---------------------------------------------
   ! Initialise ancillary arrays for output
   ! (also parallel initialisation is done here)
   !---------------------------------------------
   call init_bfm(namlst)

   !---------------------------------------------
   ! Initialise state variable names and diagnostics
   !---------------------------------------------
   call set_var_info_bfm

   !-------------------------------------------------------
   ! Allocate memory and give homogeneous initial values
   !-------------------------------------------------------
   ! the argument list is mandatory with BFM
   call init_var_bfm(namlst,'BFM_General.nml',unit,bio_setup)

   !-------------------------------------------------------
   ! Prepares the BFM 1D arrays containing the
   ! spatial informations (have to be done after allocation
   ! but temporary arrays allocated above)
   !-------------------------------------------------------
   Area2d = pack(rtmp3Db(:,:,1),SEAmask(:,:,1))
   Area   = pack(rtmp3Db,SEAmask)
   Volume = pack(rtmp3Da*rtmp3Db,SEAmask)
   !thickness at each sea gridpoint (Depth(NO_BOXES))
   Depth  = pack(rtmp3Da,SEAmask)

   !-------------------------------------------------------
   ! Initialization from analytical profiles or data
   ! Done if restart is not used
   !-------------------------------------------------------
   if (bfm_init /= 1) then
      ! this is done for compatibility with NEMO variables
      if (allocated(ln_trc_ini)) deallocate(ln_trc_ini)
      allocate(ln_trc_ini(NO_D3_BOX_STATES))
      ln_trc_ini(:) = .false.
      do m = 1,NO_D3_BOX_STATES
         if (InitVar(m) % init == 2) ln_trc_ini(m) = .true.
      end do
      ! initialize the data structure for input fields
      ! found in the top_namelist
      call trc_dta_init(NO_D3_BOX_STATES)
      if (lwp) then
         write(numout,*)
         write(numout,*) '           ---- BFM Initialization ----             '
         write(numout,*)
      endif
      do m = 1,NO_D3_BOX_STATES
         if (lwp) write(numout,*) 'BFM variable ',trim(var_names(stPelStateS+m-1)),' is initialized with:'
         select case (InitVar(m) % init)
         case (1) ! Analytical profile
            if (lwp) write(numout,*) '--> Analytical profile' 
            rtmp3Da = ZERO
            ! fsdept contains the model depth
            do j = 1,jpj
               do i = 1,jpi
                  call analytical_profile(jpk,fsdept(i,j,:),InitVar(m) % anz1, &
                    InitVar(m) % anv1,InitVar(m) % anz2,InitVar(m) % anv2,rtmp3Da(i,j,:))
               end do
            end do
            D3STATE(m,:)  = pack(rtmp3Da,SEAmask)
         case (2) ! from file
            if (lwp) write(numout,*) '--> Data file' 
            ! mapping index
            ll = n_trc_index(m)
            call trc_dta(nit000,sf_trcdta(ll),rf_trfac(ll))
            D3STATE(m,:)  = pack(sf_trcdta(ll)%fnow(:,:,:),SEAmask)
         case default 
            if (bfm_init==0) then
               if (lwp) write(numout,*) '--> Constant homogeneous value' 
            else
               if (lwp) write(numout,*) '--> Restart data file' 
            end if
         end select
      end do
   end if

   deallocate(rtmp3Da)
   deallocate(rtmp3Db)

!   IF( nb_trcdta > 0 .AND. .NOT.ln_trcdmp ) THEN
   IF( nb_trcdta > 0 ) THEN
      !==   deallocate data structure   ==! 
      !        data used only for initialisation)
      IF(lwp) WRITE(numout,*) 'trc_dta: deallocate data arrays as they are only used to initialize the run'
      DO ll = 1, ntra
                                       DEALLOCATE( sf_trcdta(ll)%fnow )     !  arrays in the structure
         IF( sf_trcdta(ll)%ln_tint )   DEALLOCATE( sf_trcdta(ll)%fdta )
      ENDDO
   ENDIF

   !-------------------------------------------------------
   ! initialise (new) and read (previous) restart file
   ! (override any previous initialisation)
   !-------------------------------------------------------
   call init_netcdf_rst_bfm(rst_fname)
   if (bfm_init == 1) call read_rst_bfm(rst_fname)

   !-------------------------------------------------------
   ! compute and report global statistics
   ! in ocean.output
   !-------------------------------------------------------
   D3STATE_tot(:) = ZERO
   do m = 1,NO_D3_BOX_STATES
         ! compute statistics (need to map to 3D shape for global sum)
         trn(:,:,:,1) = unpack(D3STATE(m,:),SEAmask,ZEROS)
         D3STATE_tot(m) = glob_sum( trn(:,:,:,1) * cvol(:,:,:)   )
         zmin  = minval( trn(:,:,:,1), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
         zmax  = maxval( trn(:,:,:,1), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
         if (lk_mpp) then
            call mpp_min( zmin )      ! min over the global domain                                                                
            call mpp_max( zmax )      ! max over the global domain
         end if
         zmean  = D3STATE_tot(m) / areatot
         if (lwp) write(numout,9000) m, trim(var_names(stPelStateS+m-1)), zmean, zmin, zmax
         if (lwp) write(numout,*)
   end do
9000  FORMAT(' tracer :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, &
      &      '    max :',e18.10)

   !-------------------------------------------------------
   ! initialise netcdf output
   !-------------------------------------------------------
   call calcmean_bfm(INIT)
   call init_netcdf_bfm(out_fname,TRIM(start_time),0,  &
             lat2d=gphit,lon2d=glamt,z=gdept_0,        &
             oceanpoint=ocepoint,                      &
             surfacepoint=surfpoint,                   &
             bottompoint=botpoint,                     &
             mask3d=tmask)
   call init_save_bfm

   !---------------------------------------------
   ! Allocate and initialise additional
   ! integration arrays
   ! Initialise prior time step for leap-frog
   !---------------------------------------------
   allocate(D3STATEB(NO_D3_BOX_STATES,NO_BOXES))
   D3STATEB = D3STATE

#ifdef INCLUDE_BEN
      allocate(D2STATEB(NO_D2_BOX_STATES,NO_BOXES))
      D2STATEB = D2STATE
#endif


   ! Initialise the array containing light bioshading
   !-------------------------------------------------------
   if ( ln_qsr_bio ) etot3(:,:,:) = ZERO

   ! Initialise the arrays containing external boundary data
   !-------------------------------------------------------
   if (allocated(ln_trc_obc)) deallocate(ln_trc_obc)
      allocate(ln_trc_obc(NO_D3_BOX_STATES))
      ln_trc_obc(:) = .false.
   if (allocated(ln_trc_sbc)) deallocate(ln_trc_sbc)
      allocate(ln_trc_sbc(NO_D3_BOX_STATES))
      ln_trc_sbc(:) = .false.
   if (allocated(ln_trc_cbc)) deallocate(ln_trc_cbc)
      allocate(ln_trc_cbc(NO_D3_BOX_STATES))
      ln_trc_cbc(:) = .false.
   do m = 1,NO_D3_BOX_STATES
      if (InitVar(m) % obc) ln_trc_obc(m) = .true.
      if (InitVar(m) % sbc) ln_trc_sbc(m) = .true.
      if (InitVar(m) % cbc) ln_trc_cbc(m) = .true.
   end do
   call trc_bc_init

   return

   end subroutine trc_ini_bfm
!EOC

