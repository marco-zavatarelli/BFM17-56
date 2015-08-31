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
                  NO_STATES,Depth,D3STATE,EPR
#ifdef INCLUDE_BEN
   use mem, only: NO_D2_BOX_STATES_BEN, D2STATE_BEN, &
                  NO_BOXES_Z_BEN, NO_BOXES_BEN, NO_STATES_BEN
   use api_bfm, only : D2STATE_BEN_tot, D2STATEB_BEN
#endif
#ifdef INCLUDE_SEAICE
   use mem, only: NO_D2_BOX_STATES_ICE, D2STATE_ICE, &
                  NO_BOXES_Z_ICE, NO_BOXES_ICE, NO_STATES_ICE
   use api_bfm, only : D2STATE_ICE_tot, D2STATEB_ICE
#endif
   use mem, only: Volume, Area, Area2d
   use mem, only: ppO2o,ppN1p,ppN3n,ppN4n,ppN5s
#ifdef INCLUDE_PELCO2
   use mem,        only: ppO3c, ppO3h, ppN6r
   use mem_CO2,    only: AtmCO2, AtmSLP, AtmTDP
#endif
   use global_mem, only: RLEN,ZERO,LOGUNIT,NML_OPEN,NML_READ, &
                         error_msg_prn,ONE, bfm_lwp
   use constants,  only: SEC_PER_DAY
   use api_bfm, only: ZEROS, SEAmask, BOTmask, SRFmask, &
        btmp1D, rtmp3Da, rtmp3Db, &
        var_names, bfm_init, &
        BOTindices,SRFindices, stPelStateS, &
        InitVar, bio_setup, in_rst_fname, out_rst_fname, save_delta, time_delta, &
        lat_len, lon_len, out_delta, out_fname, out_title, parallel_rank, &
        D3STATEB, D3STATE_tot, &
        find, update_save_delta, init_bfm, &
        ocepoint, surfpoint, botpoint
   use netcdf_bfm, only: init_netcdf_bfm,init_save_bfm
   use netcdf_bfm, only: read_rst_bfm,read_rst_bfm_glo
   use time,       only: bfmtime, julian_day, calendar_date
   use init_var_bfm_local
   ! NEMO modules
   USE trcnam_trp, only: ln_trczdf_exp,ln_trcadv_cen2,ln_trcadv_tvd
   use trc, only : ln_trc_sbc, ln_trc_ini, ln_trc_obc, ln_trc_cbc, &
        ln_top_euler, areatot, cvol, &
        trn
   use oce_trc, only : ln_qsr_bio, nn_dttrc, &
        glob_sum
   use iom_def, only:jpdom_data
   use par_oce, only: jpi, jpj, jpk, &
        jpnij, &
        jpiglo, jpjglo, jpk
   use trc_oce, only: etot3
   use trcdta, only : sf_trcdta, nb_trcdta, ntra, &
        n_trc_index, rf_trfac, trc_dta, trc_dta_init
   use trcbc, only : trc_bc_init
   use dom_oce
   use lib_mpp, only : lk_mpp, mpp_max, mpp_min
   use in_out_manager, only: numout, nitend, nit000, lwp
   use c1d, only : lk_c1d

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
   integer    :: status, yy, mm, dd, hh, nn
   integer,parameter    :: namlst=10,unit=11
   !integer,allocatable  :: ocepoint(:),surfpoint(:),botpoint(:),msktmp(:,:)
   integer,allocatable  :: msktmp(:,:)
   logical,allocatable  :: mask1d(:)
   integer              :: nc_id ! logical unit for data initialization
   character(len=40)    :: thistime
   real(RLEN)           :: julianday
   REAL(RLEN)           :: ztraf, zmin, zmax, zmean, zdrift
!EOP
!-----------------------------------------------------------------------
!BOC
   !-------------------------------------------------------
   ! Initial time
   !-------------------------------------------------------
   write(bfmtime%datestring,'(I4.4,a,I2.2,a,I2.2)') nyear,'-',nmonth,'-',nday
   write(bfmtime%date0,'(I4.4,I2.2,I2.2)') nyear,nmonth,nday
   call julian_day(nyear,nmonth,nday,0,0,julianday)
   bfmtime%time0    = julianday
   bfmtime%timeEnd  = julianday + ( ( REAL(nitend - nit000, RLEN) ) * rdt ) / SEC_PER_DAY
   bfmtime%step0    = nit000 - 1
   bfmtime%timestep = rdt
   bfmtime%stepnow  = nit000 - 1
   bfmtime%stepEnd  = nitend
   call calendar_date(bfmtime%timeEnd,yy,mm,dd,hh,nn)
   write(bfmtime%dateEnd,'(I4.4,I2.2,I2.2)') yy,mm,dd
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

   !-------------------------------------------------------
   ! assign logical Land-Sea mask
   ! in the 1D case only the central column is used by the BFM
   !-------------------------------------------------------
   if (lk_c1d) then
      where (tmask(2,2,:) > ZERO)         
        SEAmask(2,2,:) = .TRUE.
      elsewhere
        SEAmask(2,2,:) = .FALSE.
      end where
   else
      where (tmask > ZERO)         
        SEAmask = .TRUE.
      elsewhere
        SEAmask = .FALSE.
      end where
   end if

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
   NO_BOXES_Z_BEN  = 1
   NO_BOXES_BEN = NO_BOXES_XY * NO_BOXES_Z_BEN
   NO_STATES_BEN = NO_BOXES_BEN * NO_D2_BOX_STATES_BEN
#endif
#ifdef INCLUDE_SEAICE
   NO_BOXES_Z_ICE  = 1
   NO_BOXES_ICE = NO_BOXES_XY * NO_BOXES_Z_ICE
   NO_STATES_ICE = NO_BOXES_ICE * NO_D2_BOX_STATES_ICE
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

   allocate(msktmp(jpi,jpj)); msktmp = 0
   msktmp = unpack(surfpoint,SRFmask(:,:,1),0)
   allocate(botpoint(NO_BOXES_XY))
   botpoint = pack(SPREAD(msktmp,DIM=3,Ncopies=jpk),BOTmask)
   deallocate(mask1d,msktmp)

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
#ifdef INCLUDE_SEAICE
   allocate(D2STATE_ICE_tot(NO_D2_BOX_STATES_ICE))
#endif
#ifdef INCLUDE_BEN
   allocate(D2STATE_BEN_tot(NO_D2_BOX_STATES_BEN))
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
   call init_var_bfm(bio_setup)

   !---------------------------------------------
   ! Set output stepping
   !---------------------------------------------
   save_delta = bfmtime%step0
   call update_save_delta(out_delta,save_delta,time_delta)
   IF (MOD(save_delta, nn_dttrc) /= 0) THEN
      if (bfm_lwp) write(LOGUNIT,*) 'BFM : output time step must be a multiple value of the sub-stepping (nn_dttrc).'
      stop 'Mismatch of output timestep and sub-stepping. Details in bfm.log'
   ENDIF
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
   ! Water column pressure
   ! (need better approximation to convert from m to dbar)
   EPR = pack(fsdept(:,:,:),SEAmask)

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
      if (bfm_lwp) then
         LEVEL1 ' '
         LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
         LEVEL1 '              BFM INITIALIZATION               '
         LEVEL1 ' '
         write(LOGUNIT,157) 'Init', 'Unif', 'Filename    ', 'Var', 'Anal. Z1', 'Anal. V1', &
              'Anal. Z2', 'Anal. V2', 'OBC', 'SBC', 'CBC'
      endif
      do m = 1,NO_D3_BOX_STATES
         select case (InitVar(m) % init)
         case (0) 
            InitVar(m)%unif = minval( D3STATE(m,:) )
         case (1) ! Analytical profile
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
            ! mapping index
            ll = n_trc_index(m)
            call trc_dta(nit000,sf_trcdta(ll),rf_trfac(ll))
            D3STATE(m,:)  = pack(sf_trcdta(ll)%fnow(:,:,:),SEAmask)
            InitVar(m)%filename=sf_trcdta(ll)%clname
         end select
         Initvar(m)%varname=var_names(m)
         if (bfm_lwp) write(LOGUNIT, 158) InitVar(m)
      end do
      ! Initialize internal constitutents of functional groups
      call init_organic_constituents() 
   end if
   if (bfm_lwp) then 
         LEVEL1 ' '
         LEVEL1 '         BFM INITIALIZATION ... DONE!          '
         LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
         LEVEL1 ' '
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
   ! READ restart file(s) and overwrite previous initialisation
   !-------------------------------------------------------
   if (bfm_init == 1) call read_rst_bfm(in_rst_fname)
   if (bfm_init == 2) call read_rst_bfm_glo(in_rst_fname, narea, jpnij, &
        jpiglo, jpjglo, jpk, &
        nlcit, nlcjt, &
        nldit, nldjt, &
        nleit, nlejt, &
        nimppt, njmppt, &
        SEAmask )
   !-------------------------------------------------------
   ! compute and report global statistics in bfm.log
   !-------------------------------------------------------
   if (bfm_lwp) then
     LEVEL1 ' '
     LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
     LEVEL1 'VARIABLES STATISTICS AT BEGINNING OF EXPERIMENT'
     LEVEL1 ' '
   endif

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
         if (lwp) write(LOGUNIT,9000) m, trim(var_names(stPelStateS+m-1)), zmean, zmin, zmax
   end do
   if (lwp) write(LOGUNIT,*)
9000  FORMAT(' tracer :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, &
      &      '    max :',e18.10)

   if (bfm_lwp) then
     LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
     LEVEL1 ' '
   endif
   !-------------------------------------------------------
   ! Initialise DATA output netcdf file(s)
   !-------------------------------------------------------
   call calcmean_bfm(INIT)
   if (lk_c1d) then
      call init_netcdf_bfm(out_fname,TRIM(bfmtime%datestring), &
             TRIM(out_title) , 0 ,                     &
             lat=gphit(2,2),lon=glamt(2,2),            &
             roceanpoint=fsdept(1,1,1:NO_BOXES),       &
             rsurfacepoint=fsdept(1,1,1),              &
             rbottompoint=fsdept(1,1,NO_BOXES),        &
             column=.true.)
   else
      call init_netcdf_bfm(out_fname,TRIM(bfmtime%datestring), &
             TRIM(out_title) , 0 ,                     &
             lat2d=gphit,lon2d=glamt,z=fsdept(1,1,:),  &
             oceanpoint=ocepoint,                      &
             surfacepoint=surfpoint,                   &
             bottompoint=botpoint,                     &
             mask3d=tmask)
   end if
   call init_save_bfm

   !---------------------------------------------
   ! Allocate and initialise additional
   ! integration arrays
   ! Initialise prior time step for leap-frog
   !---------------------------------------------
   allocate(D3STATEB(NO_D3_BOX_STATES,NO_BOXES))
   D3STATEB = D3STATE

#ifdef INCLUDE_SEAICE
      allocate(D2STATEB_ICE(NO_D2_BOX_STATES_ICE,NO_BOXES_XY))
      D2STATEB_ICE = D2STATE_ICE
#endif


#ifdef INCLUDE_BEN
      allocate(D2STATEB_BEN(NO_D2_BOX_STATES_BEN,NO_BOXES_XY))
      D2STATEB_BEN = D2STATE_BEN
#endif

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
 
   ! Initialize boundary conditions  
   call trc_bc_init(NO_D3_BOX_STATES)

#ifdef INCLUDE_PELCO2
   ! control consistency between namelists setting
   if (AtmCO2%init .eq. 4 .and. .NOT. ln_trc_sbc(ppO3c)) then
     LEVEL1 'CO2 data from Nemo not available in surface BC &
          for O3c (check namelist_top and BFM_General).'
     stop
   endif   
   if (AtmSLP%init .eq. 4 .and. .NOT. ln_trc_sbc(ppO3h)) then
     LEVEL1 'Sea Level Pressure data from Nemo not available in surface BC &
          for O3h (check namelist_top and BFM_General).'
     stop
   endif
   if (AtmTDP%init .eq. 4 .and. .NOT. ln_trc_sbc(ppN6r)) then
     LEVEL1 'Dew Point Temperature data from Nemo not available in surface BC &
          for N6r (check namelist_top and BFM_General).'
     stop
   endif
#endif

   LEVEL1 ' '
   LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
   LEVEL1 '               EXPERIMENT START                '
   LEVEL1 '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
   LEVEL1 ' '

   return

157 FORMAT(a4, 1x, a10  , 1x, a25, 1x, a3, 1x, a10  , 1x, a10  , 1x, a10  , 1x, a10  , 3x, a3, 3x, a3, 3x, a3)
158 FORMAT(i4, 1x, E10.3, 1x, a25, 1x, a3, 1x, E10.3, 1x, E10.3, 1x, E10.3, 1x, E10.3, 3x, L3, 3x, L3, 3x, L3)
   end subroutine trc_ini_bfm
!EOC

