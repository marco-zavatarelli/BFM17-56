#include"cppdefs.h"
#define REAL_4B real(4)
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: netcdf_bfm --- Save the BFM results in NetCDF
!
! !INTERFACE:
   module netcdf_bfm
!
! !DESCRIPTION:
!  This module provides routines for saving the results using
!  NetCDF format.
!
! !USES:
   use string_functions, ONLY : replace_char
   use api_bfm, ONLY: var_names, var_units, var_long, var_ids, &
        var_ave, &
        ocepoint_len,surfpoint_len,botpoint_len, &
        lon_len, lat_len, &
        c1dim, &
        bfm_rstctl, out_dir, &
        D3ave, D2ave, &
        nc_compres,nc_shuffle,nc_deflate,nc_defllev, &
#if defined INCLUDE_SEAICE
        D2ave_ice, &
        stIceStateS, stIceDiag2dS, stIceFlux2dS, stIceStateE, &
        stIceDiag2dE, stIceFlux2dE, stIceStart, stIceEnd, &
#endif
#if defined INCLUDE_BEN
        D2ave_ben, &
        stBenStateS, stBenDiag2dS, stBenFlux2dS, stBenStateE, &
        stBenDiag2dE, stBenFlux2dE, stBenStart, stBenEnd,  &
#endif
        stStart, stEnd, stPelStateS, stPelDiagS, stPelFluxS, &
        stPelDiag2dS, stPelSurS, stPelBotS, stPelRivS, &
        stPelStateE, stPelDiagE, stPelFluxE, stPelDiag2dE, &
        stPelSurE, stPelBotE, stPelRivE, stPelStart, stPelEnd


   use mem,     only: NO_BOXES,NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z,NO_BOXES_XY,Depth
   use mem,     only: D3FLUX_FUNC
#if defined INCLUDE_SEAICE
   use mem,     only: D2FLUX_FUNC_ICE
#endif
#if defined INCLUDE_BEN
   use mem,     only: D2FLUX_FUNC_BEN
#endif

   use global_mem, only: RLEN,ZERO,LOGUNIT,bfm_lwp
   use constants, ONLY: SEC_PER_DAY
   use netcdf
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_save_bfm,save_bfm,close_ncdf,check_err
   public define_mode,new_nc_variable,set_attributes,store_data
!
! !PUBLIC DATA MEMBERS:
   !---------------------------------------------
   ! netCDF file specifications
   !---------------------------------------------
   integer,public                :: ncid_bfm
   integer,public                :: ncdf_time_unit
   ! record counter
   integer,public                :: recnum = 0
   !---------------------------------------------
   ! Dimension IDs
   !---------------------------------------------
   integer                            :: lon_dim,lat_dim,depth_dim
   integer                            :: x_dim,y_dim
   integer                            :: ocepoint_dim
   integer                            :: surfpoint_dim,botpoint_dim
   integer                            :: time_dim
   integer,allocatable                :: dims(:)
   !---------------------------------------------
   ! Coordinate variables IDs
   !---------------------------------------------
   integer          :: lon_id,lat_id,z_id,time_id
   integer          :: mask_id
   integer          :: depth_id,ocepoint_id,surfpoint_id,botpoint_id
   !---------------------------------------------
   ! Restart file ids
   !---------------------------------------------
   integer          :: lon_rid,lat_rid,z_rid,time_rid
   integer          :: mask_rid
   integer          :: depth_rid,ocepoint_rid,surfpoint_rid,botpoint_rid
   integer,public                :: ncid_rst
   integer                       :: ncid_rst_in
   integer                       :: lon_rdim,lat_rdim,depth_rdim
   integer                       :: x_rdim,y_rdim
   integer                       :: time_rdim
   integer                       :: ocepoint_rdim
   integer                       :: surfpoint_rdim,botpoint_rdim
   integer                       :: chars_rdim
   integer                       :: d3vars_rdim,d3state_rid,d3stateb_rid,d3state_name_rid,d3state_units_rid,d3state_long_rid
#ifdef INCLUDE_SEAICE
   integer                       :: d2vars_rdim_ice,d2state_rid_ice,d2stateb_rid_ice,d2state_name_rid_ice, &
        d2state_units_rid_ice,d2state_long_rid_ice
#endif
#ifdef INCLUDE_BEN
   integer                       :: d2vars_rdim_ben,d2state_rid_ben,d2stateb_rid_ben,d2state_name_rid_ben, &
        d2state_units_rid_ben,d2state_long_rid_ben
#endif
#ifdef INCLUDE_PELCO2
   integer                       :: ph_rid
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications and BFM additions: Marcello Vichi , Tomas Lovato
!
!EOP
!
! !PRIVATE DATA MEMBERS
!  variable ids
   integer, private          :: start(4),edges(4)
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the netcdf output
!
! !INTERFACE:
   subroutine init_netcdf_bfm(title,start_time,expinfo,time_unit,      &
                              lat,lon,z,dz,lat2d,lon2d,                &
                              oceanpoint,surfacepoint,bottompoint,     &
                              roceanpoint,rsurfacepoint,rbottompoint,  &
                              mask3d,column)
!
! !DESCRIPTION:
!  Prepare the netcdf output file which is finalized in init_save_bfm
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*), intent(in)                 :: title,start_time,expinfo
   integer, intent(in)                          :: time_unit
   real(RLEN), intent(in),optional                :: lat,lon
   real(RLEN), intent(in),dimension(:,:),optional :: lat2d,lon2d
   real(RLEN), intent(in),dimension(:),optional   :: z,dz
   integer, intent(in),dimension(:),optional    :: oceanpoint
   integer, intent(in),dimension(:),optional    :: surfacepoint,bottompoint
   real(RLEN), intent(in),dimension(:),optional :: roceanpoint
   real(RLEN), intent(in),optional              :: rsurfacepoint,rbottompoint
   real(RLEN),intent(in),dimension(:,:,:),optional:: mask3d
   logical, intent(in),optional                   :: column
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications: Marcello Vichi
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: ext,fname
   integer                   :: iret,ndims
   character(len=128)        :: ncdf_time_str,history
!  dimension lengths (not used yet)
   integer                   :: lon_len
   integer                   :: lat_len
   integer                   :: depth_len
!!
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'init_netcdf_bfm: create output data file(s) ...'

   !---------------------------------------------
   ! Prepare the netcdf file
   !---------------------------------------------
   ext = 'nc'
   fname = TRIM(out_dir) //'/'// TRIM(title) // '.' // ext
   LEVEL2 'Output NetCDF file is (time unit is set to seconds):'
   LEVEL2 TRIM(fname)
   call check_err(NF90_CREATE(fname,NF90_NETCDF4,ncid_bfm), fname)

   ncdf_time_unit = time_unit

   !---------------------------------------------
   ! define dimensions
   !---------------------------------------------
   if (present(lon).AND.present(lat)) then
      call check_err(NF90_DEF_DIM(ncid_bfm, 'lon', max(NO_BOXES_XY,1), lon_dim), fname)
      call check_err(NF90_DEF_DIM(ncid_bfm, 'lat', max(NO_BOXES_XY,1), lat_dim), fname)
   else if (present(lon2d).AND.present(lat2d)) then
      call check_err(NF90_DEF_DIM(ncid_bfm, 'x', NO_BOXES_X, x_dim), fname)
      call check_err(NF90_DEF_DIM(ncid_bfm, 'y', NO_BOXES_Y, y_dim), fname)
   else
      stop '### init_netcdf_bfm: lat and lon must be given'
   end if
   if (present(z).or.present(dz)) &
      call check_err(NF90_DEF_DIM(ncid_bfm, 'z', NO_BOXES_Z, depth_dim), fname)
   call check_err(NF90_DEF_DIM(ncid_bfm, 'oceanpoint', max(NO_BOXES,1), ocepoint_dim), fname)
   call check_err(NF90_DEF_DIM(ncid_bfm, 'surfacepoint', max(NO_BOXES_XY,1), surfpoint_dim), fname)
   call check_err(NF90_DEF_DIM(ncid_bfm, 'bottompoint', max(NO_BOXES_XY,1), botpoint_dim), fname)
   call check_err(NF90_DEF_DIM(ncid_bfm, 'time', NF90_UNLIMITED, time_dim), fname)   

   !---------------------------------------------
   ! define coordinate variables
   !---------------------------------------------
   ALLOCATE(dims(2))
   dims(1) = x_dim
   dims(2) = y_dim
   if (present(lon)) then
     call check_err(NF90_DEF_VAR(ncid_bfm,'lon',NF90_REAL,lon_dim,lon_id), fname)
   elseif (present(lon2d)) then
     call check_err(NF90_DEF_VAR(ncid_bfm,'lon',NF90_REAL,dims,lon_id), fname)
   end if
   if (present(lat)) then
     call check_err(NF90_DEF_VAR(ncid_bfm,'lat',NF90_REAL,lat_dim,lat_id), fname)
   elseif (present(lat2d)) then
     call check_err(NF90_DEF_VAR(ncid_bfm,'lat',NF90_REAL,dims,lat_id), fname)
   end if
   DEALLOCATE(dims)
   if (present(z)) &
      call check_err(NF90_DEF_VAR(ncid_bfm,'z',NF90_REAL,depth_dim,depth_id), fname)
   if (present(column)) then
      ! in 1D case coordinate variables are real and store the depth
      call check_err(NF90_DEF_VAR( ncid_bfm, 'oceanpoint', NF90_REAL,ocepoint_dim,ocepoint_id), fname)
      call check_err(NF90_DEF_VAR(ncid_bfm,'surfacepoint',NF90_REAL,surfpoint_dim,surfpoint_id), fname)
      call check_err(NF90_DEF_VAR(ncid_bfm,'bottompoint',NF90_REAL,botpoint_dim,botpoint_id), fname)
   else
      call check_err(NF90_DEF_VAR( ncid_bfm, 'oceanpoint', NF90_INT,ocepoint_dim,ocepoint_id), fname)
      call check_err(NF90_DEF_VAR(ncid_bfm,'surfacepoint',NF90_INT,surfpoint_dim,surfpoint_id), fname)
      call check_err(NF90_DEF_VAR(ncid_bfm,'bottompoint',NF90_INT,botpoint_dim,botpoint_id), fname)
   end if
   call check_err(NF90_DEF_VAR(ncid_bfm,'time',NF90_REAL,time_dim,time_id), fname)
   !---------------------------------------------
   ! define mask variables
   !---------------------------------------------
   if (present(mask3d)) then
      ALLOCATE(dims(3))
      dims(1) = x_dim
      dims(2) = y_dim
      dims(3) = depth_dim
      call check_err(NF90_DEF_VAR(ncid_bfm,'mask',NF90_REAL,dims,mask_id), fname)
      if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_bfm,mask_id,nc_shuffle,nc_deflate,nc_defllev))
      DEALLOCATE(dims)
   end if

   !---------------------------------------------
   ! assign attributes
   !---------------------------------------------
   !  coordinates
   call check_err(set_attributes(ncid_bfm,lon_id,units='degrees_east'), fname)
   call check_err(set_attributes(ncid_bfm,lat_id,units='degrees_north'), fname)
   if (present(z)) &
      call check_err(set_attributes(ncid_bfm,depth_id,units='meters'), fname)
#ifndef NOT_STANDALONE
   call check_err(set_attributes(ncid_bfm,ocepoint_id,formula_term='water points'), fname)
   call check_err(set_attributes(ncid_bfm,ocepoint_id,compress='none'), fname)
   call check_err(set_attributes(ncid_bfm,botpoint_id,formula_term='bottom points'), fname)
   call check_err(set_attributes(ncid_bfm,botpoint_id,compress='none'), fname)
   call check_err(set_attributes(ncid_bfm,surfpoint_id,formula_term='surface points'), fname)
   call check_err(set_attributes(ncid_bfm,surfpoint_id,compress='none'), fname)
#endif
#ifdef BFM_GOTM
   call check_err(set_attributes(ncid_bfm,ocepoint_id,formula_term='watercolumn levels'), fname)
   call check_err(set_attributes(ncid_bfm,ocepoint_id,compress='z'), fname)
   call check_err(set_attributes(ncid_bfm,surfpoint_id,formula_term='watercolumn surface'), fname)
   call check_err(set_attributes(ncid_bfm,surfpoint_id,compress='z'), fname)
   call check_err(set_attributes(ncid_bfm,botpoint_id,formula_term='watercolumn bottom'), fname)
   call check_err(set_attributes(ncid_bfm,botpoint_id,compress='z'), fname)
#endif
#ifdef BFM_NEMO
   if (present(column)) then
      call check_err(set_attributes(ncid_bfm,ocepoint_id,formula_term='watercolumn levels'), fname)
      call check_err(set_attributes(ncid_bfm,ocepoint_id,units='meters'), fname)
      call check_err(set_attributes(ncid_bfm,ocepoint_id,long_name='depth'), fname)
      call check_err(set_attributes(ncid_bfm,surfpoint_id,formula_term='watercolumn surface'), fname)
      call check_err(set_attributes(ncid_bfm,surfpoint_id,units='meters'), fname)
      call check_err(set_attributes(ncid_bfm,surfpoint_id,long_name='depth'), fname)
      call check_err(set_attributes(ncid_bfm,botpoint_id,formula_term='watercolumn bottom'), fname)
      call check_err(set_attributes(ncid_bfm,botpoint_id,units='meters'), fname)
      call check_err(set_attributes(ncid_bfm,botpoint_id,long_name='depth'), fname)
   else
      call check_err(set_attributes(ncid_bfm,ocepoint_id,formula_term='water points'), fname)
      call check_err(set_attributes(ncid_bfm,ocepoint_id,compress='x y z'), fname)
      call check_err(set_attributes(ncid_bfm,botpoint_id,formula_term='bottom points'), fname)
      call check_err(set_attributes(ncid_bfm,botpoint_id,compress='x y z'), fname)
      call check_err(set_attributes(ncid_bfm,surfpoint_id,formula_term='surface points'), fname)
      call check_err(set_attributes(ncid_bfm,surfpoint_id,compress='x y z'), fname)
   end if
#endif
   select case (ncdf_time_unit)
      case(0)                           ! seconds
         write(ncdf_time_str,100) 'seconds',trim(start_time)
      case(1)                           ! minutes
         write(ncdf_time_str,100) 'minutes',trim(start_time)
      case(2)                           ! hours
         write(ncdf_time_str,100) 'hours',trim(start_time)
      case default
         write(ncdf_time_str,100) 'seconds',trim(start_time)
   end select
100 format(A,' since ',A)
   call check_err(set_attributes(ncid_bfm,time_id,units=trim(ncdf_time_str)), fname)

   !---------------------------------------------
   !  global attributes
   !---------------------------------------------
   call check_err(NF90_PUT_ATT(ncid_bfm,NF90_GLOBAL,'Title',title), fname)
   history = RELEASE
   call check_err(NF90_PUT_ATT(ncid_bfm,NF90_GLOBAL,'history',history), fname)
   call check_err(NF90_PUT_ATT(ncid_bfm,NF90_GLOBAL,'Conventions','CF-1.0'), fname)
   call check_err(NF90_PUT_ATT(ncid_bfm,NF90_GLOBAL,'Experiment',expinfo), fname)

   !---------------------------------------------
   ! leave define mode
   !---------------------------------------------
   call check_err(NF90_ENDDEF(ncid_bfm), fname)

   !---------------------------------------------
   ! save coordinate variables
   !---------------------------------------------
   if (present(lon)) &
      call check_err(store_data(ncid_bfm,lon_id,POINT,1,scalar=lon), fname)
   if (present(lat)) &
      call check_err(store_data(ncid_bfm,lat_id,POINT,1,scalar=lat), fname) 
   if (present(lat2d)) &
      call check_err(store_data(ncid_bfm,lat_id,XY_SHAPE,NO_BOXES_Z, &
                        array2d=lat2d), fname)
   if (present(lon2d)) &
      call check_err(store_data(ncid_bfm,lon_id,XY_SHAPE,NO_BOXES_Z, &
                        array2d=lon2d), fname)
   if (present(z)) &
      call check_err(store_data(ncid_bfm,depth_id,Z_SHAPE,NO_BOXES_Z,array=z), fname)
   if (present(dz)) &
      call check_err(store_data(ncid_bfm,depth_id,Z_SHAPE,NO_BOXES_Z,array=dz), fname)
   if (present(oceanpoint)) &
      call check_err(store_data(ncid_bfm,ocepoint_id,G_SHAPE,NO_BOXES,iarray=oceanpoint), fname)
   if (present(bottompoint)) &
      call check_err(store_data(ncid_bfm,botpoint_id,G_SHAPE,NO_BOXES_XY,iarray=bottompoint), fname)
   if (present(surfacepoint)) &
      call check_err(store_data(ncid_bfm,surfpoint_id,G_SHAPE,NO_BOXES_XY,iarray=surfacepoint), fname)
   if (present(mask3d)) &
      call check_err(store_data(ncid_bfm,mask_id,XYZ_SHAPE,NO_BOXES_Z, &
                        array3d=mask3d), fname)
   if (present(roceanpoint)) &
      call check_err(store_data(ncid_bfm,ocepoint_id,G_SHAPE,NO_BOXES,array=roceanpoint), fname)
   if (present(rbottompoint)) &
      call check_err(store_data(ncid_bfm,botpoint_id,POINT,NO_BOXES_XY,scalar=rbottompoint), fname)
   if (present(rsurfacepoint)) &
      call check_err(store_data(ncid_bfm,surfpoint_id,POINT,NO_BOXES_XY,scalar=rsurfacepoint), fname)

   !---------------------------------------------
   ! syncronize
   !---------------------------------------------
   call check_err(NF90_SYNC(ncid_bfm), fname)

   ! Flush the log File
   Call FLUSH (LOGUNIT)

   end subroutine init_netcdf_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the netcdf restart
!
! !INTERFACE:
   subroutine init_netcdf_rst_bfm(title,start_time,time_unit,lat,lon,z,dz, &
                              lat2d,lon2d,oceanpoint,surfacepoint,     &
                              bottompoint,mask3d)
!
! !DESCRIPTION:
!  Prepare the netcdf restart file for the BFM
!
! !USES:

   use mem, only: NO_D3_BOX_STATES, NO_BOXES,    &
                  NO_BOXES_XY
#if defined INCLUDE_SEAICE
   use mem, only: NO_D2_BOX_STATES_ICE
#endif
#if defined INCLUDE_BEN
   use mem, only: NO_D2_BOX_STATES_BEN
#endif
   implicit none
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*), intent(in)                    :: title,start_time
   integer, intent(in)                             :: time_unit
   real(RLEN), intent(in),optional                 :: lat,lon
   real(RLEN), intent(in),dimension(:,:),optional  :: lat2d,lon2d
   real(RLEN), intent(in),dimension(:),optional    :: z,dz
   integer, intent(in),dimension(:),optional       :: oceanpoint
   integer, intent(in),dimension(:),optional       :: surfacepoint,bottompoint
   real(RLEN),intent(in),dimension(:,:,:),optional :: mask3d
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications: Marcello Vichi
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: ext,fname
   integer                   :: iret,ndims
   character(len=128)        :: ncdf_time_str,history
!!
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'init_netcdf_rst_bfm: define output restart file(s) ...'

   !---------------------------------------------
   ! Prepare the netcdf file
   !---------------------------------------------
   ext = 'nc'
   fname = './'// TRIM(title) // '.' // ext
   LEVEL2 'Restart NetCDF file is :'
   LEVEL2 TRIM(fname)
   call check_err(NF90_CREATE(fname,NF90_NETCDF4,ncid_rst), fname)

   ncdf_time_unit = time_unit

   !---------------------------------------------
   ! define 3D dimensions and variables
   !---------------------------------------------
   if (present(lon).AND.present(lat)) then
      call check_err(NF90_DEF_DIM(ncid_rst, 'lon', max(NO_BOXES_XY,1), lon_rdim), fname)
      call check_err(NF90_DEF_DIM(ncid_rst, 'lat', max(NO_BOXES_XY,1), lat_rdim), fname)
   else if (present(lon2d).AND.present(lat2d)) then
      call check_err(NF90_DEF_DIM(ncid_rst, 'x', NO_BOXES_X, x_rdim), fname)
      call check_err(NF90_DEF_DIM(ncid_rst, 'y', NO_BOXES_Y, y_rdim), fname)
   else
      stop '### init_netcdf_rst_bfm: lat and lon must be given'
   end if
   call check_err(NF90_DEF_DIM(ncid_rst, 'z', NO_BOXES_Z, depth_rdim), fname)
   call check_err(NF90_DEF_DIM(ncid_rst, 'oceanpoint', max(NO_BOXES,1), ocepoint_rdim), fname)
   call check_err(NF90_DEF_DIM(ncid_rst, 'surfacepoint', max(NO_BOXES_XY,1), surfpoint_rdim), fname)
   call check_err(NF90_DEF_DIM(ncid_rst, 'bottompoint', max(NO_BOXES_XY,1), botpoint_rdim), fname)
   call check_err(NF90_DEF_DIM(ncid_rst, 'time', NF90_UNLIMITED, time_rdim), fname)
   call check_err(NF90_DEF_DIM(ncid_rst, 'char_max', LEN(var_names), chars_rdim), fname)
   call check_err(NF90_DEF_DIM(ncid_rst, 'd3vars', NO_D3_BOX_STATES, d3vars_rdim), fname)
   ALLOCATE(dims(2))
   dims(1) = d3vars_rdim
   dims(2) = ocepoint_rdim
   call check_err(NF90_DEF_VAR(ncid_rst,'D3STATE',NF90_DOUBLE,dims,d3state_rid), fname)
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_rst,d3state_rid,nc_shuffle,nc_deflate,nc_defllev))
   call check_err(NF90_DEF_VAR(ncid_rst,'D3STATE_NAME',NF90_CHAR,(/chars_rdim, d3vars_rdim/),d3state_name_rid), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'D3STATE_UNITS',NF90_CHAR,(/chars_rdim, d3vars_rdim/),d3state_units_rid), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'D3STATE_LONG',NF90_CHAR,(/chars_rdim, d3vars_rdim/),d3state_long_rid), fname)
#ifdef INCLUDE_PELCO2
   call check_err(NF90_DEF_VAR(ncid_rst,'pH',NF90_DOUBLE,ocepoint_rdim,ph_rid), fname)
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_rst,ph_rid,nc_shuffle,nc_deflate,nc_defllev))
#endif
#ifdef BFM_POM
   call check_err(NF90_DEF_VAR(ncid_rst,'D3STATEB',NF90_DOUBLE,dims,d3stateb_rid), fname)
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_rst,d3stateb_rid,nc_shuffle,nc_deflate,nc_defllev))
#endif

   !---------------------------------------------
   ! define coordinate variables
   !---------------------------------------------
   dims(1) = x_rdim
   dims(2) = y_rdim
   if (present(lon)) then
     call check_err(NF90_DEF_VAR(ncid_rst,'lon',NF90_REAL,lon_rdim,lon_rid), fname)
   elseif (present(lon2d)) then
     call check_err(NF90_DEF_VAR(ncid_rst,'lon',NF90_REAL,dims,lon_rid), fname)
   end if
   if (present(lat)) then
     call check_err(NF90_DEF_VAR(ncid_rst,'lat',NF90_REAL,lat_rdim,lat_rid), fname)
   elseif (present(lat2d)) then
     call check_err(NF90_DEF_VAR(ncid_rst,'lat',NF90_REAL,dims,lat_rid), fname)
   end if
   call check_err(NF90_DEF_VAR(ncid_rst,'z',NF90_REAL,depth_rdim,depth_rid), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'oceanpoint', NF90_INT,ocepoint_rdim,ocepoint_rid), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'surfacepoint',NF90_INT,surfpoint_rdim,surfpoint_rid), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'bottompoint',NF90_INT,botpoint_rdim,botpoint_rid), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'time',NF90_REAL,time_rdim,time_rid), fname)

#if defined INCLUDE_SEAICE
   !---------------------------------------------
   ! define 2D Seaice dimensions and variables
   !---------------------------------------------
   call check_err(NF90_DEF_DIM(ncid_rst,'d2vars_ice', NO_D2_BOX_STATES_ICE, d2vars_rdim_ice), fname)
   dims(1) = d2vars_rdim_ice
   dims(2) = surfpoint_rdim
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_ICE',NF90_DOUBLE,dims,d2state_rid_ice), fname)
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_rst,d2state_rid_ice,nc_shuffle,nc_deflate,nc_defllev))
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_ICE_NAME',NF90_CHAR,(/chars_rdim, d2vars_rdim_ice/),d2state_name_rid_ice), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_ICE_UNITS',NF90_CHAR,(/chars_rdim, d2vars_rdim_ice/),d2state_units_rid_ice), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_ICE_LONG',NF90_CHAR,(/chars_rdim, d2vars_rdim_ice/),d2state_long_rid_ice), fname)
#ifdef BFM_POM
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATEB_ICE',NF90_DOUBLE,dims,d2stateb_rid_ice), fname)
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_rst,d2stateb_rid_ice,nc_shuffle,nc_deflate,nc_defllev))
#endif
#endif


#if defined INCLUDE_BEN
   !---------------------------------------------
   ! define 2D Benthic dimensions and variables
   !---------------------------------------------
   call check_err(NF90_DEF_DIM(ncid_rst,'d2vars_ben', NO_D2_BOX_STATES_BEN, d2vars_rdim_ben), fname)
   dims(1) = d2vars_rdim_ben
   dims(2) = botpoint_rdim
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_BEN',NF90_DOUBLE,dims,d2state_rid_ben), fname)
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_rst,d2state_rid_ben,nc_shuffle,nc_deflate,nc_defllev))
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_BEN_NAME',NF90_CHAR,(/chars_rdim, d2vars_rdim_ben/),d2state_name_rid_ben), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_BEN_UNITS',NF90_CHAR,(/chars_rdim, d2vars_rdim_ben/),d2state_units_rid_ben), fname)
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATE_BEN_LONG',NF90_CHAR,(/chars_rdim, d2vars_rdim_ben/),d2state_long_rid_ben), fname)
#ifdef BFM_POM
   call check_err(NF90_DEF_VAR(ncid_rst,'D2STATEB_BEN',NF90_DOUBLE,dims,d2stateb_rid_ben), fname)
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid_rst,d2stateb_rid_ben,nc_shuffle,nc_deflate,nc_defllev))
#endif
#endif

   DEALLOCATE(dims)

   !---------------------------------------------
   ! define mask variables
   !---------------------------------------------
   if (present(mask3d)) then
      ALLOCATE(dims(3))
      dims(1) = x_rdim
      dims(2) = y_rdim
      dims(3) = depth_rdim
      call check_err(NF90_DEF_VAR(ncid_rst,'mask',NF90_REAL,dims,mask_rid), fname)
      DEALLOCATE(dims)
   end if

   !---------------------------------------------
   ! assign attributes
   !---------------------------------------------
   call check_err(set_attributes(ncid_rst,lon_rid,units='degrees_east'), fname)
   call check_err(set_attributes(ncid_rst,lat_rid,units='degrees_north'), fname)
   call check_err(set_attributes(ncid_rst,depth_rid,units='meters'), fname)
#ifndef NOT_STANDALONE
   call check_err(set_attributes(ncid_rst,ocepoint_rid,formula_term='water points'), fname)
   call check_err(set_attributes(ncid_rst,ocepoint_rid,compress='none'), fname)
   call check_err(set_attributes(ncid_rst,surfpoint_rid,formula_term='surface points'), fname)
   call check_err(set_attributes(ncid_rst,surfpoint_rid,compress='none'), fname)
   call check_err(set_attributes(ncid_rst,botpoint_rid,formula_term='bottom points'), fname)
   call check_err(set_attributes(ncid_rst,botpoint_rid,compress='none'), fname)
#endif
#ifdef BFM_GOTM
   call check_err(set_attributes(ncid_rst,ocepoint_rid,formula_term='watercolumn levels'), fname)
   call check_err(set_attributes(ncid_rst,ocepoint_rid,compress='z'), fname)
   call check_err(set_attributes(ncid_rst,surfpoint_rid,formula_term='watercolumn surface'), fname)
   call check_err(set_attributes(ncid_rst,surfpoint_rid,compress='z'), fname)
   call check_err(set_attributes(ncid_rst,botpoint_rid,formula_term='watercolumn bottom'), fname)
   call check_err(set_attributes(ncid_rst,botpoint_rid,compress='z'), fname)
#endif
#ifdef BFM_NEMO
   call check_err(set_attributes(ncid_rst,ocepoint_rid,formula_term='water points'), fname)
   call check_err(set_attributes(ncid_rst,ocepoint_rid,compress='x y z'), fname)
   call check_err(set_attributes(ncid_rst,surfpoint_rid,formula_term='surface points'), fname)
   call check_err(set_attributes(ncid_rst,surfpoint_rid,compress='x y z'), fname)
   call check_err(set_attributes(ncid_rst,botpoint_rid,formula_term='bottom points'), fname)
   call check_err(set_attributes(ncid_rst,botpoint_rid,compress='x y z'), fname)
#endif

   select case (ncdf_time_unit)
      case(0)                           ! seconds
         write(ncdf_time_str,100) 'seconds',trim(start_time)
      case(1)                           ! minutes
         write(ncdf_time_str,100) 'minutes',trim(start_time)
      case(2)                           ! hours
         write(ncdf_time_str,100) 'hours',trim(start_time)
      case default
         write(ncdf_time_str,100) 'seconds',trim(start_time)
   end select
100 format(A,' since ',A)
   call check_err(set_attributes(ncid_rst,time_rid,units=trim(ncdf_time_str)), fname)

   !---------------------------------------------
   !  global attributes
   !---------------------------------------------
   call check_err(NF90_PUT_ATT(ncid_rst,NF90_GLOBAL,'Title',title), fname)
   history = RELEASE
   call check_err(NF90_PUT_ATT(ncid_rst,NF90_GLOBAL,'history',history), fname)
   call check_err(NF90_PUT_ATT(ncid_rst,NF90_GLOBAL,'Conventions','CF-1.0'), fname)

   !---------------------------------------------
   ! leave define mode
   !---------------------------------------------
   call check_err(NF90_ENDDEF(ncid_rst), fname)

   !---------------------------------------------
   ! save coordinate variables
   !---------------------------------------------
   if (present(lon)) &
      call check_err(store_data(ncid_rst,lon_rid,POINT,1,scalar=lon), fname)
   if (present(lat)) &
      call check_err(store_data(ncid_rst,lat_rid,POINT,1,scalar=lat), fname) 
   if (present(lat2d)) &
      call check_err(store_data(ncid_rst,lat_rid,XY_SHAPE,NO_BOXES_Z, &
                        array2d=lat2d), fname)
   if (present(lon2d)) &
      call check_err(store_data(ncid_rst,lon_rid,XY_SHAPE,NO_BOXES_Z, &
                        array2d=lon2d), fname)
   if (present(z)) &
      call check_err(store_data(ncid_rst,depth_rid,Z_SHAPE,NO_BOXES_Z,array=z), fname)
   if (present(dz)) &
      call check_err(store_data(ncid_rst,depth_rid,Z_SHAPE,NO_BOXES_Z,array=dz), fname)
   if (present(oceanpoint)) &
      call check_err(store_data(ncid_rst,ocepoint_rid,G_SHAPE,NO_BOXES,iarray=oceanpoint), fname)
   if (present(bottompoint)) &
      call check_err(store_data(ncid_rst,botpoint_rid,G_SHAPE,NO_BOXES_XY,iarray=bottompoint), fname)
   if (present(surfacepoint)) &
      call check_err(store_data(ncid_rst,surfpoint_rid,G_SHAPE,NO_BOXES_XY,iarray=surfacepoint), fname)
   if (present(mask3d)) &
      call check_err(store_data(ncid_rst,mask_rid,XYZ_SHAPE,NO_BOXES_Z, &
                        array3d=mask3d), fname)

   !---------------------------------------------
   ! syncronize
   !---------------------------------------------
   call check_err(NF90_SYNC(ncid_rst), fname)

   ! Flush the log File
   Call FLUSH (LOGUNIT)

end subroutine init_netcdf_rst_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Store the restart file
!
! !INTERFACE:
  subroutine save_rst_bfm(time)
!
! !DESCRIPTION:
! output restart file of BFM variables 
!
! !USES:
   use mem, only: D3STATE, NO_D3_BOX_STATES, NO_BOXES
#ifdef INCLUDE_PELCO2
   use mem, only: D3DIAGNOS,pppH
#endif
#ifdef BFM_POM
   use api_bfm, only: D3STATEB
#endif

#if defined INCLUDE_SEAICE
   use mem, only: D2STATE_ICE, NO_D2_BOX_STATES_ICE
#ifdef BFM_POM
   use api_bfm, only: D2STATEB_ICE
#endif
#endif

#if defined INCLUDE_BEN
   use mem, only: D2STATE_BEN, NO_D2_BOX_STATES_BEN
#ifdef BFM_POM
   use api_bfm, only: D2STATEB_BEN
#endif
#endif

   implicit none
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in)     :: time
! !LOCAL VARIABLES:
   integer                   :: iret, iter
   character(len=80)         :: restfile
   real(RLEN)                :: temp_time
   character(len=LEN(var_names)), dimension(NO_D3_BOX_STATES) :: tmp_d3names,tmp_d3units,tmp_d3long
#ifdef INCLUDE_SEAICE
   character(len=LEN(var_names)), dimension(NO_D2_BOX_STATES_ICE) :: tmp_d2names_ice,tmp_d2units_ice,tmp_d2long_ice
#endif
#ifdef INCLUDE_BEN
   character(len=LEN(var_names)), dimension(NO_D2_BOX_STATES_BEN) :: tmp_d2names_ben,tmp_d2units_ben,tmp_d2long_ben
#endif
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV) 
!
!EOP
!-----------------------------------------------------------------------
!BOC

!  Storing the time - both the coordinate and later a time string.
     select case (ncdf_time_unit)
        case(0)                           ! seconds
           temp_time = time
        case(1)                           ! minutes
           temp_time = time/60.
        case(2)                           ! hours
           temp_time = time/3600.
        case default
           temp_time = time
     end select
     iret = store_data(ncid_rst,time_rid,POINT,1,scalar=temp_time)


     restfile="out_restart"
     start(1) = 1;   edges(1) = NO_D3_BOX_STATES
     start(2) = 1;   edges(2) = NO_BOXES
     tmp_d3names(:) = var_names(stPelStateS:stPelStateE)
     do iter=1, SIZE(tmp_d3names)
        call replace_char(str=tmp_d3names(iter), Tar='()', rep='_')
     end do
     tmp_d3units(:) = var_units(stPelStateS:stPelStateE)
     tmp_d3long(:)  = var_long(stPelStateS:stPelStateE)
     call check_err(NF90_PUT_VAR(ncid_rst,d3state_rid,D3STATE(:,:),start,edges), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d3state_name_rid, tmp_d3names, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d3names), NO_D3_BOX_STATES /)), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d3state_units_rid, tmp_d3units, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d3units), NO_D3_BOX_STATES /)), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d3state_long_rid, tmp_d3long, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d3long), NO_D3_BOX_STATES /)), restfile)
#ifdef INCLUDE_PELCO2
     call check_err(NF90_PUT_VAR(ncid_rst,ph_rid,D3DIAGNOS(pppH,:),start=(/1/),count=(/NO_BOXES/)), restfile)
#endif
#ifdef BFM_POM
     call check_err(NF90_PUT_VAR(ncid_rst,d3stateb_rid,D3STATEB(:,:),start,edges), restfile)
#endif

#if defined INCLUDE_SEAICE
     start(1) = 1;   edges(1) = NO_D2_BOX_STATES_ICE
     start(2) = 1;   edges(2) = NO_BOXES_XY
     tmp_d2names_ice(:) = var_names(stIceStateS:stIceStateE)
     do iter=1, SIZE(tmp_d2names_ice)
        call replace_char(str=tmp_d2names_ice(iter), tar='()', rep='_')
     end do
     tmp_d2units_ice(:) = var_units(stIceStateS:stIceStateE)
     tmp_d2long_ice(:)  = var_long(stIceStateS:stIceStateE)
     call check_err(NF90_PUT_VAR(ncid_rst,d2state_rid_ice,D2STATE_ICE(:,:),start,edges), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d2state_name_rid_ice, tmp_d2names_ice, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d2names_ice), NO_D2_BOX_STATES_ICE /)), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d2state_units_rid_ice, tmp_d2units_ice, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d2units_ice), NO_D2_BOX_STATES_ICE /)), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d2state_long_rid_ice, tmp_d2long_ice, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d2long_ice), NO_D2_BOX_STATES_ICE /)), restfile)
#ifdef BFM_POM
     call check_err(NF90_PUT_VAR(ncid_rst,d2state_rid_ice,D2STATEB_ICE(:,:),start,edges), restfile)
#endif
#endif

#if defined INCLUDE_BEN
     start(1) = 1;   edges(1) = NO_D2_BOX_STATES_BEN
     start(2) = 1;   edges(2) = NO_BOXES_XY
     tmp_d2names_ben(:) = var_names(stBenStateS:stBenStateE)
     do iter=1, SIZE(tmp_d2names_ben)
        call replace_char(str=tmp_d2names_ben(iter), tar='()', rep='_')
     end do
     tmp_d2units_ben(:) = var_units(stBenStateS:stBenStateE)
     tmp_d2long_ben(:)  = var_long(stBenStateS:stBenStateE)
     call check_err(NF90_PUT_VAR(ncid_rst,d2state_rid_ben,D2STATE_BEN(:,:),start,edges), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d2state_name_rid_ben, tmp_d2names_ben, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d2names_ben), NO_D2_BOX_STATES_BEN /)), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d2state_units_rid_ben, tmp_d2units_ben, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d2units_ben), NO_D2_BOX_STATES_BEN /)), restfile)
     call check_err(NF90_PUT_VAR( ncid_rst, d2state_long_rid_ben, tmp_d2long_ben, &
          start=(/ 1, 1 /), count=(/ LEN(tmp_d2long_ben), NO_D2_BOX_STATES_BEN /)), restfile)
#ifdef BFM_POM
     call check_err(NF90_PUT_VAR(ncid_rst,d2state_rid_ben,D2STATEB_BEN(:,:),start,edges), restfile)
#endif
#endif
     LEVEL2 'save_rst_bfm: Restart data has been written'
! the file is closed in the main (in case of more restart files)

  end subroutine save_rst_bfm 
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read the restart file
!
! !INTERFACE:
  subroutine read_rst_bfm(title)
!
! !DESCRIPTION:
! Read restart file of BFM variables 
!
! !USES:
   use mem, only: D3STATE, NO_D3_BOX_STATES, NO_BOXES
#ifdef INCLUDE_PELCO2
   use mem, only: D3DIAGNOS,pppH
#endif
#ifdef BFM_POM
   use api_bfm, only: D3STATEB
#endif

#if defined INCLUDE_SEAICE
   use mem, only: D2STATE_ICE, NO_D2_BOX_STATES_ICE
#ifdef BFM_POM
   use api_bfm, only: D2STATEB_ICE
#endif
#endif

#if defined INCLUDE_BEN
   use mem, only: D2STATE_BEN, NO_D2_BOX_STATES_BEN
#ifdef BFM_POM
   use api_bfm, only: D2STATEB_BEN
#endif
#endif
   implicit none
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)                 :: title
   character(len=NF90_MAX_NAME)                :: namedimt 
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: ext,fname
   integer                   :: iret
   integer                   :: nstate_id,nstate_len
   integer                   :: ncomp_id,ncomp_len
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV) 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'read_rst_bfm: READ initial conditions from multiple restart files ...'

   !---------------------------------------------
   ! open the netcdf restart file
   !---------------------------------------------
   ext = 'nc'
   fname = './'// TRIM(title) // '.' // ext
   LEVEL2 'Restart NetCDF file is:'
   LEVEL2 TRIM(fname)
   call check_err(NF90_OPEN(fname,NF90_NOWRITE,ncid_rst_in), fname)
   !---------------------------------------------
   ! Check 3D dimensions 
   !---------------------------------------------
   call check_err(NF90_INQ_DIMID(ncid_rst_in,"d3vars",nstate_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_in,nstate_id,namedimt,nstate_len), fname)
   call check_err(NF90_INQ_DIMID(ncid_rst_in,"oceanpoint",ncomp_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_in,ncomp_id,namedimt,ncomp_len), fname)
   if (nstate_len/=NO_D3_BOX_STATES .OR. ncomp_len/=NO_BOXES .AND. ncomp_len > 1) then
      LEVEL1 "3D Dimension mismatch in restart file:"
      LEVEL2 TRIM(fname)
      LEVEL3 "NO_D3_BOX_STATES in model:", NO_D3_BOX_STATES
      LEVEL3 "NO_D3_BOX_STATES in file:", nstate_len
      LEVEL3 "NO_BOXES in model:", NO_BOXES
      LEVEL3 "NO_BOXES in file:", ncomp_len
      stop "STOP in read_rst_bfm contained in netcdf_bfm.F90"
   end if
   !---------------------------------------------
   ! Initialize 3D variable
   !---------------------------------------------
   call check_err(NF90_INQ_VARID(ncid_rst_in,"D3STATE",nstate_id), fname)
   call check_err(NF90_GET_VAR(ncid_rst_in,nstate_id,D3STATE(:,:)), fname)
#ifdef INCLUDE_PELCO2
   call check_err(NF90_INQ_VARID(ncid_rst_in,"pH",nstate_id), fname)
   call check_err(NF90_GET_VAR(ncid_rst_in,nstate_id,D3DIAGNOS(pppH,:)), fname)
#endif 
#ifdef BFM_POM
   call check_err(NF90_INQ_VARID(ncid_rst_in,"D3STATEB",nstate_id), fname)
   call check_err(NF90_GET_VAR(ncid_rst_in,nstate_id,D3STATEB(:,:)), fname)
#endif

#if defined INCLUDE_SEAICE
   !---------------------------------------------
   ! Check Seaice 2D dimensions 
   !---------------------------------------------
   call check_err(NF90_INQ_DIMID(ncid_rst_in,"d2vars_ice",nstate_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_in,nstate_id,namedimt,nstate_len), fname)
   call check_err(NF90_INQ_DIMID(ncid_rst_in,"surfacepoint",ncomp_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_in,ncomp_id,namedimt,ncomp_len), fname)
   if (nstate_len/=NO_D2_BOX_STATES_ICE .OR. ncomp_len/=NO_BOXES_XY) then
      LEVEL1 '2D Seaice Dimension mismatch in restart file:'
      LEVEL2 TRIM(fname)
      LEVEL3 "NO_D2_BOX_STATES_ICE in model:", NO_D2_BOX_STATES_ICE
      LEVEL3 "NO_D2_BOX_STATES_ICE in file:", nstate_len
      LEVEL3 "NO_BOXES_XY in model:", NO_BOXES_XY
      LEVEL3 "NO_BOXES_XY in file:", ncomp_len
      stop 'STOP in read_rst_bfm contained in netcdf_bfm.F90'
   end if
   !---------------------------------------------
   ! Initialize 2D variable
   !---------------------------------------------
   call check_err(NF90_INQ_VARID(ncid_rst_in,"D2STATE_ICE",nstate_id), fname)
   call check_err(NF90_GET_VAR(ncid_rst_in,nstate_id,D2STATE_ICE(:,:)), fname)
#ifdef BFM_POM
   call check_err(NF90_INQ_VARID(ncid_rst_in,"D2STATEB_ICE",nstate_id), fname)
   call check_err(NF90_GET_VAR(ncid_rst_in,nstate_id,D2STATEB_ICE(:,:)), fname)
#endif
#endif


#if defined INCLUDE_BEN
   !---------------------------------------------
   ! Check Benthic 2D dimensions 
   !---------------------------------------------
   call check_err(NF90_INQ_DIMID(ncid_rst_in,"d2vars_ben",nstate_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_in,nstate_id,namedimt,nstate_len), fname)
   call check_err(NF90_INQ_DIMID(ncid_rst_in,"bottompoint",ncomp_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_in,ncomp_id,namedimt,ncomp_len), fname)
   if (nstate_len/=NO_D2_BOX_STATES_BEN .OR. ncomp_len/=NO_BOXES_XY) then
      LEVEL1 '2D Benthic Dimension mismatch in restart file:'
      LEVEL2 TRIM(fname)
      LEVEL3 "NO_D2_BOX_STATES_BEN in model:", NO_D2_BOX_STATES_BEN
      LEVEL3 "NO_D2_BOX_STATES_BEN in file:", nstate_len
      LEVEL3 "NO_BOXES_XY in model:", NO_BOXES_XY
      LEVEL3 "NO_BOXES_XY in file:", ncomp_len
      stop 'STOP in read_rst_bfm contained in netcdf_bfm.F90'
   end if
   !---------------------------------------------
   ! Initialize 2D variable
   !---------------------------------------------
   call check_err(NF90_INQ_VARID(ncid_rst_in,"D2STATE_BEN",nstate_id), fname)
   call check_err(NF90_GET_VAR(ncid_rst_in,nstate_id,D2STATE_BEN(:,:)), fname)
#ifdef BFM_POM
   call check_err(NF90_INQ_VARID(ncid_rst_in,"D2STATEB_BEN",nstate_id), fname)
   call check_err(NF90_GET_VAR(ncid_rst_in,nstate_id,D2STATEB_BEN(:,:)), fname)
#endif
#endif

   LEVEL2 'read_rst_bfm: READ multiple restart files ... DONE! '
   LEVEL2 ' '
   call check_err(NF90_CLOSE(ncid_rst_in), fname)

  end subroutine read_rst_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read the restart file
!
! !INTERFACE:
  subroutine read_rst_bfm_glo(title, narea, jpnij, &
        jpiglo, jpjglo, jpkglo, &
        nlcit, nlcjt, &
        nldit, nldjt, &
        nleit, nlejt, &
        nimppt, njmppt, &
        SEAmask )
!
! !DESCRIPTION:
! Read restart file of BFM variables from one merged file in 3d 
!
! !USES:
   use mem, only: D3STATE
#ifdef INCLUDE_PELCO2
   use mem, only: D3DIAGNOS,pppH
#endif

   implicit none
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)          :: title                  ! name of the file to open
   integer, intent(in)                   :: narea                  ! index of subdomain
   integer, intent(in)                   :: jpnij                  ! nb of local domain = nb of processors ( <= jpni x jpnj )
   integer, intent(in)                   :: jpiglo, jpjglo, jpkglo ! i-, j-, k-dimensions of the global domain
   integer, intent(in), dimension(:)     :: nlcit,  nlcjt          !: dimensions of every subdomain
   integer, intent(in), dimension(:)     :: nldit,  nldjt          !: first, last indoor index for each i-domain
   integer, intent(in), dimension(:)     :: nleit,  nlejt          !: first, last indoor index for each j-domain
   integer, intent(in), dimension(:)     :: nimppt, njmppt         !: i-, j-indexes for each processor
   logical, intent(in), dimension(:,:,:) :: SEAmask !: 3D boolean Land-sea mask
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX) :: fname
   integer                 :: ncid_rst_3d, ncomp_id, ncomp_len, iret

   integer :: idx_var, idx_var_array, vid

   integer,dimension(4)                      :: array_3d_start, array_3d_count, array_3d_end
   real(RLEN),allocatable,dimension(:,:,:,:) :: array_3d


   integer :: iniI, iniJ, cntI, cntJ, cntK

   integer :: idx_i, idx_j, idx_k, noce

   character(len=PATH_MAX) :: fname_ph
   integer :: ncid_ph, IDx, IDy, IDz, IDtime, IDboxes, IDtarget, IDtarget_mask, IDtarget_box

   character(len=NF90_MAX_NAME) :: string
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV) 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'read_rst_bfm_glo: READ initial conditions from single restart file ...'

   !---------------------------------------------
   ! open the netcdf restart file
   !---------------------------------------------
   fname = TRIM(title) // '.' // 'nc'
   LEVEL2 'Reading Restart file in NetCDF 3D'
   LEVEL2 TRIM(fname)
   call check_err(NF90_OPEN(fname,NF90_NOWRITE,ncid_rst_3d), fname)
   !---------------------------------------------
   ! Check 3D dimensions 
   !---------------------------------------------
   call check_err(NF90_INQ_DIMID(ncid_rst_3d,"x",ncomp_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_3d,ncomp_id, len=ncomp_len), fname)
   if (ncomp_len/=jpiglo) then
      LEVEL1 "Global X Dimension mismatch in 3d restart file:"
      LEVEL2 TRIM(fname)
      LEVEL3 "DIM X in model:", jpiglo
      LEVEL3 "DIM X in file:",  ncomp_len
      stop "STOP in read_rst_bfm_glo contained in netcdf_bfm.F90"
   end if
   call check_err(NF90_INQ_DIMID(ncid_rst_3d,"y",ncomp_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_3d,ncomp_id, len=ncomp_len), fname)
   if (ncomp_len/=jpjglo) then
      LEVEL1 "Global Y Dimension mismatch in 3d restart file:"
      LEVEL2 TRIM(fname)
      LEVEL3 "DIM Y in model:", jpjglo
      LEVEL3 "DIM Y in file:",  ncomp_len
      stop "STOP in read_rst_bfm_glo contained in netcdf_bfm.F90"
   end if
   call check_err(NF90_INQ_DIMID(ncid_rst_3d,"depth",ncomp_id), fname)
   call check_err(NF90_INQUIRE_DIMENSION(ncid_rst_3d,ncomp_id, len=ncomp_len), fname)
   if (ncomp_len/=jpkglo) then
      LEVEL1 "Global Z Dimension mismatch in 3d restart file:"
      LEVEL2 TRIM(fname)
      LEVEL3 "DIM Z in model:", jpkglo
      LEVEL3 "DIM Z in file:",  ncomp_len
      stop "STOP in read_rst_bfm_glo contained in netcdf_bfm.F90"
   end if

   !---------------------------------------------
   ! get the coordinates of the sub-domain
   !---------------------------------------------
   iniI = nimppt(narea)
   iniJ = njmppt(narea)
   cntI = nlcit(narea)
   cntJ = nlcjt(narea)
   cntK = jpkglo


   !allocate subdomain array and global mask
   allocate( array_3d(cntI,cntJ,cntK,1) )
   array_3d_start = (/ iniI       , iniJ       , 1   , 1 /)
   array_3d_count = (/ cntI       , cntJ       , cntk, 1 /)
   array_3d_end   = (/ iniI+cntI-1, iniJ+cntJ-1, cntk, 1 /)

   !---------------------------------------------
   ! Initialize 3D Pelagic variables
   !---------------------------------------------
   do idx_var=stPelStateS, stPelStateE
      idx_var_array = idx_var - stPelStateS + 1
      string = var_names(idx_var)
      iret = NF90_INQ_VARID(ncid_rst_3d, string, vid)
      if( iret /= NF90_NOERR ) then
         ! in new BFM version netcdf output remove '(' and ')' 
         call replace_char(str=string, tar='()', rep='_')
         iret = NF90_INQ_VARID(ncid_rst_3d, string, vid)
      end if
      if( iret == NF90_NOERR ) then
         call check_err(nf90_get_var(ncid_rst_3d, vid, array_3d, &
              start=array_3d_start, count=array_3d_count), fname)

         noce = 0
         do idx_k=1,cntK
            do idx_j=1,cntJ
               do idx_i=1,cntI
                  if( SEAmask(idx_i,idx_j,idx_k) ) then
                     noce = noce+1
                     D3STATE(idx_var_array,noce) = array_3d(idx_i,idx_j,idx_k,1)
                  end if
               end do
            end do
         end do

      end if
   end do

   !---------------------------------------------
   ! Initialize Ph variable
   !---------------------------------------------
#ifdef INCLUDE_PELCO2
   call check_err(NF90_INQ_VARID(ncid_rst_3d,"pH", vid), fname)
   call check_err(nf90_get_var(ncid_rst_3d, vid, array_3d, &
        start=array_3d_start, count=array_3d_count), fname)

   noce = 0
   do idx_k=1,cntK
      do idx_j=1,cntJ
         do idx_i=1,cntI
            if( SEAmask(idx_i,idx_j,idx_k) ) then
               noce = noce+1
               D3DIAGNOS(pppH,noce) = array_3d(idx_i,idx_j,idx_k,1)
            end if
         end do
      end do
   end do

#endif
   !---------------------------------------------
   ! Initialize 2D Benthic variables
   !---------------------------------------------
#if defined INCLUDE_BEN
      LEVEL2 "read_rst_bfm_glo: READ 2D benthic initial conditions"
      stop "STOP in read_rst_bfm_glo: read benthic initial conditions not yet implemented"
#endif

   call check_err(nf90_close(ncid_rst_3d),fname)

   LEVEL2 'read_rst_bfm_glo: READ single restart file ... DONE! '
   LEVEL2 ' '

   if(allocated(array_3d)) deallocate(array_3d)

  end subroutine read_rst_bfm_glo
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Intialise the storage of results in NetCDF
!
! !INTERFACE:
   subroutine init_save_bfm()
!
! !DESCRIPTION:
! Preparation of the netcdf output.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:

!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  Adapted to BFM: Marcello Vichi (INGV) & Piet Ruardij (NIOZ)
!
! !LOCAL VARIABLES:
   integer, save             :: nn       ! number pel.var to be saved 
   integer, save             :: nnb      ! number ben.var to be saved 
   integer                   :: iret,rc
   real(RLEN)                :: ltime
   integer                   :: out_unit=67
   integer                   :: i,j,n
!EOP
!-----------------------------------------------------------------------
!BOC

   !---------------------------------------------
   ! Enter define mode
   !---------------------------------------------
   iret = define_mode(ncid_bfm,.true.)

      ALLOCATE(dims(2))
      dims(1) = ocepoint_dim
      dims(2) = time_dim

      do n=stPelStateS,stPelFluxE
         if ( var_ids(n) /= 0 )  then
            iret = new_nc_variable(ncid_bfm,var_names(n),NF90_REAL, &
                 dims,var_ids(n))
            iret = set_attributes(ncid_bfm,var_ids(n),            &
                 units=var_units(n),         &
                 long_name=var_long(n))
         end if 
      end do

      do n=stPelDiag2dS,stPelRivE
         dims(1) = botpoint_dim 
         if ( n >= stPelDiag2dS .AND. n <= stPelSurE ) dims(1) = surfpoint_dim

         if ( var_ids(n) /= 0 )  then 
            iret = new_nc_variable(ncid_bfm,var_names(n),NF90_REAL, &
                 dims,var_ids(n))
            iret = set_attributes(ncid_bfm,var_ids(n),            &
                 units=var_units(n),         &
                 long_name=var_long(n)) 
         endif
      end do

#if defined INCLUDE_SEAICE
      do n=stIceStart,stIceEnd

         dims(1) = surfpoint_dim

         if ( var_ids(n) /= 0 )  then 
            iret = new_nc_variable(ncid_bfm,var_names(n),NF90_REAL, &
                 dims,var_ids(n))
            iret = set_attributes(ncid_bfm,var_ids(n),            &
                 units=var_units(n),         &
                 long_name=var_long(n)) 
         endif
      end do
#endif

#if defined INCLUDE_BEN
      do n=stBenStart,stBenEnd
         dims(1) = botpoint_dim

         j=0
#ifdef INCLUDE_BENPROFILES
         ! this is a special part for the variable with
         ! alternative dimensions (eg benthic profiles)
         if ( var_ids(n) /= 0 )  &
            j=special_dims(2,ncid_bfm,NO_BOXES_Z,var_names(n),var_long(n),  &
                           var_units(n),time_dim,var_ids(n))
#endif
         if ( j.eq.0 .and. var_ids(n) /= 0 )  then 
            iret = new_nc_variable(ncid_bfm,var_names(n),NF90_REAL, &
                 dims,var_ids(n))
            iret = set_attributes(ncid_bfm,var_ids(n),            &
                 units=var_units(n),         &
                 long_name=var_long(n)) 
         endif
      end do
#endif


   DEALLOCATE(dims)
   iret = define_mode(ncid_bfm,.false.)
   LEVEL1 'init_save_bfm: output data file(s) creation ... DONE!'
   LEVEL1 ' '

   !---------------------------------------------
   ! Store the initial conditions
   !---------------------------------------------
   if (bfm_rstctl) then
      ltime=0.0
      call save_bfm(ltime)
      LEVEL1 'init_save_bfm: BFM data at initial step saved into the output file.'
      LEVEL1 ' '
   endif
   return
   end subroutine init_save_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Store the results
!
! !INTERFACE:
  subroutine save_bfm(time)
!
! !DESCRIPTION:
! output of BFM variables 
!
! !USES:
   use mem, only: D3STATE,D3DIAGNOS,D3FLUX_FUNC
   use mem, only: D2DIAGNOS
#if defined INCLUDE_SEAICE
   use mem, only: D2STATE_ICE,D2DIAGNOS_ICE,D2DIAGNOS_ICE,D2FLUX_FUNC_ICE
#endif
#if defined INCLUDE_BEN
   use mem, only: D2STATE_BEN,D2DIAGNOS_BEN,D2DIAGNOS_BEN,D2FLUX_FUNC_BEN
#endif
   implicit none
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in)     :: time
! !LOCAL VARIABLES:
   integer                   :: iret
   integer                   :: k,n,idx_tmp
   real(RLEN)                :: temp_time
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  Adapted to BFM: Marcello Vichi (INGV) & Piet Ruardij (NIOZ)
!  Rev. 2012 : Tomas Lovato (CMCC)
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef BFM_STANDALONE
   LEVEL1 'save_bfm: SAVE bfm output data at day ',time/SEC_PER_DAY
#endif
! increase the time record number
   recnum = recnum + 1

!  Storing the time - both the coordinate and later a time string.
   select case (ncdf_time_unit)
      case(0)                           ! seconds
         temp_time = time
      case(1)                           ! minutes
         temp_time = time/60.
      case(2)                           ! hours
         temp_time = time/3600.
      case default
         temp_time = time
   end select
   iret = store_data(ncid_bfm,time_id,T_SHAPE,1,scalar=temp_time)

   !---------------------------------------------
   ! Pelagic 3D variables
   !---------------------------------------------
   k = 0
   do n = stPelStart , stPelFluxE
      if ( var_ids(n) > 0 ) then
        IF ( .not. var_ave(n) ) THEN
         !-- Store snapshot of pelagic state variables
         if ( n >= stPelStateS .AND. n <= stPelStateE ) then
            idx_tmp=n-stPelStateS+1
            iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=D3STATE(idx_tmp,:))     
         end if
         !-- Store snapshot of pelagic diagnostics
         if ( n >= stPelDiagS .AND. n <= stPelDiagE ) then
            idx_tmp=n-stPelDiagS+1
            iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=D3DIAGNOS(idx_tmp,:))
         end if
         !-- Store snapshot of pelagic fluxes
         if ( n >= stPelFluxS .AND. n <= stPelFluxE ) then 
            idx_tmp=n-stPelFluxS+1
            call correct_flux_output(1,idx_tmp,1,NO_BOXES,c1dim)
            iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=c1dim)
         endif

         ! Store mean values of (any) 3D entity
        ELSE
           if (temp_time /= ZERO ) then
              k=k+1
              iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=D3ave(k,:))
           endif
        ENDIF
 
      endif
   enddo

   !---------------------------------------------
   ! Pelagic 2D variables
   !---------------------------------------------
   k=0
   do n = stPelDiag2dS , stPelRivE
      if ( var_ids(n) > 0 ) then   
         IF ( .not. var_ave(n) ) THEN
            ! Store snapshot of pelagic 2D diagnostics
            if ( n >= stPelDiag2dS .AND. n <= stPelDiag2dE ) then
               idx_tmp=n-stPelDiag2dS+1
               iret = store_data(ncid_bfm,var_ids(n),SURFT_SHAPE,NO_BOXES_XY,garray=D2DIAGNOS(idx_tmp,:))
            end if
            ! Store snapshot of pelagic 2D diagnostics at surface
            if ( n >= stPelSurS .AND. n <= stPelSurE) then
               idx_tmp=n-stPelDiag2dS+1
               iret = store_data(ncid_bfm,var_ids(n),SURFT_SHAPE,NO_BOXES_XY,garray=D2DIAGNOS(idx_tmp,:))
            end if
            ! Store snapshot of pelagic 2D diagnostics at bottom
            if ( n >= stPelBotS .AND. n <= stPelRivE) then
               idx_tmp=n-stPelDiag2dS+1
               iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2DIAGNOS(idx_tmp,:))
            end if
         ELSE
            ! Store mean values of (pelagic) 2D entity
            if ( temp_time /= ZERO ) then
               k=k+1
               iret = store_data(ncid_bfm,var_ids(n),SURFT_SHAPE,NO_BOXES_XY,garray=D2ave(k,:))
            end if
         ENDIF
      end if
   end do

#if defined INCLUDE_SEAICE
   !---------------------------------------------
   ! 2D Seaice variables
   !---------------------------------------------
   k=0
   do n = stIceStart , stIceEnd
      if ( var_ids(n) > 0 ) then   
         IF ( .not. var_ave(n) ) THEN

            ! Store snapshot of seaice 2D state
            if ( n >= stIceStateS .AND. n <= stIceStateE) then
               idx_tmp=n-stIceStateS+1
               iret = store_data(ncid_bfm,var_ids(n),SURFT_SHAPE,NO_BOXES_XY,garray=D2STATE_ICE(idx_tmp,:))
            end if
            ! Store snapshot of seaice 2D diagnostics
            if ( n >= stIceDiag2dS .AND. n <= stIceDiag2dE ) then
               idx_tmp=n-stIceDiag2dS+1
               iret = store_data(ncid_bfm,var_ids(n),SURFT_SHAPE,NO_BOXES_XY,garray=D2DIAGNOS_ICE(idx_tmp,:))
            end if
            ! Store snapshot of seaice 2D flux
            if ( n >= stIceFlux2dS .AND. n <= stIceFlux2dE ) then
               idx_tmp=n-stIceFlux2dS+1
               iret = store_data(ncid_bfm,var_ids(n),SURFT_SHAPE,NO_BOXES_XY,garray=D2FLUX_FUNC_ICE(idx_tmp))
            end if
         ELSE
            ! Store mean values of (any) 2D entity
            if ( temp_time /= ZERO ) then
               k=k+1
               iret = store_data(ncid_bfm,var_ids(n),SURFT_SHAPE,NO_BOXES_XY,garray=D2ave_ice(k,:))
            end if
         ENDIF
      end if
   end do
#endif

#if defined INCLUDE_BEN
   !---------------------------------------------
   ! 2D Benthic variables
   !---------------------------------------------
   k=0
   do n = stBenStart , stBenEnd
      if ( var_ids(n) > 0 ) then   
         IF ( .not. var_ave(n) ) THEN

            ! Store snapshot of benthic 2D state
            if ( n >= stBenStateS .AND. n <= stBenStateE) then
               idx_tmp=n-stBenStateS+1
               iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2STATE_BEN(idx_tmp,:))
            end if
            ! Store snapshot of benthic 2D diagnostics
            if ( n >= stBenDiag2dS .AND. n <= stBenDiag2dE ) then
               idx_tmp=n-stBenDiag2dS+1
               iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2DIAGNOS_BEN(idx_tmp,:))
            end if
            ! Store snapshot of benthic 2D flux
            if ( n >= stBenFlux2dS .AND. n <= stBenFlux2dE ) then
               idx_tmp=n-stBenFlux2dS+1
               iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2FLUX_FUNC_BEN(idx_tmp))
            end if
         ELSE
            ! Store mean values of (any) 2D entity
            if ( temp_time /= ZERO ) then
               k=k+1
               iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2ave_ben(k,:))
            end if
         ENDIF
      end if
   end do
#endif


   iret = NF90_SYNC(ncid_bfm)
   call check_err(iret, 'Save_bfm: writing output')

   LEVEL1 'save_bfm: SAVE bfm output ... DONE! '
   LEVEL1 ' '

   ! Flush the log File
   Call FLUSH (LOGUNIT)

   return
   end subroutine save_bfm
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Close files used for saving model results
!
! !INTERFACE:
   subroutine close_ncdf(ncid)
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Closes the NetCDF file.
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: ncid
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret
!
!-------------------------------------------------------------------------
!BOC

   iret = NF90_CLOSE(ncid)
   call check_err(iret, 'close_ncdf')

   return
   end subroutine close_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Begin or end define mode
!
! !INTERFACE:
   integer function define_mode(ncid,action)
!
! !DESCRIPTION:
!  Depending on the value of the argument {\tt action},
!  this routine put NetCDF in the `define' mode or not.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: ncid
   logical, intent(in)       :: action
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer         :: iret
!
!-----------------------------------------------------------------------
!BOC
   if(action) then
      iret = NF90_REDEF(ncid)
   else
      iret = NF90_ENDDEF(ncid)
   end if
   define_mode = 0
   return
   end function define_mode
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define a new NetCDF variable
!
! !INTERFACE:
   integer function new_nc_variable(ncid,name,data_type,dimids,id)
!
! !DESCRIPTION:
!  This routine is used to define a new variable to store in a NetCDF file.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid
   character(len=*), intent(in)        :: name
   integer, intent(in)                 :: data_type
   integer, intent(in)                 :: dimids(:)
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: id
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  TOM
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret
   character(LEN=LEN(name))  :: string
!
!-----------------------------------------------------------------------
!BOC
   !Replace the name if have "(" or ")"
   string = name
   call replace_char(str=string, tar='()', rep='_')
   iret = NF90_DEF_VAR(ncid,string,data_type,dimids,id)
   call check_err(iret, ('caller: new_nc_variable with input '//trim(string)))
   if (nc_compres) call check_err(NF90_DEF_VAR_DEFLATE(ncid,id,nc_shuffle,nc_deflate,nc_defllev))
   new_nc_variable = iret
   return
   end function new_nc_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set attributes for a NetCDF variable.
!
! !INTERFACE:
   integer function set_attributes(ncid,id,                         &
                                   units,long_name,                 &
                                   valid_min,valid_max,valid_range, &
                                   scale_factor,add_offset,         &
                                   FillValue,missing_value,         &
                                   C_format,FORTRAN_format,         &
                                   compress,formula_term)
!
! !DESCRIPTION:
!  This routine is used to set a number of attributes for
!  variables. The routine makes heavy use of the {\tt optional} keyword.
!  The list of recognized keywords is very easy to extend. 
!  The CF-1.0 convention is used.
!
! !USES:
!  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid,id
   character(len=*), optional          :: units,long_name
   real(RLEN), optional                  :: valid_min,valid_max
   real(RLEN), optional                  :: valid_range(2)
   real(RLEN), optional                  :: scale_factor,add_offset
   real(RLEN), optional                  :: FillValue,missing_value
   character(len=*), optional          :: C_format,FORTRAN_format
   character(len=*), optional          :: compress,formula_term
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!
! !LOCAL VARIABLES:
   integer                   :: len,iret
   REAL_4B                   :: vals(2)
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if(present(units)) then
!      len = len_trim(units)
      iret = NF90_PUT_ATT(ncid,id,'units',units)
   end if

   if(present(long_name)) then
      iret = NF90_PUT_ATT(ncid,id,'long_name',long_name)
   end if

   if(present(C_format)) then
      iret = NF90_PUT_ATT(ncid,id,'C_format',C_format)
   end if

   if(present(FORTRAN_format)) then
      iret = NF90_PUT_ATT(ncid,id,'FORTRAN_format',FORTRAN_format)
   end if

   if(present(compress)) then
      iret = NF90_PUT_ATT(ncid,id,'compress',compress)
   end if

   if(present(formula_term)) then
      iret = NF90_PUT_ATT(ncid,id,'formula_term',formula_term)
   end if

   if(present(valid_min)) then
      iret = NF90_PUT_ATT(ncid,id,'valid_min',valid_min)
   end if

   if(present(valid_max)) then
      iret = NF90_PUT_ATT(ncid,id,'valid_max',valid_max)
   end if

   if(present(valid_range)) then
      vals(1) = valid_range(1)
      vals(2) = valid_range(2)
      iret = NF90_PUT_ATT(ncid,id,'valid_range',vals)
   end if

   if(present(scale_factor)) then
      iret = NF90_PUT_ATT(ncid,id,'scale_factor',scale_factor)
   end if

   if(present(add_offset)) then
      iret = NF90_PUT_ATT(ncid,id,'add_offset',add_offset)
   end if

   if(present(FillValue)) then
      iret = NF90_PUT_ATT(ncid,id,'_FillValue',FillValue)
   end if

   if(present(missing_value)) then
      iret = NF90_PUT_ATT(ncid,id,'missing_value',missing_value)
   end if

   set_attributes = 0
   return
   end function set_attributes
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Store values in a NetCDF file
!
! !INTERFACE:
   integer function store_data(ncid,id,var_shape,nbox,             &
                               iscalar,iarray,scalar,array,garray, &
                               array2d,array3d)
!
! !DESCRIPTION:
!  This routine is used to store a variable in the NetCDF file.
!  The subroutine uses {\tt optional} parameters to find out which data
!  type to save.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid,id,var_shape,nbox
   integer, optional                   :: iscalar
   integer, optional                   :: iarray(1:nbox)
   real(RLEN), optional                  :: scalar
   real(RLEN), optional                  :: array(1:nbox)
   real(RLEN), optional                  :: garray(1:nbox)
   real(RLEN), optional                  :: array2d(:,:)
   real(RLEN), optional                  :: array3d(:,:,:)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications: Marcello Vichi
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret,n=0
   integer                   :: idum(1:nbox)
   REAL_4B                   :: r4,dum(1:nbox)
!
!-----------------------------------------------------------------------
!BOC
   if (.not. present(iscalar) .and. .not. present(iarray)  .and. &
       .not. present(scalar)  .and. .not. present(array)   .and. &
       .not. present(garray)  .and. .not. present(array2d) .and. &
       .not. present(array3d)) then
      FATAL 'At least one optional argument has to be passed to - store_data()'
      stop 'store_data'
   end if
   n = 0
   if(present(iscalar)) n = n+1
   if(present(iarray))  n = n+1
   if(present(scalar))  n = n+1
   if(present(array))   n = n+1
   if(present(garray))  n = n+1
   if(present(array2d)) n = n+1
   if(present(array3d)) n = n+1
   if(n .ne. 1) then
      FATAL 'Only one optional argument must be passed to - store_data()'
      stop 'store_data'
   end if

   if (present(iscalar)) then
      select case (var_shape)
         case(POINT)
            iret = NF90_PUT_VAR(ncid,id,iscalar)
         case(T_SHAPE)
            start(1) = recnum; edges(1) = 1
            idum(1)=iscalar
            iret = NF90_PUT_VAR(ncid,id,idum,start,edges)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(iarray)) then
      select case (var_shape)
         case(Z_SHAPE,G_SHAPE)
            start(1) = 1;   edges(1) = nbox
            idum(1:nbox)=iarray(1:nbox)
            iret = NF90_PUT_VAR(ncid,id,idum,start,edges)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(scalar)) then
      select case (var_shape)
         case(POINT)
            r4 = scalar
            iret = NF90_PUT_VAR(ncid,id,r4)
         case(T_SHAPE)
            start(1) = recnum; edges(1) = 1
            dum(1)=scalar
            iret = NF90_PUT_VAR(ncid,id,dum,start,edges)
         case(XYT_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = recnum; edges(3) = 1
            dum(1)=scalar
            iret = NF90_PUT_VAR(ncid,id,dum,start,edges)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(array) .OR. present(array2d) .OR. present(array3d) ) then
      select case (var_shape)
         case(Z_SHAPE,G_SHAPE)
            start(1) = 1;   edges(1) = nbox
            dum(1:nbox) = array(1:nbox)
            iret = NF90_PUT_VAR(ncid,id,dum(1:nbox),start,edges)
         case(XY_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            iret = NF90_PUT_VAR(ncid,id,real(array2d(:,:),4),start,edges)
         case(XYZ_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = 1;   edges(3) = nbox
            iret = NF90_PUT_VAR(ncid,id,real(array3d(:,:,:),4),start,edges)
         case(XYZT_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = 1;   edges(3) = nbox
            start(4) = recnum; edges(4) = 1
            iret = NF90_PUT_VAR(ncid,id,real(array3d(:,:,:),4),start,edges)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(garray)) then
      select case (var_shape)
         case(OCET_SHAPE)
            start(1) = 1;   edges(1) = ocepoint_len
            start(2) = recnum; edges(2) = 1
         case(SURFT_SHAPE)
            start(1) = 1;   edges(1) = surfpoint_len
            start(2) = recnum; edges(2) = 1
         case(BOTT_SHAPE)
            start(1) = 1;   edges(1) = botpoint_len
            start(2) = recnum; edges(2) = 1
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
      dum(1:nbox)=garray(1:nbox)
      iret = NF90_PUT_VAR(ncid,id,dum,start,edges)
   else
   end if
   call check_err(iret, 'store_data')
   store_data = iret
   return
   end function store_data
!EOC

#if defined INCLUDE_BEN && defined INCLUDE_BENPROFILES
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Definine extra dimension variables
!
! !INTERFACE:
   integer function special_dims(mode,ncid,nlev,name,extname,units, &
                                 time_dim,vars_id)
!
! !DESCRIPTION:
! This is a spcialized routine for the storage of diagnostic variables 
! with  alternative dimensions.
! The typical example are the benthic profiles, which have a sigma
! layer grid with nlev levels.
! 2 additional dimension variables, one with the sigma levels and 
! one for the data points, which is a compressed coordinate
!
! !USES:
   use mem, only: seddepth
   use netcdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: mode
   integer, intent(in)                 :: ncid
   integer, intent(in)                 :: nlev
   character(*), intent(in)            :: name
   character(*), intent(in)            :: extname
   character(*), intent(in)            :: units
   integer, intent(inout)              :: vars_id
   integer, intent(in)                 :: time_dim
!
!  Generic BFM version: Marcello Vichi
!
! !LOCAL VARIABLES:
   real(RLEN)                  :: zz,r,s
   integer                   :: i,j,n,status,altZ_id,dim_altZ
   integer                   :: benprofpoint_dim,benprofpoint_id
   real(RLEN)                  :: arr(0:nlev)
   character(len=30)         :: altZ,altZ_longname
   character(len=6)          :: dum,alt_unit
   integer                   :: dims(2)
!EOP
       if ( index(extname,'__Z' ) ==1 ) then
          j=index(extname,':')-1
          read(extname(1:j),*) dum,altZ, zz,alt_unit, altZ_longname
          ! check is done on the compressed dimension
          status = NF90_INQ_DIMID(ncid, 'benprofpoint', benprofpoint_dim)
          if (status.ne.NF90_NOERR) then
            ! define additional dimensions
            status=NF90_DEF_DIM(ncid,altZ,nlev,dim_altZ)
            if (status.eq.NF90_NOERR) then
               status = NF90_DEF_VAR(ncid,altZ,NF90_REAL,dim_altZ,altZ_id)
               if (status.eq.NF90_NOERR) then
                  i=len_trim(altZ_longname);
                  i=index(extname(1:j),altZ_longname(1:i))
                  status= set_attributes(ncid,altZ_id,long_name=extname(i:j))
                  status= set_attributes(ncid,altZ_id,units=alt_unit)
                  status = NF90_ENDDEF(ncid)
                  status = store_data(ncid,altZ_id,Z_SHAPE,nlev,array=seddepth)
                  status = NF90_REDEF(ncid)
               endif
            endif
            ! define coordinate dimension for benthic profile data
            status=NF90_DEF_DIM(ncid,'benprofpoint',max(NO_BOXES,1),benprofpoint_dim)
            if (status.eq.NF90_NOERR) then
               status = NF90_DEF_VAR(ncid,'benprofpoint',NF90_INT, &
                       benprofpoint_dim, benprofpoint_id)
               if (status.eq.NF90_NOERR) then
                  status = set_attributes(ncid,benprofpoint_id, &
                                          formula_term='alternative sigma point')
                  status = set_attributes(ncid,benprofpoint_id, &
                                          compress='x y ' // altZ)
                  status = NF90_ENDDEF(ncid)
                  !MAV this should be improved
                  status = store_data(ncid,benprofpoint_id,G_SHAPE,NO_BOXES, &
                                      iarray=(/(i,i=1,NO_BOXES)/))
                  status = NF90_REDEF(ncid)
               endif
            endif
          endif
          if ( mode.eq.1) return
          dims(1) = benprofpoint_dim
          dims(2) = time_dim
          status = NF90_DEF_VAR(ncid,name,NF90_REAL,dims,vars_id)
          status= set_attributes(ncid,vars_id,long_name=trim(extname(j+2:)))
          status= set_attributes(ncid,vars_id,units=units)
          special_dims=1
       else
          special_dims=0
       endif
   end function special_dims
#endif

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: check_err() - error reporting on netcdf operations 

   subroutine check_err(iret,filename)
!
! !DESCRIPTION:
!
! !USES
   use netcdf
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer iret
   character(len=*),optional,intent(IN) :: filename
!
!-----------------------------------------------------------------------
!BOC
   if (iret .ne. NF90_NOERR) then

     FATAL "Access to file (or caller): ", trim(filename)
     FATAL "NetCDF Error Message : ",   NF90_STRERROR(iret)
     stop

   endif

   return
   end subroutine check_err
!-----------------------------------------------------------------------

   end module netcdf_bfm


!-----------------------------------------------------------------------
! Copyright 2013 BFM System Team (bfm_st@lists.cmcc.it)
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

