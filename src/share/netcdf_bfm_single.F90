#include"cppdefs.h"
#define REAL_4B real(4)
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: netcdf_bfm_single --- Save the BFM results in NetCDF
!
! !INTERFACE:
   module netcdf_bfm
!
! !DESCRIPTION:
!  This module provides routines for saving the results using
!  NetCDF format.
!  This is an alternative file to store restarts in single precision
!
! !USES:
   use api_bfm
   use mem,     only: NO_BOXES,NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z,NO_BOXES_XY,Depth
   use global_mem, only: RLEN
   implicit none

   include 'netcdf.inc'
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
   integer                       :: time_len=NF_UNLIMITED
   !---------------------------------------------
   ! Dimension IDs
   !---------------------------------------------
   integer                            :: lon_dim,lat_dim,depth_dim
   integer                            :: x_dim,y_dim
   integer                            :: ocepoint_dim
   integer                            :: surfpoint_dim,botpoint_dim
   integer                            :: time_dim
   integer, parameter                 :: dim1=1,dim4=4
   integer                            :: dims(dim4)
   !---------------------------------------------
   ! Coordinate variables IDs
   !---------------------------------------------
   integer          :: lon_id,lat_id,z_id,z1_id,time_id
   integer          :: zeta_id, mask_id
   integer          :: depth_id,ocepoint_id,surfpoint_id,botpoint_id
   !---------------------------------------------
   ! Restart file ids
   !---------------------------------------------
   integer,public                :: ncid_rst
   integer                       :: ncid_rst_in
   integer                       :: ocepoint_rdim
   integer                       :: surfpoint_rdim,botpoint_rdim
   integer                       :: d3vars_rdim,d3state_rid,d3stateb_rid
   integer                       :: d2vars_rdim,d2state_rid,d2stateb_rid

!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!  Modifications and BFM additions: Marcello Vichi
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
   subroutine init_netcdf_bfm(title,start_time,time_unit,lat,lon,z,dz, &
                              lat2d,lon2d,oceanpoint,surfacepoint,     &
                              bottompoint,mask3d)
!
! !DESCRIPTION:
!  Prepare the netcdf output file which is finalized in init_save_bfm
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*), intent(in)                 :: title,start_time
   integer, intent(in)                          :: time_unit
   real(RLEN), intent(in),optional                :: lat,lon
   real(RLEN), intent(in),dimension(:,:),optional :: lat2d,lon2d
   real(RLEN), intent(in),dimension(:),optional   :: z,dz
   integer, intent(in),dimension(:),optional    :: oceanpoint
   integer, intent(in),dimension(:),optional    :: surfacepoint,bottompoint
   real(RLEN),intent(in),dimension(:,:,:),optional:: mask3d
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
   integer                   :: time_len=NF_UNLIMITED
!!
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'init_netcdf_bfm'

   !---------------------------------------------
   ! Prepare the netcdf file
   !---------------------------------------------
   ext = 'nc'
   fname = TRIM(out_dir) //'/'// TRIM(title) // '.' // ext
   LEVEL2 'Output in NetCDF (time unit is set to seconds):'
   LEVEL2 TRIM(fname)
   iret = nf_create(fname,NF_CLOBBER,ncid_bfm)
   call check_err(iret)

   ncdf_time_unit = time_unit

   !---------------------------------------------
   ! define dimensions
   !---------------------------------------------
   if (present(lon).AND.present(lat)) then
      iret = nf_def_dim(ncid_bfm, 'lon', NO_BOXES_XY, lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid_bfm, 'lat', NO_BOXES_XY, lat_dim)
      call check_err(iret)
   else if (present(lon2d).AND.present(lat2d)) then
      iret = nf_def_dim(ncid_bfm, 'x', NO_BOXES_X, x_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid_bfm, 'y', NO_BOXES_Y, y_dim)
      call check_err(iret)
   else
      stop '### init_netcdf_bfm: lat and lon must be given'
   end if
   iret = nf_def_dim(ncid_bfm, 'z', NO_BOXES_Z, depth_dim)
   call check_err(iret)
   iret = nf_def_dim(ncid_bfm, 'oceanpoint', NO_BOXES, ocepoint_dim)
   call check_err(iret)
   iret = nf_def_dim(ncid_bfm, 'surfacepoint', NO_BOXES_XY, surfpoint_dim)
   call check_err(iret)
   iret = nf_def_dim(ncid_bfm, 'bottompoint', NO_BOXES_XY, botpoint_dim)
   call check_err(iret)
   iret = nf_def_dim(ncid_bfm, 'time', NF_UNLIMITED, time_dim)
   call check_err(iret)

   !---------------------------------------------
   ! define coordinate variables
   !---------------------------------------------
   if (present(lon)) then
      dims(1) = lon_dim
      ndims = 1
   elseif (present(lon2d)) then
      dims(1) = x_dim
      dims(2) = y_dim
      ndims = 2
   end if
   iret = nf_def_var(ncid_bfm,'lon',NF_REAL,ndims,dims,lon_id)
   call check_err(iret)
   if (present(lat)) then
      dims(1) = lat_dim
      ndims = 1
   elseif (present(lat2d)) then
      dims(1) = x_dim
      dims(2) = y_dim
      ndims = 2
   end if
   iret = nf_def_var(ncid_bfm,'lat',NF_REAL,ndims,dims,lat_id)
   call check_err(iret)
   dims(1) = depth_dim
   iret = nf_def_var(ncid_bfm,'z',NF_REAL,1,dims,depth_id)
   call check_err(iret)
   dims(1) = ocepoint_dim
   iret = nf_def_var(ncid_bfm,'oceanpoint',NF_INT,1,dims,ocepoint_id)
   call check_err(iret)
   dims(1) = surfpoint_dim
   iret = nf_def_var(ncid_bfm,'surfacepoint',NF_INT,1,dims,surfpoint_id)
   call check_err(iret)
   dims(1) = botpoint_dim
   iret = nf_def_var(ncid_bfm,'bottompoint',NF_INT,1,dims,botpoint_id)
   call check_err(iret)
   dims(1) = time_dim
   iret = nf_def_var(ncid_bfm,'time',NF_REAL,1,dims,time_id)
   call check_err(iret)
   !---------------------------------------------
   ! define mask variables
   !---------------------------------------------
   if (present(mask3d)) then
      dims(1) = x_dim
      dims(2) = y_dim
      dims(3) = depth_dim
      iret = nf_def_var(ncid_bfm,'mask',NF_REAL,3,dims,mask_id)
      call check_err(iret)
   end if

   !---------------------------------------------
   ! assign attributes
   !---------------------------------------------
   !  coordinates
   iret = set_attributes(ncid_bfm,lon_id,units='degrees_east')
   iret = set_attributes(ncid_bfm,lat_id,units='degrees_north')
   iret = set_attributes(ncid_bfm,depth_id,units='meters')
#ifndef NOT_STANDALONE
   iret = set_attributes(ncid_bfm,ocepoint_id,formula_term='water points')
   iret = set_attributes(ncid_bfm,ocepoint_id,compress='none')
   iret = set_attributes(ncid_bfm,botpoint_id,formula_term='bottom points')
   iret = set_attributes(ncid_bfm,botpoint_id,compress='none')
   iret = set_attributes(ncid_bfm,surfpoint_id,formula_term='surface points')
   iret = set_attributes(ncid_bfm,surfpoint_id,compress='none')
#endif
#ifdef BFM_GOTM
   iret = set_attributes(ncid_bfm,ocepoint_id,formula_term='watercolumn levels')
   iret = set_attributes(ncid_bfm,ocepoint_id,compress='z')
   iret = set_attributes(ncid_bfm,surfpoint_id,formula_term='watercolumn surface')
   iret = set_attributes(ncid_bfm,surfpoint_id,compress='z')
   iret = set_attributes(ncid_bfm,botpoint_id,formula_term='watercolumn bottom')
   iret = set_attributes(ncid_bfm,botpoint_id,compress='z')
#endif
#ifdef BFM_NEMO
   iret = set_attributes(ncid_bfm,ocepoint_id,formula_term='water points')
   iret = set_attributes(ncid_bfm,ocepoint_id,compress='x y z')
   iret = set_attributes(ncid_bfm,botpoint_id,formula_term='bottom points')
   iret = set_attributes(ncid_bfm,botpoint_id,compress='x y z')
   iret = set_attributes(ncid_bfm,surfpoint_id,formula_term='surface points')
   iret = set_attributes(ncid_bfm,surfpoint_id,compress='x y z')
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
   iret = set_attributes(ncid_bfm,time_id,units=trim(ncdf_time_str))

   !---------------------------------------------
   !  global attributes
   !---------------------------------------------
   iret = nf_put_att_text(ncid_bfm,NF_GLOBAL,  &
                          'Title',LEN_TRIM(title),title)
   history = 'Created by BFM v. '//RELEASE
   iret = nf_put_att_text(ncid_bfm,NF_GLOBAL,'history',  &
                          LEN_TRIM(history),history)
   iret = nf_put_att_text(ncid_bfm,NF_GLOBAL,'Conventions',6,'CF-1.0')
   call check_err(iret)

   !---------------------------------------------
   ! leave define mode
   !---------------------------------------------
   iret = nf_enddef(ncid_bfm)
   call check_err(iret)

   !---------------------------------------------
   ! save coordinate variables
   !---------------------------------------------
   if (present(lon)) &
      iret = store_data(ncid_bfm,lon_id,POINT,1,scalar=lon)
   if (present(lat)) &
      iret = store_data(ncid_bfm,lat_id,POINT,1,scalar=lat)
   if (present(lat2d)) &
      iret = store_data(ncid_bfm,lat_id,XY_SHAPE,NO_BOXES_Z, &
                        array2d=lat2d)
   if (present(lon2d)) &
      iret = store_data(ncid_bfm,lon_id,XY_SHAPE,NO_BOXES_Z, &
                        array2d=lon2d)
   if (present(z)) &
      iret = store_data(ncid_bfm,depth_id,Z_SHAPE,NO_BOXES_Z,array=z)
   if (present(dz)) &
      iret = store_data(ncid_bfm,depth_id,Z_SHAPE,NO_BOXES_Z,array=dz)
   if (present(oceanpoint)) &
      iret = store_data(ncid_bfm,ocepoint_id,G_SHAPE,NO_BOXES,iarray=oceanpoint)
   if (present(bottompoint)) &
      iret = store_data(ncid_bfm,botpoint_id,G_SHAPE,NO_BOXES_XY,iarray=bottompoint)
   if (present(surfacepoint)) &
      iret = store_data(ncid_bfm,surfpoint_id,G_SHAPE,NO_BOXES_XY,iarray=surfacepoint)
   if (present(mask3d)) &
      iret = store_data(ncid_bfm,mask_id,XYZ_SHAPE,NO_BOXES_Z, &
                        array3d=mask3d)

   !---------------------------------------------
   ! syncronize
   !---------------------------------------------
   iret = nf_sync(ncid_bfm)
   call check_err(iret)

   end subroutine init_netcdf_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the netcdf restart
!
! !INTERFACE:
   subroutine init_netcdf_rst_bfm(title)
!
! !DESCRIPTION:
!  Prepare the netcdf restart file for the BFM
!
! !USES:
   use mem, only: NO_D3_BOX_STATES, NO_BOXES,    &
                  NO_BOXES_XY, NO_D2_BOX_STATES
   implicit none
!
! !INPUT/OUTPUT PARAMETERS:
   character(len=*), intent(in)                 :: title
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
!!
!-------------------------------------------------------------------------
!BOC
   !---------------------------------------------
   ! Prepare the netcdf file
   !---------------------------------------------
   ext = 'nc'
   fname = './out_'// TRIM(title) // '.' // ext
   LEVEL2 'Restart file in NetCDF'
   LEVEL2 TRIM(fname)
   iret = nf_create(fname,NF_CLOBBER,ncid_rst)
   call check_err(iret)

   !---------------------------------------------
   ! define 3D dimensions and variables
   !---------------------------------------------
   iret = nf_def_dim(ncid_rst, 'd3vars', NO_D3_BOX_STATES, d3vars_rdim)
   call check_err(iret)
   iret = nf_def_dim(ncid_rst, 'oceanpoint', NO_BOXES, ocepoint_rdim)
   call check_err(iret)
   iret = nf_def_dim(ncid_rst, 'surfacepoint', NO_BOXES_XY, surfpoint_rdim)
   call check_err(iret)
   dims(1) = d3vars_rdim
   dims(2) = ocepoint_rdim
   ndims = 2
   iret = nf_def_var(ncid_rst,'D3STATE',NF_FLOAT,ndims,dims,d3state_rid)
   call check_err(iret)
#ifdef BFM_POM
   iret = nf_def_var(ncid_rst,'D3STATEB',NF_FLOAT,ndims,dims,d3stateb_rid)
   call check_err(iret)
#endif

#ifdef INCLUDE_BEN
   !---------------------------------------------
   ! define 2D dimensions and variables
   !---------------------------------------------
   iret = nf_def_dim(ncid_rst, 'd2vars', NO_D2_BOX_STATES, d2vars_rdim)
   call check_err(iret)
   iret = nf_def_dim(ncid_rst, 'bottompoint', NO_BOXES_XY, botpoint_rdim)
   call check_err(iret)
   dims(1) = d2vars_rdim
   dims(2) = botpoint_rdim
   ndims = 2
   iret = nf_def_var(ncid_rst,'D2STATE',NF_FLOAT,ndims,dims,d2state_rid)
   call check_err(iret)
#ifdef BFM_POM
   iret = nf_def_var(ncid_rst,'D2STATEB',NF_FLOAT,ndims,dims,d2stateb_rid)
   call check_err(iret)
#endif
#endif

   !---------------------------------------------
   ! leave define mode
   !---------------------------------------------
   iret = nf_enddef(ncid_rst)
   call check_err(iret)

end subroutine init_netcdf_rst_bfm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Store the restart file
!
! !INTERFACE:
  subroutine save_rst_bfm()
!
! !DESCRIPTION:
! output restart file of BFM variables 
!
! !USES:
   use mem, only: D3STATE, NO_D3_BOX_STATES, NO_BOXES
#ifdef BFM_POM
   use api_bfm, only: D3STATEB
#endif
#ifdef INCLUDE_BEN
   use mem, only: D2STATE, NO_D2_BOX_STATES, NO_BOXES_XY
#ifdef BFM_POM
   use api_bfm, only: D2STATEB
#endif
#endif
   implicit none
!
! !INPUT PARAMETERS:

! !LOCAL VARIABLES:
   integer                   :: iret
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi (INGV) 
!
!EOP
!-----------------------------------------------------------------------
!BOC
     start(1) = 1;   edges(1) = NO_D3_BOX_STATES
     start(2) = 1;   edges(2) = NO_BOXES
     iret = nf_put_vara_real(ncid_rst,d3state_rid,start,edges,D3STATE(:,:))
     call check_err(iret)
#ifdef BFM_POM
     iret = nf_put_vara_real(ncid_rst,d3stateb_rid,start,edges,D3STATEB(:,:))
     call check_err(iret)
#endif
#ifdef INCLUDE_BEN
     start(1) = 1;   edges(1) = NO_D2_BOX_STATES
     start(2) = 1;   edges(2) = NO_BOXES_XY
     iret = nf_put_vara_real(ncid_rst,d2state_rid,start,edges,D2STATE(:,:))
     call check_err(iret)
#ifdef BFM_POM
     iret = nf_put_vara_real(ncid_rst,d2state_rid,start,edges,D2STATEB(:,:))
     call check_err(iret)
#endif
#endif
     LEVEL1 'Restart has been written in NetCDF'
     ! the file is closed in the main (in case of more restart files)
     !iret = nf_close(ncid_rst)
     !call check_err(iret)

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
#ifdef BFM_POM
   use api_bfm, only: D3STATEB
#endif
#ifdef INCLUDE_BEN
   use mem, only: D2STATE, NO_D2_BOX_STATES, NO_BOXES_XY
#ifdef BFM_POM
   use api_bfm, only: D2STATEB
#endif
#endif
   implicit none
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)                 :: title
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
   !---------------------------------------------
   ! open the netcdf restart file
   !---------------------------------------------
   ext = 'nc'
   fname = './in_'// TRIM(title) // '.' // ext
   LEVEL2 'Reading Restart file in NetCDF'
   LEVEL2 TRIM(fname)
   iret = nf_open(fname,NF_NOWRITE,ncid_rst_in)
   call check_err(iret)
   !---------------------------------------------
   ! Check 3D dimensions 
   !---------------------------------------------
   iret = nf_inq_dimid(ncid_rst_in,"d3vars",nstate_id)
   call check_err(iret)
   iret = nf_inq_dimlen(ncid_rst_in,nstate_id,nstate_len)
   call check_err(iret)
   iret = nf_inq_dimid(ncid_rst_in,"oceanpoint",ncomp_id)
   call check_err(iret)
   iret = nf_inq_dimlen(ncid_rst_in,ncomp_id,ncomp_len)
   call check_err(iret)
   if (nstate_len/=NO_D3_BOX_STATES .OR. ncomp_len/=NO_BOXES) then
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
   iret = nf_inq_varid(ncid_rst_in,"D3STATE",nstate_id)
   call check_err(iret)
   iret = nf_get_var_real(ncid_rst_in,nstate_id,D3STATE(:,:))
   call check_err(iret)
#ifdef BFM_POM
   iret = nf_inq_varid(ncid_rst_in,"D3STATEB",nstate_id)
   call check_err(iret)
   iret = nf_get_var_real(ncid_rst_in,nstate_id,D3STATEB(:,:))
   call check_err(iret)
#endif

#ifdef INCLUDE_BEN
   !---------------------------------------------
   ! Check 2D dimensions 
   !---------------------------------------------
   iret = nf_inq_dimid(ncid_rst_in,"d2vars",nstate_id)
   call check_err(iret)
   iret = nf_inq_dimlen(ncid_rst_in,nstate_id,nstate_len)
   call check_err(iret)
   iret = nf_inq_dimid(ncid_rst_in,"bottompoint",ncomp_id)
   call check_err(iret)
   iret = nf_inq_dimlen(ncid_rst_in,ncomp_id,ncomp_len)
   call check_err(iret)
   if (nstate_len/=NO_D2_BOX_STATES .OR. ncomp_len/=NO_BOXES_XY) then
      LEVEL1 '2D Dimension mismatch in restart file:'
      LEVEL2 TRIM(fname)
      LEVEL3 "NO_D2_BOX_STATES in model:", NO_D2_BOX_STATES
      LEVEL3 "NO_D2_BOX_STATES in file:", nstate_len
      LEVEL3 "NO_BOXES_XY in model:", NO_BOXES_XY
      LEVEL3 "NO_BOXES_XY in file:", ncomp_len
      stop 'STOP in read_rst_bfm contained in netcdf_bfm.F90'
   end if
   !---------------------------------------------
   ! Initialize 2D variable
   !---------------------------------------------
   iret = nf_inq_varid(ncid_rst_in,"D2STATE",nstate_id)
   call check_err(iret)
   iret = nf_get_var_real(ncid_rst_in,nstate_id,D2STATE(:,:))
   call check_err(iret)
#ifdef BFM_POM
   iret = nf_inq_varid(ncid_rst_in,"D2STATEB",nstate_id)
   call check_err(iret)
   iret = nf_get_var_real(ncid_rst_in,nstate_id,D2STATEB(:,:))
   call check_err(iret)
#endif
#endif

   LEVEL1 'Finished reading Restart file'
   LEVEL2 TRIM(fname)
   iret = nf_close(ncid_rst_in)
   call check_err(iret)

  end subroutine read_rst_bfm
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

   LEVEL1 'init_save_bfm'
   !---------------------------------------------
   ! Enter define mode
   !---------------------------------------------
   iret = define_mode(ncid_bfm,.true.)

!MAV   if (bio_setup/=2) then
      dims(1) = ocepoint_dim
      dims(2) = time_dim
      do n=stPelStateS,stPelFluxE
         j=0
#ifdef INCLUDE_BENPROFILES
         ! this is a special part for the variable with
         ! alternative dimensions (eg benthic profiles)
         if ( var_ids(n) /= 0 )  &
            j=special_dims(2,ncid_bfm,NO_BOXES_Z,var_names(n),var_long(n),  &
                           var_units(n),time_dim,var_ids(n))
#endif
         if ( j.eq.0 .and. var_ids(n) /= 0 )  then
            iret = new_nc_variable(ncid_bfm,var_names(n),NF_REAL, &
                                 2,dims,var_ids(n))
            iret = set_attributes(ncid_bfm,var_ids(n),            &
                                 units=var_units(n),              &
                                 long_name=var_long(n))
         end if 
      end do
!MAV   endif

!MAV   if (bio_setup>1) then
      dims(1) = botpoint_dim
      dims(2) = time_dim
      do n=stBenStateS,stBenFluxE
         if ( var_ids(n) /= 0 )  then 
            iret = new_nc_variable(ncid_bfm,var_names(n),NF_REAL, &
                                    2,dims,var_ids(n))
            iret = set_attributes(ncid_bfm,var_ids(n),            &
                                 units=var_units(n),         &
                                 long_name=var_long(n)) 
         endif
      end do 
!MAV    end if

   iret = define_mode(ncid_bfm,.false.)
   LEVEL2 'NetCDF definitions completed.'

   !---------------------------------------------
   ! Store the initial conditions
   !---------------------------------------------
   ltime=0.0
   call save_bfm(ltime)

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
   use mem, only: D3STATE,D3DIAGNOS,D2DIAGNOS
#ifdef INCLUDE_BEN
   use mem, only: D2STATE
#endif
   implicit none
!
! !INPUT PARAMETERS:
   real(RLEN),intent(in)     :: time
! !LOCAL VARIABLES:
   integer                   :: iret
   integer                   :: i,j,k,n
   real(RLEN)                  :: temp_time
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  Adapted to BFM: Marcello Vichi (INGV) & Piet Ruardij (NIOZ)
!
!EOP
!-----------------------------------------------------------------------
!BOC

#ifdef DEBUG
   LEVEL1 'save_bfm',time
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
   ! Store snapshot of pelagic variables
   !---------------------------------------------
   do n=stPelStateS,stPelStateE
     if ( var_ids(n) > 0 ) & 
       iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=D3STATE(n,:))
   end do
   !---------------------------------------------
   ! Store snapshot of pelagic diagnostics
   !---------------------------------------------
   i=0
   do n=stPelDiagS,stPelDiagE
      i=i+1
      if ( var_ids(n) > 0 .and. (.not.var_ave(n))) & 
         iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=D3DIAGNOS(i,:))
   end do
#ifndef D1SOURCE
   !---------------------------------------------
   ! Store snapshot of pelagic fluxes
   !---------------------------------------------
   i=0
   do n=stPelFluxS,stPelFluxE
      i=i+1
      if ( var_ids(n) > 0 ) then
         call make_flux_output(1,i,1,NO_BOXES,c1dim)
         iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=c1dim)
      end if
   end do
#endif
   !---------------------------------------------
   ! Store mean values of (any) pelagic entity
   !---------------------------------------------
   k=0
   do n=stPelStateS,stPelFluxE
      if ( var_ids(n) > 0 .and.var_ave(n) ) then
         k=k+1
         iret = store_data(ncid_bfm,var_ids(n),OCET_SHAPE,NO_BOXES,garray=D3ave(k,:))
      endif
   end do

! MAV: we need to solve the storage of 2D diagnostics
!      going through all these loops is probably too expensive
!   if (bio_setup>1) then
#ifdef INCLUDE_BEN
      !---------------------------------------------
      ! Store snapshot of benthic variables
      !---------------------------------------------
      i=0
      do n=stBenStateS,stBenStateE
         i=i+1
         if ( var_ids(n) > 0 ) &
            iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2STATE(i,:))
      end do
#endif
      !---------------------------------------------
      ! Store snapshot of benthic diagnostics
      !---------------------------------------------
      i=0
      do n=stBenDiagS,stBenDiagE
         i=i+1
         if ( var_ids(n) > 0 ) &
            iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2DIAGNOS(i,:))
      end do
#ifdef INCLUDE_BEN
      !---------------------------------------------
      ! Store mean values of (any) benthic entity
      !---------------------------------------------
      k=0
      do n=stBenStateS,stBenFluxE
         if ( var_ids(n) > 0 .and. var_ave(n)) then
            k=k+1
            iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=D2ave(k,:))
         endif
      end do 
#endif
!   end if
#ifdef INCLUDE_BEN
#ifndef D1SOURCE
   !---------------------------------------------
   ! Store snapshot of benthic fluxes and pel. fluxes per square meter!
   !---------------------------------------------
   i=0
   do n=stBenFluxS,stBenFluxE
      i=i+1
      if ( var_ids(n) > 0 ) then
          call make_flux_output(2,i,1,NO_BOXES_XY, c1dim)
          iret = store_data(ncid_bfm,var_ids(n),BOTT_SHAPE,NO_BOXES_XY,garray=c1dim)
      endif
   end do
#endif
#endif


   iret = nf_sync(ncid_bfm)
   call check_err(iret)

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

   iret = nf_close(ncid)
   call check_err(iret)

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
      iret = nf_redef(ncid)
   else
      iret = nf_enddef(ncid)
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
   integer function new_nc_variable(ncid,name,data_type,n,dims,id)
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
   integer, intent(in)                 :: data_type,n
   integer, intent(in)                 :: dims(:)
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: id
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret
!
!-----------------------------------------------------------------------
!BOC
   iret = nf_def_var(ncid,name,data_type,n,dims,id)
   call check_err(iret)
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
      len = len_trim(units)
      iret = nf_put_att_text(ncid,id,'units',len,units)
   end if

   if(present(long_name)) then
      len = len_trim(long_name)
      iret = nf_put_att_text(ncid,id,'long_name',len,long_name)
   end if

   if(present(C_format)) then
      len = len_trim(C_format)
      iret = nf_put_att_text(ncid,id,'C_format',len,C_format)
   end if

   if(present(FORTRAN_format)) then
      len = len_trim(FORTRAN_format)
      iret = nf_put_att_text(ncid,id,'FORTRAN_format',len,FORTRAN_format)
   end if

   if(present(compress)) then
      len = len_trim(compress)
      iret = nf_put_att_text(ncid,id,'compress',len,compress)
   end if

   if(present(formula_term)) then
      len = len_trim(formula_term)
      iret = nf_put_att_text(ncid,id,'formula_term',len,formula_term)
   end if

   if(present(valid_min)) then
      vals(1) = valid_min
      iret = nf_put_att_real(ncid,id,'valid_min',NF_FLOAT,1,vals)
   end if

   if(present(valid_max)) then
      vals(1) = valid_max
      iret = nf_put_att_real(ncid,id,'valid_max',NF_FLOAT,1,vals)
   end if

   if(present(valid_range)) then
      vals(1) = valid_range(1)
      vals(2) = valid_range(2)
      iret = nf_put_att_real(ncid,id,'valid_range',NF_FLOAT,2,vals)
   end if

   if(present(scale_factor)) then
      vals(1) = scale_factor
      iret = nf_put_att_real(ncid,id,'scale_factor',NF_FLOAT,1,vals)
   end if

   if(present(add_offset)) then
      vals(1) = add_offset
      iret = nf_put_att_real(ncid,id,'add_offset',NF_FLOAT,1,vals)
   end if

   if(present(FillValue)) then
      vals(1) = FillValue
      iret = nf_put_att_real(ncid,id,'_FillValue',NF_FLOAT,1,vals)
   end if

   if(present(missing_value)) then
      vals(1) = missing_value
      iret = nf_put_att_real(ncid,id,'missing_value',NF_FLOAT,1,vals)
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
            iret = nf_put_var_int(ncid,id,iscalar)
         case(T_SHAPE)
            start(1) = recnum; edges(1) = 1
            idum(1)=iscalar
            iret = nf_put_vara_int(ncid,id,start,edges,idum)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(iarray)) then
      select case (var_shape)
         case(Z_SHAPE,G_SHAPE)
            start(1) = 1;   edges(1) = nbox
            idum(1:nbox)=iarray(1:nbox)
            iret = nf_put_vara_int(ncid,id,start,edges,idum)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(scalar)) then
      select case (var_shape)
         case(POINT)
            r4 = scalar
            iret = nf_put_var_real(ncid,id,r4)
         case(T_SHAPE)
            start(1) = recnum; edges(1) = 1
            dum(1)=scalar
            iret = nf_put_vara_real(ncid,id,start,edges,dum)
         case(XYT_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = recnum; edges(3) = 1
            dum(1)=scalar
            iret = nf_put_vara_real(ncid,id,start,edges,dum)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(array) .OR. present(array2d) .OR. present(array3d) ) then
      select case (var_shape)
         case(Z_SHAPE,G_SHAPE)
            start(1) = 1;   edges(1) = nbox
            dum(1:nbox) = array(1:nbox)
            iret = nf_put_vara_real(ncid,id,start,edges,dum(1:nbox))
         case(XY_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            iret = nf_put_vara_real(ncid,id,start,edges,real(array2d(:,:),4))
         case(XYZ_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = 1;   edges(3) = nbox
            iret = nf_put_vara_real(ncid,id,start,edges,real(array3d(:,:,:),4))
         case(XYZT_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = 1;   edges(3) = nbox
            start(4) = recnum; edges(4) = 1
            iret = nf_put_vara_real(ncid,id,start,edges,real(array3d(:,:,:),4))
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
      iret = nf_put_vara_real(ncid,id,start,edges,dum)
   else
   end if
   call check_err(iret)
   store_data = iret
   return
   end function store_data
!EOC

#ifdef INCLUDE_BENPROFILES
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
   IMPLICIT NONE
#include "netcdf.inc"
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
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!  Generic BFM version: Marcello Vichi
!
! !LOCAL VARIABLES:
   real(RLEN)                  :: zz,r,s
   integer                   :: i,j,n,status,altZ_id,dim_altZ
   integer                   :: benprofpoint_dim,benprofpoint_id
   real(RLEN)                  :: arr(0:nlev)
   character(len=30)         :: altZ,altZ_longname
   character(len=6)          :: dum,alt_unit
   integer                   :: dims(4)
!EOP
       if ( index(extname,'__Z' ) ==1 ) then
          j=index(extname,':')-1
          read(extname(1:j),*) dum,altZ, zz,alt_unit, altZ_longname
          ! check is done on the compressed dimension
          status = nf_inq_dimid(ncid, 'benprofpoint', benprofpoint_dim)
          if (status.ne.NF_NOERR) then
            ! define additional dimensions
            status=nf_def_dim(ncid,altZ,nlev,dim_altZ)
            if (status.eq.NF_NOERR) then
               dims(1)=dim_altZ
               status = nf_def_var(ncid,altZ,NF_REAL,1,dims,altZ_id)
               if (status.eq.NF_NOERR) then
                  i=len_trim(altZ_longname);
                  i=index(extname(1:j),altZ_longname(1:i))
                  status= set_attributes(ncid,altZ_id,long_name=extname(i:j))
                  status= set_attributes(ncid,altZ_id,units=alt_unit)
                  status = nf_enddef(ncid)
                  status = store_data(ncid,altZ_id,Z_SHAPE,nlev,array=seddepth)
                  status = nf_redef(ncid)
               endif
            endif
            ! define coordinate dimension for benthic profile data
            status=nf_def_dim(ncid,'benprofpoint',NO_BOXES,benprofpoint_dim)
            if (status.eq.NF_NOERR) then
               dims(1) = benprofpoint_dim
               status = nf_def_var(ncid,'benprofpoint',NF_INT,1,dims, &
                                   benprofpoint_id)
               if (status.eq.NF_NOERR) then
                  status = set_attributes(ncid,benprofpoint_id, &
                                          formula_term='alternative sigma point')
                  status = set_attributes(ncid,benprofpoint_id, &
                                          compress='x y ' // altZ)
                  status = nf_enddef(ncid)
                  !MAV this should be improved
                  status = store_data(ncid,benprofpoint_id,G_SHAPE,NO_BOXES, &
                                      iarray=(/(i,i=1,NO_BOXES)/))
                  status = nf_redef(ncid)
               endif
            endif
          endif
          if ( mode.eq.1) return
          dims(1) = benprofpoint_dim
          dims(2) = time_dim
          status = nf_def_var(ncid,name,NF_REAL,2,dims,vars_id)
          status= set_attributes(ncid,vars_id,long_name=trim(extname(j+2:)))
          status= set_attributes(ncid,vars_id,units=units)
          special_dims=1
       else
          special_dims=0
       endif
   end function special_dims
#endif

!-----------------------------------------------------------------------

   subroutine check_err(iret)
     integer iret
     include 'netcdf.inc'
     if (iret .ne. NF_NOERR) then
       print *, 'NetCDF Error: ',nf_strerror(iret)
     stop
     endif
   end subroutine check_err
!-----------------------------------------------------------------------

   end module netcdf_bfm


!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
! Additional copyright by M. Vichi (vichi@bo.ingv.it)
!-----------------------------------------------------------------------

