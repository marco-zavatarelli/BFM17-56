#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Store the BFM results
!
! !INTERFACE:
   subroutine bio_save_bfm(nlev)
!
! !DESCRIPTION:
! Write NetCDF output files for BFM
!
! !USES:
   use bio_var
   use output, only: out_fmt,ts
   use ncdfout, only: ncid
   use ncdfout, only: lon_dim,lat_dim,z_dim,time_dim,dims
   use ncdfout, only: define_mode,new_nc_variable,set_attributes,store_data
   IMPLICIT NONE
#include "netcdf.inc"
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!  Modified for BFM by: Piet Ruardj and Marcello Vichi
!
! !LOCAL VARIABLES:
   logical, save             :: first=.true.
   integer                   :: iret
   integer                   :: i,j,n
   REALTYPE                  :: c1dim(0:nlev)

   integer,external          :: special_dims
!EOP
!-----------------------------------------------------------------------
!BOC

   ! compute average values
   ! MAV: to be upgraded soon to standard BFM function calcmean_bfm
   call prepare_bio_output(1,nlev,_ZERO_)
   if (first) then
      first = .false.
      iret = define_mode(ncid,.true.)
      dims(1) = lon_dim
      dims(2) = lat_dim
      dims(3) = z_dim
      dims(4) = time_dim
      do n=stPelStateS,stPelFluxE
         j=0
#ifdef INCLUDE_BENPROFILES
         ! this is a special part for the variables with
         ! alternative dimensions (eg benthic profiles)
         if ( var_ids(n) /= 0 )  &
            j=special_dims(2,ncid,nlev,var_names(n),var_long(n),var_units(n), &
                           lon_dim,lat_dim,time_dim,var_ids(n))
#endif
         if ( j.eq.0 .and. var_ids(n) /= 0 )  then
          iret = new_nc_variable(ncid,var_names(n),NF_REAL, &
                            4,dims,var_ids(n))
          iret = set_attributes(ncid,var_ids(n),       &
                           units=var_units(n),    &
                           long_name=var_long(n))
         end if
      end do
      dims(1) = lon_dim
      dims(2) = lat_dim
      dims(3) = time_dim
      do n=stBenStateS,stBenFluxE
         if ( var_ids(n) /= 0 )  then 
            iret = new_nc_variable(ncid,var_names(n),NF_REAL, &
                                   3,dims,var_ids(n))
            iret = set_attributes(ncid,var_ids(n),       &
                                  units=var_units(n),    &
                                  long_name=var_long(n))
         endif
      end do
      iret = define_mode(ncid,.false.)
   end if ! first
   !---------------------------------------------
   ! Store snapshot of pelagic variables
   !---------------------------------------------
   do n=stPelStateS,stPelStateE
      if (var_ids(n) > 0)  &
         iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=cc(n,:))
   end do
   !---------------------------------------------
   ! Store snapshot of pelagic diagnostics
   !---------------------------------------------
   i=0
   do n=stPelDiagS,stPelDiagE
      i=i+1
      if (var_ids(n) > 0) & 
         iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=diag(i,:))
   end do
   !---------------------------------------------
   ! Store snapshot of pelagic fluxes
   !---------------------------------------------
   i=0
   do n=stPelFluxS,stPelFluxE
      i=i+1
      if (var_ids(n) > 0) then
         call make_flux_output(1,i,0,nlev, c1dim)
         iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=c1dim)
      end if
   end do
   !---------------------------------------------
   ! Store mean values of (any) pelagic entity
   !---------------------------------------------
   j=0
   do n=stPelStateS,stPelFluxE
      if ( (var_ids(n) > 0) .and.var_ave(n) ) then
         j=j+1
         iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=cc_ave(j,:))
      endif
   end do

! MAV: we need to solve the storage of 2D diagnostics
!      going through all these loops is probably too expensive
!   if (bio_setup>1) then
      !---------------------------------------------
      ! Store snapshot of benthic variables
      !---------------------------------------------
      i=0
      do n=stBenStateS,stBenStateE
         i=i+1
         if (var_ids(n) > 0) &
            iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=ccb(i,1))
      end do
      !---------------------------------------------
      ! Store snapshot of 2D diagnostics
      !---------------------------------------------
      i=0
      do n=stBenDiagS,stBenDiagE
         i=i+1
         if (var_ids(n) > 0)   &
            iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=diagb(i,1))
      end do
      !---------------------------------------------
      ! Store mean values of (any) benthic entity
      !---------------------------------------------
      j=0
      do n=stBenStateS,stBenFluxE
         if ( (var_ids(n) > 0) .and. var_ave(n)) then
            j=j+1
            iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=ccb_ave(j,1))
         endif
      end do
      !---------------------------------------------
      ! Store snapshot of benthic fluxes and pel. fluxes per square meter!
      !---------------------------------------------
      i=0
      do n=stBenFluxS,stBenFluxE
         i=i+1
         if (var_ids(n) > 0) then
            call make_flux_output(2,i,nlev, h, c1dim)
            iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=c1dim(1))
         endif
      end do 

   return
   end subroutine bio_save_bfm
!EOC

#ifdef INCLUDE_BENPROFILES
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Definine extra dimension variables
!
! !INTERFACE:
   integer function special_dims(mode,ncid,nlev,name,extname,units, &
                    lon_dim,lat_dim,time_dim,vars_id)
!
! !DESCRIPTION:
! This is a spcialized routine for the storage of diagnostic variables 
! with  alternative dimensions.
! The typical example are the benthic profiles, which have a sigma
! layer grid with the same number of levels as NO_BOXES_Z
!
! !USES:
   use mem, only: seddepth
   use ncdfout, only: set_attributes,store_data
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
   integer, intent(in)                 :: lon_dim
   integer, intent(in)                 :: lat_dim
   integer, intent(in)                 :: time_dim
   integer, intent(inout)              :: vars_id
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!
! !LOCAL VARIABLES:
   REALTYPE,parameter        :: ddu=2.0
   REALTYPE                  :: zz,r,s
   integer                   :: dims(4)
   integer                   :: i,j,n,status,altZ_id,dim_altZ
   REALTYPE                   :: arr(0:nlev)
   character(len=30)         :: altZ,altZ_longname
   character(len=6)          :: dum,alt_unit
!EOP
!-----------------------------------------------------------------------
!BOC
       if ( index(extname,'__Z' ) ==1 ) then
          j=index(extname,':')-1
          read(extname(1:j),*) dum,altZ, zz,alt_unit, altZ_longname
          status = nf_inq_dimid(ncid, altZ, dim_altZ)
          if (status.ne.NF_NOERR) then
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
                  ! seddepth is computed in init_benthic_bfm
                  status = store_data(ncid,altZ_id,Z_SHAPE,nlev,array=seddepth)
                  status = nf_redef(ncid)
               endif
            endif
          endif
          if ( mode.eq.1) return
          dims(1)=lon_dim;dims(2)=lat_dim
          dims(3)=dim_altZ;dims(4)=time_dim
          status = nf_def_var(ncid,name,NF_REAL,4,dims,vars_id)
          status= set_attributes(ncid,vars_id,long_name=trim(extname(j+2:)))
          status= set_attributes(ncid,vars_id,units=units)
          special_dims=1
       else
          special_dims=0
       endif
   end function special_dims
!EOC
#endif
