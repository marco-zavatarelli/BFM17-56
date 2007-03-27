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
!EOP
!-----------------------------------------------------------------------
!BOC
            ! compute average values
            ! MAV: to be upgraded soon to standard BFM function calcmean_bfm
            call prepare_bio_output(1,nlev,_ZERO_)

            if(first) then
               first = .false.
               iret = define_mode(ncid,.true.)
               if (bio_setup/=2) then 
               dims(1) = lon_dim
               dims(2) = lat_dim
               dims(3) = z_dim
               dims(4) = time_dim
               do n=stPelStateS,stPelFluxE
                  if ( var_ids(n) /= 0 )  then 
                     iret = new_nc_variable(ncid,var_names(n),NF_REAL, &
                                            4,dims,var_ids(n))
                     iret = set_attributes(ncid,var_ids(n),       &
                                           units=var_units(n),    &
                                           long_name=var_long(n))
                  end if
               end do
            end if !first
            if (bio_setup>1) then ! define benthic variables
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
            end if   
            iret = define_mode(ncid,.false.)
         end if

         do n=stPelStateS,stPelStateE
            if ( (var_ids(n) > 0) .and. (.not.var_ave(n) )) &
              iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=cc(n,:))
         end do
         i=0
         do n=stPelDiagS,stPelDiagE
            i=i+1
            if ( (var_ids(n) > 0).and. (.not.var_ave(n) ) ) & 
               iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=diag(i,:))
         end do

         i=0
         do n=stPelFluxS,stPelFluxE
            i=i+1
            if ( (var_ids(n) > 0)  .and. (.not.var_ave(n))) then
               call make_flux_output(1,i,nlev, h, c1dim)
               iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=c1dim)
            end if
         end do
         j=0
         do n=stPelStateS,stPelFluxE
            if ( (var_ids(n) > 0) .and.var_ave(n) ) then
              j=j+1
              iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev,array=cc_ave(j,:))
            endif
         end do

! storage of benthic variables
! stored as scalar: to be modified if benvar are arrays
         if (bio_setup>1) then
           i=0
           do n=stBenStateS,stBenStateE
              i=i+1
              if ( (var_ids(n) > 0)  .and. (.not.var_ave(n))) &
               iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=ccb(i,1))
           end do
           i=0
           do n=stBenDiagS,stBenDiagE
             i=i+1
             if ( (var_ids(n) > 0)  .and. (.not.var_ave(n))) &
               iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=diagb(i,1))
           end do
           i=0
           do n=stBenFluxS,stBenFluxE
             i=i+1
             if ( (var_ids(n) > 0)  .and. (.not.var_ave(n))) then
               call make_flux_output(2,i,nlev, h, c1dim)
               iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=c1dim(1))
             endif
           end do 
           j=0
           do n=stBenStateS,stBenFluxE
              if ( (var_ids(n) > 0) .and. var_ave(n)) then
                 j=j+1
                 iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,scalar=ccb_ave(j,1))
              endif
           end do
         end if                

   return
   end subroutine bio_save_bfm
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License 
! www.gnu.org
!-----------------------------------------------------------------------
