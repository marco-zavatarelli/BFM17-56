!$Id: bio.F90,v 1.43 2008-03-26 08:56:53 kb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio --- biological model \label{sec:bio}
!
! !INTERFACE:
   module bio
!
! !DESCRIPTION:
! This is the central module for coupling the BFM model into GOTM.
! This file is directly derived from the original bio.F90 in the GOTM
! source directory and all the other biomodels have been removed.
! From here, after reading the namelist file {\tt bio.nml},
! the biogeochemical model is initialised, the memory
! is allocated, the advection and diffusion is called, the ODE solvers
! for the right hand sides are called
! 
! !USES:
   use bio_var

   use bfm_solver
   use bio_bfm, only : init_bio_bfm,var_info_bfm,settling_vel_bfm
   use bio_bfm, only : envforcing_bfm
   use bio_bfm, only : reset_diagonal, allocate_memory_bfm, &
                       test_on_negative_states
   use trace_bdy, only:init_trace_bdy,init_var_trace

   use output, only: out_fmt,write_results,ts,close_output

   use util
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio, set_env_bio, do_bio, get_bio_updates, clean_bio
   logical, public                     :: bio_calc=.false.
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !PRIVATE DATA MEMBERS:
!  from a namelist
   logical                   :: bio_eulerian=.true.
   REALTYPE                  :: cnpar=0.5
   integer                   :: w_adv_discr=6
   integer                   :: ode_method=1
   integer                   :: split_factor=1
   logical                   :: bioshade_feedback=.true.
   logical                   :: bio_lagrange_mean=.true.
   integer                   :: bio_npar=10000
   REALTYPE                  :: depth
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio(namlst,fname,unit,nlev,h)
!
! !DESCRIPTION:
! Here, the bio namelist {\tt bio.nml} is read and memory for the
! Lagrangian part of the model is allocated (note that the
! Lagrangian model up to now only works for the simple suspended matter model).
! If a Lagrangian particle method is chosen, particles are 
! equidistantly distributed. 
! The initial  Furthermore, information on the specific settings are
! written to standard output.
! Finally, the mussel module is called for initialisation.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: h(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: rc,i,j,n
   namelist /bio_nml/ bio_calc,bio_model,bio_eulerian, &
                      bio_setup, &
                      cnpar,w_adv_discr,ode_method,split_factor, &
                      bioshade_feedback,bio_lagrange_mean,bio_npar
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_bio'

   bio_setup=1
   depth=sum(h)

!  Open and read the namelist
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_nml,err=99)
   close(namlst)

   if (bio_calc) then

!     a sanity check (only temporarely)
      if (.not. bio_eulerian) then
         if (bio_model .ne. 3) then
            FATAL "Lagrangian simulations only tested/works with bio_model=3"
         end if
      end if

      select case (bio_model)

      case (101)  ! The BFM model

         call init_bio_bfm(nlev,unit)

         call allocate_memory(nlev)
         call allocate_memory_bfm(nlev)

         call var_info_bfm()

         call init_var_bfm(namlst,'bio_bfm.nml',unit,bio_setup)
         ! this call is needed because benthic initialisation requires 
         ! water-column physical conditions
         call envforcing_bfm(nlev,bioshade_feedback,h)
         call init_benthic_bfm(namlst,'bio_bfm.nml',unit,bio_setup)
         ! initialize averaging of BFM output variables
         ! MAV: to be upgraded soon to standard BFM function calcmean_bfm
         call prepare_bio_output(0,nlev,_ZERO_)

       case (102)  ! The trace model

         call init_trace_bdy('bio_trace.nml',nlev)
         call allocate_memory(nlev)
         call init_var_trace(bio_model)
      case default
         stop "bio: Using the BFM model without a valid biomodel type in bio.nml (101,102)!"
      end select


         LEVEL3 "Using Eulerian solver"
         select case (ode_method)
            case (1)
               LEVEL2 'Using euler_forward()'
            case (2)
               LEVEL2 'Using runge_kutta_2()'
            case (3)
               LEVEL2 'Using runge_kutta_4()'
            case (4)
               LEVEL2 'Using patankar()'
            case (5)
               LEVEL2 'Using patankar_runge_kutta_2()'
            case (6)
               LEVEL2 'Using patankar_runge_kutta_4()'
            case (7)
               LEVEL2 'Using modified_patankar()'
            case (8)
               LEVEL2 'Using modified_patankar_2()'
            case (9)
               LEVEL2 'Using modified_patankar_4()'
            case (10)
               LEVEL2 'Using emp_1()'
            case (11)
               LEVEL2 'Using emp_2()'
            case default
               stop "bio: no valid ode_method specified in bio.nml!"
         end select

   end if


   return

98 LEVEL2 'I could not open bio.nml'
   LEVEL2 'If thats not what you want you have to supply bio.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio.nml'
   bio_calc = .false.
   return
99 FATAL 'I could not read bio.nml'
   stop 'init_bio'
   end subroutine init_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set external variables used by the BIO
! modules
!
! !INTERFACE: 
   subroutine set_env_bio(nlev,h_,t_,s_,rho_,nuh_,rad_,wind_,I_0_, &
                          w_,w_adv_ctr_,u_taub_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
   REALTYPE, intent(in)                :: h_(0:nlev)
   REALTYPE, intent(in)                :: nuh_(0:nlev)
   REALTYPE, intent(in)                :: t_(0:nlev)
   REALTYPE, intent(in)                :: s_(0:nlev)
   REALTYPE, intent(in)                :: rho_(0:nlev)
   REALTYPE, intent(in)                :: rad_(0:nlev)
   REALTYPE, intent(in)                :: wind_
   REALTYPE, intent(in)                :: I_0_
   REALTYPE, optional, intent(in)      :: w_(0:nlev)
   integer, optional, intent(in)       :: w_adv_ctr_
   REALTYPE, optional, intent(in)      :: u_taub_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES
!EOP
!-----------------------------------------------------------------------
!BOC

   h         = h_
   t         = t_
   s         = s_
   rho       = rho_
   nuh       = nuh_
   rad       = rad_
   wind      = wind_
   I_0       = I_0_
   if (present(w_)) w = w_
   if (present(w_adv_ctr_)) w_adv_ctr = w_adv_ctr_
   if (present(u_taub_)) u_taub = u_taub

   return
   end subroutine set_env_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the bio model \label{sec:do-bio}
!
! !INTERFACE:
   subroutine do_bio(nlev,dt)
!
! !DESCRIPTION:
! This is the main loop for the biogeochemical model. Basically 
! an operational split method is used, with first calculating the
! transport part, and than the reaction part.
! During the transport part, all sinks and sources are set to zero,
! and the surface fluxes are computed by calling the
! model specific surface flux subroutine. Then the mussel module
! is called.  For the Eulerian calculation, vertical advection
! (due to settling or rising or vertical migration), vertical advection due
! to physical velocity and vertical
! diffusion (due to mixing) and afterwards the light 
! calculation (for the PAR) and the ODE solver for the right
! hand sides are called. The vertical advection due to settling and
! rising must be conservative,
! which is ensured by setting the local variable {\tt adv\_mode\_1=1},
! see section \ref{sec:advectionMean} on page \pageref{sec:advectionMean}.
! In contrast to this, the vertical advection due to physical velocities must be
! non-conservative, such that for that the local variable {\tt adv\_mode\_0}
! is set to 0, see  see section \ref{sec:advectionMean} on page
! \pageref{sec:advectionMean}.
! It should be noted here that the PAR and the selfshading effect
! is calculated in a similar way for all biogeochemical models
! implemented in GOTM so far. In the temperature equation the
! absorption of solar radiation, $I(z)$, is the only source term,
! see equation (\ref{Iz}) section \ref{sec:temperature}.
! In (\ref{Iz}), a term $B(z)$ due to bioturbidity is used, which 
! is calculated as a function of the biogeochemical particulate
! matter in the water column:
! \begin{equation}\label{B}
! B(z)=\exp\left(-k_c\int_z^0\left(\sum C_{turb}(\xi)\right)\,d\xi\right),
! \end{equation}
! where $k_c$ is the attenuation constant for self shading and 
! $\sum C_{turb}$ is the sum of the biogeochemical particulate 
! matter concentrations.
! The photosynthetically
! available radiation, $I_{PAR}$, follows from
! \begin{equation}
!   \label{light}
!   I_{PAR}(z)=I_0
! (1-a)\exp\left(\frac{z}{\tilde\eta_2}\right)
!   B(z).
! \end{equation}
! 
! For Lagrangian particle calculations, 
! the Lagrangian advection-diffusion routine {\tt lagrange} is called,
! and afterwards, if chosen, the removal of particles due to benthic
! filter feeders (mussels) is done.
! Finally, the calculation of Eulerian concentrations are calculated
! from Lagrangian counts per grid cell for output.
! 
!
! !USES:
   use bio_var, only: I_0_local => I_0
   use gotm_error_msg, only:gotm_error,set_warning_for_getm , &
                            get_parallel_flag_from_getm
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
   REALTYPE, intent(in)                :: dt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer, parameter        :: adv_mode_0=0
   integer, parameter        :: adv_mode_1=1
   REALTYPE                  :: Qsour(0:nlev),Lsour(0:nlev)
   REALTYPE                  :: RelaxTau(0:nlev)
   REALTYPE                  :: dt_eff
   integer                   :: j,n
   integer                   :: split
   integer                   :: i,np
   REALTYPE                  :: filter_depth
   integer, save             :: count=0
   logical, save             :: set_C_zero=.true.
   integer                   :: k
   integer                   :: kt=0
   logical                   :: parallel
   logical                   :: llsumh=.TRUE.
   REALTYPE                  :: c1dim(0:nlev)
!EOP
!-----------------------------------------------------------------------
!BOC
   if (bio_calc) then

      I_0_local = I_0

      Qsour    = _ZERO_
      Lsour    = _ZERO_
      RelaxTau = 1.e15

      select case (bio_model)
         case (101)
            !MAV create a surface_fluxes_bfm 
      end select
 
      ! Store the value of the state variables before transport 
      cc_before_transport=cc

#ifdef BFM_GETM
      call get_parallel_flag_from_getm(parallel)
      if (parallel .and.bio_eulerian.and.bio_model==6) then
        ! Sometimes it happens that a concentration becomes negative
        ! after calculations of the 3d-transport
        ! this check is only done when bottom depth is greater than a certain depth
!MAV: move this constant out of the code asap!
        llsumh=(sum(h(1:nlev)) .gt. 5.0D+00) 
          do j=1,numcc
            if (bio_setup /= 2 ) then
              if (pelvar_type(j)>=ALLTRANSPORT) then
                c1dim = cc(j,:)
                call test_on_negative_states( j,llsumh,c1dim,h,nlev, "GETM advection", kt )
                if ( kt.gt.0) cc(j,:) = c1dim
              endif
            endif
          enddo
      endif
#endif

      ! transfer the particle settling velocities to the BFM
      call settling_vel_bfm

      do j=1,numcc
         if (bio_setup /= 2 ) then
           if (pelvar_type(j)>=ALLTRANSPORT) then
            call test_on_negative_states(j,.TRUE.,cc(j,:),h,nlev,"biology",kt)
            ! MAV: this part could be put behind a flag: stop_on_negative
            if ( kt.lt.0) then
               call bio_save(nlev,_ZERO_)
               call close_output()
               call gotm_error('do_bio', 'negative state value');
               return
            endif

!           do advection step due to settling or rising
!           NOTE: in BFM_GETM this is done only when depth is greater than a
!           certain value and if there is settling
!           llsumh is always true in GOTM
            if ( llsumh .and. llws(j)) then
               call adv_center(nlev,dt,h,h,ws(j,:),flux,           &
                    flux,_ZERO_,_ZERO_,w_adv_discr,adv_mode_1,cc(j,:))
            end if

!           do advection step due to vertical velocity
            if(w_adv_ctr .ne. 0) then
               call adv_center(nlev,dt,h,h,w,flux,                   &
                    flux,_ZERO_,_ZERO_,w_adv_ctr,adv_mode_0,cc(j,:))
            end if
            
!           do diffusion step
            call diff_center(nlev,dt,cnpar,posconc(j),h,Neumann,Neumann,&
                sfl(j),bfl(j),nuh,Lsour,Qsour,RelaxTau,cc(j,:),cc(j,:))
            ! MAV: this part could be put behind a flag: stop_on_negative
            call test_on_negative_states(j,llsumh,cc(j,:),h,nlev,"GOTM physics",kt)
            if ( kt.lt.0) then
              call bio_save(nlev,_ZERO_)
              call close_output()
              call gotm_error('do_bio', 'negative state value');
              return
            endif
          end if  ! BFM ALLTRANSPORT
         end if ! BFM bio_setup
        end do

      do split=1,split_factor
         dt_eff=dt/float(split_factor)

!        Very important for 3D models to save extra 3D field:
         bioshade_=_ONE_

         select case (bio_model)
            case (101)
               ! this flag is only used with GETM
               ! it is always true for GOTM
               if ( llsumh ) then
                  call envforcing_bfm(nlev,bioshade_feedback,h)
                  call ode_solver_bfm(ode_method,nlev,dt_eff)
                  ! accumulate BFM output variables
                  ! MAV: to be upgraded soon to standard BFM function calcmean_bfm
                  call prepare_bio_output(10,nlev,_ZERO_)
                  call prepare_bio_output(11,nlev,_ZERO_)
                  call prepare_bio_output(12,nlev,_ZERO_)
               end if
         end select

      end do

   end if
   return
   end subroutine do_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: return updated variables in the bio modules
! modules
!
! !INTERFACE: 
   subroutine get_bio_updates(nlev,bioshade)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: nlev
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: bioshade(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES
!EOP
!-----------------------------------------------------------------------
!BOC

   if (bioshade_feedback) then
      bioshade = bioshade_
   end if

   return
   end subroutine get_bio_updates
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine clean_bio
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'clean_bio'

   if (allocated(par))            deallocate(par)
   if (allocated(cc))             deallocate(cc)
   if (allocated(ws))             deallocate(ws)
   if (allocated(sfl))            deallocate(sfl)
   if (allocated(bfl))            deallocate(bfl)
   if (allocated(posconc))        deallocate(posconc)
   if (allocated(var_ids))        deallocate(var_ids)
   if (allocated(var_names))      deallocate(var_names)
   if (allocated(var_units))      deallocate(var_units)
   if (allocated(var_long))       deallocate(var_long)

!  The external provide arrays
   if (allocated(h))              deallocate(h)
   if (allocated(nuh))            deallocate(nuh)
   if (allocated(t))              deallocate(t)
   if (allocated(s))              deallocate(s)
   if (allocated(rho))            deallocate(rho)
   if (allocated(rad))            deallocate(rad)
   if (allocated(w))              deallocate(w)
   if (allocated(bioshade_))      deallocate(bioshade_)
   if (allocated(abioshade_))     deallocate(abioshade_)

   if (allocated(llws))           deallocate(llws)
   if (allocated(pp))             deallocate(pp)
   if (allocated(dd))             deallocate(dd)
   if (allocated(ccb))            deallocate(ccb)
   if (allocated(ppb))            deallocate(ppb)
   if (allocated(ddb))            deallocate(ddb)
   if (allocated(pelvar_type))    deallocate(pelvar_type)

   init_saved_vars=.true.

   LEVEL1 'done.'

   return
   end subroutine clean_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory for biological variables
!
! !INTERFACE:
   subroutine allocate_memory(nlev)
!
! !DESCRIPTION:
! Here, the memory for the global biogeochemical parameters
! such as concentrations, settling velocities, surface and bottom
! boundary fluxes, and various other parameters is allocated.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: rc
   integer                   :: numsave
!EOP
!-----------------------------------------------------------------------
!BOC

   allocate(par(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (par)'

   allocate(cc(1:numc,0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (cc)'
   cc=_ZERO_

   allocate(ws(1:numc,0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (ws)'
   ws=_ZERO_

   allocate(sfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (sfl)'
   sfl=_ZERO_

   allocate(bfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (bfl)'
   bfl=_ZERO_

   allocate(posconc(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (posconc)'
   posconc=1

   allocate(llws(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (llws)'
   llws=.false.

   numsave = numc
   allocate(pp(1:numc,1:numc,0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (pp)'
   pp=_ZERO_

   allocate(dd(1:numc,1:numc,0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (dd)'
   dd=_ZERO_

   ! allocate variable holding type and save attributes
   allocate(pelvar_type(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (pelvar_type)'
   pelvar_type = 0

   numsave=numc+numc_diag+numc_flux+numbc+numbc_diag+numbc_flux

   allocate(var_ave(1:numsave),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_ave)'
   var_ave=.false.

   allocate(var_ids(1:numsave),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_ids)'
   var_ids=0;

   allocate(var_names(1:numsave),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_names)'

   allocate(var_units(1:numsave),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_units)'

   allocate(var_long(1:numsave),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating var_long)'

!  The external provide arrays
   allocate(h(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (h)'

   allocate(nuh(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (nuh)'

   allocate(t(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (t)'

   allocate(s(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (s)'

   allocate(rho(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (rho)'

   allocate(rad(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (rad)'

   allocate(w(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (w)'

   allocate(bioshade_(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (bioshade)'

   allocate(abioshade_(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_bio(): Error allocating (abioshade)'

   return
   end subroutine allocate_memory
!EOC
!-----------------------------------------------------------------------

   end module bio

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------
