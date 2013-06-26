#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:routine to calculate sum of fluxes of D*SINK's and D*SOURCE's
!    into diagnostic variables defined as flux-variables in GlobalDefsBFM.model 
!
! !INTERFACE:
      subroutine make_flux_output(mode, nr0,zlev,nlev,out)
!
! !DESCRIPTION:
!     The information provided in the flux definition as given in GlobalDefsBFM.model
!     is transferred into a number of array and scalars (flx_*). All
!     elements of these arrays are defined in AllocateMem.F90. This
!     information is used by this routines to calculates sums of fluxes
!     from D3SINK, D2SINK, D3SOURCE, D2SOURCE, D3STATE and D2STATE ( See
!     ModuleMem.F90 and AllocateMem.F90)
!
!     At input the value of zlev is 0 or 1. The value depend to which
!      model BFM is coupled. The arrays in GOTM start at index 0. In most
!     other coupled models at 1.
!
! !USES:
      use constants, only: RLEN, ZERO, SEC_PER_DAY
      USE global_mem, ONLY: LOGUNIT
#ifdef NOPOINTERS
      use mem
#else
      use mem, only: NO_BOXES_XY,NO_BOXES,NO_BOXES_X,NO_BOXES_Y, &
            NO_BOXES_Z,BoxNumberX,BoxNumberY,BoxNumberZ,BoxNumberXY , &
            BoxNumber
      use mem, only: iiBen
#ifndef D1SOURCE
      use mem, only: flx_calc_nr,flx_CalcIn,flx_t,flx_states, &
              flx_ostates,flx_SS,flx_cal_ben_start,flx_option
#else
      use mem, only: D3FLUX_FUNC, D3FLUX_MATRIX
#ifdef BFM_GOTM
      use bio_var, only: stPelStateS, stPelStateE
#else
       use api_bfm, only: stPelStateS, stPelStateE
#endif
#endif
      use mem, only: D3SOURCE,D3STATE
#ifndef ONESOURCE
      use mem, only: D3SINK
#endif
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
      use mem, only: D2SOURCE,D2STATE
#ifndef ONESOURCE
      use mem, only: D2SINK
#endif
#endif
      use mem, only: PELBOTTOM,PELSURFACE,Depth
#endif

#ifdef BFM_GOTM
      use bio_var, ONLY: SRFindices,BOTindices
#else
      use api_bfm, ONLY: SRFindices,BOTindices
#endif


!
! !INPUT PARAMETERS:
      implicit none
      integer,intent(IN)                  ::mode
      integer,intent(IN)                  ::nr0
      integer,intent(IN)                  ::zlev
      integer,intent(IN)                  ::nlev
#ifdef BFM_GOTM
      real(RLEN),intent(OUT),dimension(0:nlev)  :: out
#else
      real(RLEN),intent(OUT),dimension(nlev)    :: out
#endif

!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij (NIOZ), Marcello Vichi (INGV)
!

!
! !LOCAL VARIABLES:
      integer      ::nr
      integer      ::i
      integer      ::k
#ifdef D1SOURCE
      integer      ::idx_i, idx_j, origin, destination
#endif
#ifndef NOPOINTERS
      integer      ::idummy
#endif
      integer      ::klev
#ifdef BFM_GOTM
      real(RLEN),dimension(0:NO_BOXES) :: hulp
#else
      real(RLEN),dimension(NO_BOXES)   :: hulp
#endif
      real(RLEN)   :: r

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! user defined external functions
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      integer, external  :: D3toD1
      integer, external  :: D2toD1


!EOP
!-----------------------------------------------------------------------
!BOC
! This subroutine is compiled only if sources are stored
      CHARACTER(LEN=80) :: format1, format2
      format1 = '(A,I2,A,1E20.10,A,I2,A,1E20.10)'
      format2 = '(A,I2,A,A,1E20.10,A,1E20.10)'
      write(LOGUNIT,*) "-------MAKE FLUX BEGIN--------------"

#ifndef D1SOURCE

      nr=nr0;if ( mode == 2 ) nr=nr0+flx_cal_ben_start
      klev=NO_BOXES ; if ( flx_CalcIn(nr) == iiBen)  klev=NO_BOXES_XY 
      hulp(:)=ZERO
      if ( flx_CalcIn(nr) == iiBen) then
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
        do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
          ! check if the variable is defined (e.g. in case of O3c)
          ! if not, the flux becomes diagonal
          idummy = flx_ostates(i)
          if (idummy == 0) idummy=flx_states(i)
          if (flx_SS(i) ==1 ) then
#ifdef ONESOURCE
             ! notice that the negative sign is already included in FluxFunctions.F90
             ! (l. 80) thus there is a further change of sign here
             hulp(1:klev)= hulp(1:klev) &
                    - flx_t(i) * D2SOURCE(idummy,flx_states(i),:)
#else
             hulp(1:klev)= hulp(1:klev) &
                    + flx_t(i) * D2SINK(flx_states(i),idummy,:)
#endif
          else
#ifdef ONESOURCE
             ! notice that the negative sign is already included in FluxFunctions.F90
             ! (l. 80) thus there is a further change of sign here
             hulp(1:klev)= hulp(1:klev) &
                    - flx_t(i) * D2SOURCE(idummy,flx_states(i),:)
#else
             hulp(1:klev)= hulp(1:klev) &
                   + flx_t(i) * D2SOURCE(flx_states(i),idummy,:)
#endif
          endif
        enddo
#endif
      else
        do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
          ! check if the variable is defined (e.g. in case of O3c)
          ! if not, the flux becomes diagonal
          idummy = flx_ostates(i)
          if (idummy == 0) idummy=flx_states(i)
          !write(LOGUNIT,*) "FUNC_INI     : ",flx_states(i), flx_ostates(i), "HULP", hulp(:)
          if (flx_SS(i) ==1 ) then
#ifdef ONESOURCE
             ! notice that the negative sign is already included in FluxFunctions.F90
             ! (l. 80) thus there is a further change of sign here
             hulp(1:klev)= hulp(1:klev) &
                    - flx_t(i) * D3SOURCE(idummy,flx_states(i),:)
#else
             hulp(1:klev)= hulp(1:klev) &
                    + flx_t(i) * D3SINK(flx_states(i),idummy,:)
#endif
             !write(LOGUNIT,*) "FUNC_CORR_SS1: ",nr0, flx_states(i), idummy, "HULP", hulp(:)
             ! correcting for fluxes  to other systems
             if ( flx_states(i) == idummy) then
                !write(*,*) "Func=", nr0, " Dir=", flx_SS(i), " Compo1=", flx_states(i), " Compo2=", idummy 
                write(LOGUNIT,format1) "FUNC_SS1(", nr0 , ")->" , hulp(BOTindices), &
                     " PELBOTTOM (", flx_states(i), "): ",  PELBOTTOM(flx_states(i),:)
                write(LOGUNIT,format1) "FUNC_SS1(", nr0 ,  ")->" ,hulp(BOTindices), &
                     " PELSURFACE(", flx_states(i), "): ", PELSURFACE(flx_states(i),:)
                hulp(BOTindices)=hulp(BOTindices)+flx_t(i) *min(ZERO,&
                       PELBOTTOM(flx_states(i),:))/Depth(BOTindices)/SEC_PER_DAY
                hulp(SRFindices)=hulp(SRFindices)+flx_t(i) *min(ZERO,&
                    PELSURFACE(flx_states(i),:))/Depth(SRFindices)/SEC_PER_DAY
                write(LOGUNIT,format2) "FUNC_SS1(", nr0 ,  ")-> " ,"BOT: ", hulp(BOTindices), &
                     " -- SUR: ", hulp(SRFindices)
             endif

          else
#ifdef ONESOURCE
             ! notice that the negative sign is already included in FluxFunctions.F90
             ! (l. 80) thus there is a further change of sign here
             hulp(1:klev)= hulp(1:klev) &
                    - flx_t(i) * D3SOURCE(idummy,flx_states(i),:)
#else
             hulp(1:klev)= hulp(1:klev) &
                  + flx_t(i) * D3SOURCE(flx_states(i),idummy,:)
#endif
             !write(LOGUNIT,*) "FUNC_CORR_SS0: ",nr0, flx_states(i), idummy, "HULP", hulp(:)
             ! correcting for fluxes  to other systems
             if ( flx_states(i) ==idummy) then
                !write(*,*) "Func=", nr0, " Dir=", flx_SS(i), " Compo1=", flx_states(i), " Compo2=", idummy 
                write(LOGUNIT,format1) "FUNC_SS0(", nr0 , ")->" , hulp(BOTindices), &
                     " PELBOTTOM (", flx_states(i), "): ",  PELBOTTOM(flx_states(i),:)
                write(LOGUNIT,format1) "FUNC_SS0(", nr0 , ")->" , hulp(BOTindices), &
                     " PELSURFACE(", flx_states(i), "): ", PELSURFACE(flx_states(i),:)
                hulp(BOTindices)=hulp(BOTindices)-flx_t(i) *max(ZERO,&
                       PELBOTTOM(flx_states(i),:))/Depth(BOTindices)/SEC_PER_DAY
                hulp(SRFindices)=hulp(SRFindices)-flx_t(i) *max(ZERO, &
                    PELSURFACE(flx_states(i),:))/Depth(SRFindices)/SEC_PER_DAY
                write(LOGUNIT,format2) "FUNC_SS0(", nr0 , ")-> " , "BOT: ", hulp(BOTindices), &
                     " -- SUR: ", hulp(SRFindices)
             endif
          endif
          !write(LOGUNIT,*) "FUNC_FIN     : ",nr0, flx_states(i), flx_ostates(i), "HULP", hulp(:)
        enddo
      endif

      hulp(1:klev)=hulp(1:klev)*SEC_PER_DAY

      select case ( flx_option(nr) )
        case(0) !Raw rate
          out(1:klev)=hulp(1:klev);
        case(1) !Specific rate
          out(1:klev)=1.0D-80
          k=0
          ! Sum the total mass
          if ( flx_CalcIn(nr) == iiBen) then
#if defined INCLUDE_BEN || defined INCLUDE_SEAICE
             do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
               if ( k.ne. flx_states(i) ) then
                  k=flx_states(i)
                  out(1:klev)=out(1:klev) + D2STATE(k,:)
               endif
              enddo
#endif
          else 
             do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
               if ( k.ne. flx_states(i) ) then
                 k=flx_states(i)
                 out(1:klev)=out(1:klev) + D3STATE(k,:)
               endif
             enddo
          endif
          ! Compute specific rate
          out(1:klev)=hulp(1:klev)/out(1:klev)
#ifdef BFM_GOTM
        ! This part was only tested with GOTM
        case(2) ! summing the column :perm2
          ! d3 -->d2 var.
          hulp(1:klev)=hulp(1:klev) *Depth(1:klev)
          if (NO_BOXES_X==1 .and. NO_BOXES_Y==1 ) then
            out(1)=sum(hulp(1:klev))
          else
            DO BoxNumberY=1,NO_BOXES_Y
              DO BoxNumberX=1,NO_BOXES_X
                 r=0.0
                 DO BoxNumberZ=1,NO_BOXES_Z
                   BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
                   r=r+hulp(BoxNumber) 
                 enddo
                 BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)
                 out(BoxNumberXY)=r
              enddo
            enddo
          endif
#endif
      end select
#else
      !if D1SOURCE

      !nr=nr0;if ( mode == 2 ) nr=nr0+flx_cal_ben_start
      nr=nr0
      !klev=NO_BOXES ; if ( flx_CalcIn(nr) == iiBen)  klev=NO_BOXES_XY 
      klev=NO_BOXES

      !out(:)  = ZERO
      out(:)  = D3FLUX_FUNC(nr0,:)

      ! correcting for fluxes  to other systems along the diagonal (origin == destination)
      do idx_i=stPelStateS, stPelStateE
         origin      = idx_i
         destination = idx_i

         if( allocated( D3FLUX_MATRIX(origin,destination)%p ) ) then

            do idx_j=1, SIZE(D3FLUX_MATRIX(origin,destination)%p)

               if( ABS(D3FLUX_MATRIX(origin,destination)%p(idx_j)) .eq. nr0 ) then

                  !write(*,*) " FUNC:  ", ABS(D3FLUX_MATRIX(origin,destination)%p(idx_j)), " Compo: ", origin 
                  if( D3FLUX_MATRIX(origin,destination)%dir(idx_j) == 0  ) then ! "A->B" => (out flow) => flux < ZERO => D3SINK
                     write(LOGUNIT,format1) "FUNC_SS0(", nr0 , ")->", out(BOTindices)/SEC_PER_DAY, &
                          " PELBOTTOM (", origin, "): ", PELBOTTOM(origin,:)
                     write(LOGUNIT,format1) "FUNC_SS0(", nr0 , ")->" ,out(BOTindices)/SEC_PER_DAY, &
                          " PELSURFACE(", origin, "): ", PELSURFACE(origin,:)
                     out(BOTindices) = out(BOTindices) - &
                          SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) * &
                          max(ZERO, PELBOTTOM(origin,:)) / Depth(BOTindices)
                     out(SRFindices) = out(SRFindices) - &
                          SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) * &
                          max(ZERO, PELSURFACE(origin,:)) / Depth(SRFindices)
                     write(LOGUNIT,format2) "FUNC_SS0(", nr0 , ")-> ", "BOT: ", out(BOTindices)/SEC_PER_DAY, &
                          " -- SUR: ", out(SRFindices)/SEC_PER_DAY
                  else ! "A<-B" => (in flow) => flux > ZERO => D3SOURCE
                     write(LOGUNIT,format1) "FUNC_SS1(", nr0 , ")->", out(BOTindices)/SEC_PER_DAY, &
                          " PELBOTTOM (", origin, "): ", PELBOTTOM(origin,:)
                     write(LOGUNIT,format1) "FUNC_SS1(", nr0 , ")->" ,out(BOTindices)/SEC_PER_DAY, &
                          " PELSURFACE(", origin, "): ", PELSURFACE(origin,:)
                     out(BOTindices) = out(BOTindices) + &
                          SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) * &
                          min(ZERO, PELBOTTOM(origin,:)) / Depth(BOTindices)
                     out(SRFindices) = out(SRFindices) + &
                          SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) * &
                          min(ZERO, PELSURFACE(origin,:)) / Depth(SRFindices)
                     write(LOGUNIT,format2) "FUNC_SS1(", nr0 , ")-> ", "BOT: ", out(BOTindices)/SEC_PER_DAY, &
                          " -- SUR: ", out(SRFindices)/SEC_PER_DAY
                  endif

               end if

            end do

         end if

      end do



#endif
      write(LOGUNIT,*) "-------MAKE FLUX END--------------"
      return
      end subroutine make_flux_output
!EOC
!-----------------------------------------------------------------------
