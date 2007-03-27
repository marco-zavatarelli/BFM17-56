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
      use constants, only: RLEN, SEC_PER_DAY
      use mem, only: NO_BOXES_XY,NO_BOXES,NO_BOXES_X,NO_BOXES_Y, &
            NO_BOXES_Z,BoxNumberX,BoxNumberY,BoxNumberZ,BoxNumberXY , &
            BoxNumber
      use mem, only: iiBen
      use mem, only: flx_calc_nr,flx_CalcIn,flx_t,flx_states, &
              flx_ostates,flx_SS,flx_cal_ben_start,flx_option
      use mem, only: D2SINK,D3SINK,D2SOURCE,D3SOURCE,D2STATE,D3STATE
      use mem, only: PELBOTTOM,PELSURFACE,Depth


!
! !INPUT PARAMETERS:
      implicit none
      integer,intent(IN)                  ::mode
      integer,intent(IN)                  ::nr0
      integer,intent(IN)                  ::zlev
      integer,intent(IN)                  ::nlev
      real(RLEN),intent(OUT),dimension(zlev:nlev)  :: out

!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij (NIOZ)
!

!
! !LOCAL VARIABLES:
      integer      ::nr
      integer      ::i
      integer      ::k
      integer      ::klev
      real(RLEN),dimension(zlev:NO_BOXES):: hulp
      real(RLEN)   :: r

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! user defined external functions
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      integer, external  :: D3toD1
      integer, external  :: D2toD1

!EOP
!-----------------------------------------------------------------------
!BOC

      nr=nr0;if ( mode == 2 ) nr=nr0+flx_cal_ben_start
      klev=NO_BOXES ; if ( flx_CalcIn(nr) == iiBen)  klev=NO_BOXES_XY 
      hulp(:)=0.0
      if ( flx_CalcIn(nr) == iiBen) then
        do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
          if (flx_SS(i) ==1 ) then
             hulp(1:klev)= hulp(1:klev) &
                    + flx_t(i) * D2SINK(flx_states(i),flx_ostates(i),:)
          else
             hulp(1:klev)= hulp(1:klev) &
                   + flx_t(i) * D2SOURCE(flx_states(i),flx_ostates(i),:)
          endif
        enddo
      else
        do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
          if (flx_SS(i) ==1 ) then
             hulp(1:klev)= hulp(1:klev) &
                    + flx_t(i) * D3SINK(flx_states(i),flx_ostates(i),:)
             ! correcting for fluxes  to other systems
             if ( flx_states(i) ==flx_ostates(i)) then
                hulp(1)=hulp(1)+flx_t(i) *min(0.0D+00,&
                       PELBOTTOM(flx_states(i),1))/Depth(1)/SEC_PER_DAY
                hulp(klev)=hulp(klev)+flx_t(i) *min(0.0D+00,&
                    PELSURFACE(flx_states(i),1))/Depth(klev)/SEC_PER_DAY
             endif

          else
             hulp(1:klev)= hulp(1:klev) &
                  + flx_t(i) * D3SOURCE(flx_states(i),flx_ostates(i),:)
             ! correcting for fluxes  to other systems
             if ( flx_states(i) ==flx_ostates(i)) then
                hulp(1)=hulp(1)-flx_t(i) *max(0.0D+00,&
                       PELBOTTOM(flx_states(i),1))/Depth(1)/SEC_PER_DAY
                hulp(nlev)=hulp(nlev)-flx_t(i) *max(0.0D+00, &
                    PELSURFACE(flx_states(i),1))/Depth(nlev)/SEC_PER_DAY
             endif
          endif
        enddo
      endif

      hulp(1:klev)=hulp(1:klev)*SEC_PER_DAY        

      select case ( flx_option(nr) )
        case(0) !Specific rate
          out(1:klev)=hulp(1:klev);
        case(1) !Specific rate
          out(1:klev)=1.0D-80
          k=0
          if ( flx_CalcIn(nr) == iiBen) then
             do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
               if ( k.ne. flx_states(i) ) then
                  k=flx_states(i)
                  out(1:klev)=out(1:klev) + D2STATE(k,:)
               endif
              enddo
          else 
             do i=flx_calc_nr(nr-1)+1,flx_calc_nr(nr)
               if ( k.ne. flx_states(i) ) then
                 k=flx_states(i)
                 out(1:klev)=out(1:klev) + D3STATE(k,:)
               endif
             enddo
          endif
          out(1:klev)=hulp(1:klev)/out(1:klev)
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
      end select

      return
      end subroutine make_flux_output
!EOC
!-----------------------------------------------------------------------




