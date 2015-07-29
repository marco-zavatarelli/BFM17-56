!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Global helper functions
!
! !INTERFACE:
!  Included in ModuleMem.F90
!
! !DESCRIPTION:
!  This file contains the helper routines for the assignment
!  and retrieval of rates from the main source/sink term arrays
!  (D3SOURCE and D3SINK).
!  It also includes all the routines for NaN checking (only IFORT)
!  This is a special file included when the EXPLICIT_SINK key is activated.
!
! !REVISION HISTORY:
!  Author(s): Piet Ruardij & Marcello Vichi
!
!EOP
!-----------------------------------------------------------------------
!BOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! flux functions
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          subroutine flux_vector(iiSub,origin,destination,flux)
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            use constants, only: RLEN, ZERO,  SEC_PER_DAY, DAY_PER_SEC
            use global_mem, only: LOGUNIT
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Implicit typing is never allowed
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            implicit none
            integer,intent(IN) :: iiSub
            integer :: origin
            integer :: destination
            real(RLEN) :: flux(:)
            integer :: i, j
            character(len=8) :: D23
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !BEGIN compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            if ( destination ==0 ) then
               ! call (iiSub,origin,origin,-flux)
               destination = origin
               flux        = -flux
            elseif ( origin ==0 ) then
              ! call (iiSub,destination,destination,flux)
               origin = destination
            endif

#ifdef IFORT
            TESTNANVECTOR(flux,iiSub,origin,destination)
            CHECKFLUX(-1,iiSub,origin,destination)
#endif

            if ( origin /= destination )  then
#ifdef DEBUG
              if ( minval(flux) < ZERO) then
                do i=1,size(flux)
                  if (flux(i)< 0.0D+00) write(LOGUNIT,'(''at level:'',I4)') i
                enddo
                D23="Pelagic"
                if ( iiSub == iiBen) D23="Benthic"
                write(LOGUNIT,'(''In '',A,'':origin='',i4,'' destination='',i4)') &
                  D23, origin,destination
                write(LOGUNIT,'(''flux='',(G16.8))') flux
                STDERR  "Error in flux_vector function: negative flux !"
                do i=1,size(flux)
                  if (flux(i)< 0.0D+00) then
                    if( iiSub == iiPel ) then
                      write(LOGUNIT,*) "state value origin:",D3STATE(origin,i)
                      write(LOGUNIT,*) "state value destination:",D3STATE(destination,i)
#if defined INCLUDE_SEAICE
                    elseif ( iiSub == iiIce ) then
                      write(LOGUNIT,*) "state value origin:",D2STATE_ICE(origin,i)
                      write(LOGUNIT,*) "state value destination:",D2STATE_ICE(destination,i)
#endif
#if defined INCLUDE_BEN
                    elseif ( iiSub == iiBen ) then
                      write(LOGUNIT,*) "state value origin:",D2STATE_BEN(origin,i)
                      write(LOGUNIT,*) "state value destination:",D2STATE_BEN(destination,i)
#endif
                    endif
                  endif
                enddo
                call BFM_ERROR("flux_vector","negative flux")
              endif ! minval<0
#endif
              select case ( iiSub )
                case (iiPel)
                   D3SINK(origin,destination,:)  =  D3SINK(origin,destination,:) + & 
                        flux*DAY_PER_SEC
                   D3SOURCE(destination,origin,:)=  D3SOURCE(destination,origin,:) + &
                        flux*DAY_PER_SEC
                   if( allocated( D3FLUX_MATRIX ) .AND. &
                        allocated( D3FLUX_MATRIX(origin,destination)%p ) ) then
                      do j=1, SIZE(D3FLUX_MATRIX(origin,destination)%p)
                         D3FLUX_FUNC( ABS(D3FLUX_MATRIX(origin,destination)%p(j)), : ) =      &
                              D3FLUX_FUNC( ABS(D3FLUX_MATRIX(origin,destination)%p(j)), : ) + &
                              (SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(j) ) * flux )
                      end do
                   endif
#if defined INCLUDE_SEAICE
                case (iiIce)
                   D2SINK_ICE(origin,destination,:)  =  D2SINK_ICE(origin,destination,:) + & 
                        flux*DAY_PER_SEC
                   D2SOURCE_ICE(destination,origin,:)=  D2SOURCE_ICE(destination,origin,:) + &
                        flux*DAY_PER_SEC
                   if( allocated( D2FLUX_MATRIX_ICE ) .AND. &
                        allocated( D2FLUX_MATRIX_ICE(origin,destination)%p ) ) then
                      do j=1, SIZE(D2FLUX_MATRIX_ICE(origin,destination)%p)
                         D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) =      &
                              D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) + &
                              (SIGN( 1, D2FLUX_MATRIX_ICE(origin,destination)%p(j) ) * flux(1) )
                      end do
                   end if
#endif
#if defined INCLUDE_BEN
                case (iiBen)
                   D2SINK_BEN(origin,destination,:)  =  D2SINK_BEN(origin,destination,:) + & 
                        flux*DAY_PER_SEC
                   D2SOURCE_BEN(destination,origin,:)=  D2SOURCE_BEN(destination,origin,:) + &
                        flux*DAY_PER_SEC
                   if( allocated( D2FLUX_MATRIX_BEN ) .AND. &
                        allocated( D2FLUX_MATRIX_BEN(origin,destination)%p ) ) then
                      do j=1, SIZE(D2FLUX_MATRIX_BEN(origin,destination)%p)
                         D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) =      &
                              D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) + &
                              (SIGN( 1, D2FLUX_MATRIX_BEN(origin,destination)%p(j) ) * flux(1) )
                      end do
                   end if
#endif
              end select
            else ! origin==destination
              select case ( iiSub )
                case (iiPel)
                  where ( flux > ZERO )
                    D3SOURCE(origin,destination,:) =D3SOURCE(origin,destination,:) + &
                      flux*DAY_PER_SEC
                  elsewhere
                    D3SINK(destination,origin,:) =D3SINK(destination,origin,:) - &
                      flux*DAY_PER_SEC
                  end where
                  if( allocated( D3FLUX_MATRIX ) .AND. &
                       allocated( D3FLUX_MATRIX(origin,destination)%p ) ) then
                     do j=1, SIZE(D3FLUX_MATRIX(origin,destination)%p)
                        if( D3FLUX_MATRIX(origin,destination)%dir(j) == 1 ) then ! "A->B" => (out flow) => flux < ZERO => D3SINK
                           where( flux < ZERO )
                              D3FLUX_FUNC( ABS(D3FLUX_MATRIX(origin,destination)%p(j)), : ) =     &
                                   D3FLUX_FUNC( ABS(D3FLUX_MATRIX(origin,destination)%p(j)), : ) - flux
                           end where
                        else ! "A<-B" => (in flow) => flux > ZERO => D3SOURCE
                           where( flux > ZERO )
                              D3FLUX_FUNC( ABS(D3FLUX_MATRIX(origin,destination)%p(j)), : ) =     &
                                   D3FLUX_FUNC( ABS(D3FLUX_MATRIX(origin,destination)%p(j)), : ) + flux
                           end where
                        end if
                     end do
                  endif
#if defined INCLUDE_SEAICE
                case (iiIce)
                  where ( flux > ZERO )
                    D2SOURCE_ICE(origin,destination,:) =D2SOURCE_ICE(origin,destination,:) + &
                      flux*DAY_PER_SEC
                  elsewhere
                    D2SINK_ICE(destination,origin,:) =D2SINK_ICE(destination,origin,:) - &
                      flux*DAY_PER_SEC
                  end where
                  if( allocated( D2FLUX_MATRIX_ICE ) .AND. &
                       allocated( D2FLUX_MATRIX_ICE(origin,destination)%p ) ) then
                     do j=1, SIZE(D2FLUX_MATRIX_ICE(origin,destination)%p)
                        if( D2FLUX_MATRIX_ICE(origin,destination)%dir(j) == 1 ) then ! "A->B" => (out flow) => flux < ZERO => D2SINK_ICE
                           if( flux(1) < ZERO ) then 
                              D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) =     &
                                   D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) - flux(1)
                           end if
                        else ! "A<-B" => (in flow) => flux > ZERO => D2SOURCE_ICE
                           if( flux(1) > ZERO ) then
                              D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) =     &
                                   D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) + flux(1)
                           end if
                        end if
                     end do
                  endif
#endif
#if defined INCLUDE_BEN
                case (iiBen)
                  where ( flux > ZERO )
                    D2SOURCE_BEN(origin,destination,:) =D2SOURCE_BEN(origin,destination,:) + &
                      flux*DAY_PER_SEC
                  elsewhere
                    D2SINK_BEN(destination,origin,:) =D2SINK_BEN(destination,origin,:) - &
                      flux*DAY_PER_SEC
                  end where
                  if( allocated( D2FLUX_MATRIX_BEN ) .AND. &
                       allocated( D2FLUX_MATRIX_BEN(origin,destination)%p ) ) then
                     do j=1, SIZE(D2FLUX_MATRIX_BEN(origin,destination)%p)
                        if( D2FLUX_MATRIX_BEN(origin,destination)%dir(j) == 1 ) then ! "A->B" => (out flow) => flux < ZERO => D2SINK_BEN
                           if( flux(1) < ZERO ) then 
                              D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) =     &
                                   D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) - flux(1)
                           end if
                        else ! "A<-B" => (in flow) => flux > ZERO => D2SOURCE_BEN
                           if( flux(1) > ZERO ) then
                              D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) =     &
                                   D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) + flux(1)
                           end if
                        end if
                     end do
                  endif
#endif
              end select

            endif !origin <> destination

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !END compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            return
          end subroutine flux_vector
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          subroutine flux(grid_nr,iiSub,origin,destination,flow,error)
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO, SEC_PER_DAY, DAY_PER_SEC
          use global_mem, only: LOGUNIT
          implicit none
          integer,intent(IN)                 :: grid_nr
          integer,intent(IN)                 :: iiSub
          integer,intent(IN)                 :: origin
          integer,intent(IN)                 :: destination
          real(RLEN),intent(IN)              :: flow
          integer,intent(INOUT),optional     :: error
          character(len=8)                   :: D23
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !BEGIN compute
#ifdef IFORT
          TESTNAN(flow,grid_nr,iiSub,origin,destination)
          CHECKFLUX(grid_nr,iiSub,origin,destination)
#endif

          if ( origin /= destination ) then
#ifdef DEBUG
            if ( flow < ZERO) then
              D23="Pelagic"
              if ( iiSub == iiBen) D23="Benthic"
              write(LOGUNIT,'(''In '',A,'':origin='',i4,'' destination='',i4)') &
                D23, origin,destination
              write(LOGUNIT,*) "Error in (scalar) vector  function: negative flux!"
              write(LOGUNIT,*) "origin,destination:", origin,destination
              write(LOGUNIT,*) flow
              if( iiSub == iiPel ) 
                 write(LOGUNIT,*) "state value origin:",D3STATE(origin,grid_nr)
                 write(LOGUNIT,*) "state value destination:",D3STATE(destination,grid_nr)
#if defined INCLUDE_SEAICE
              elseif ( iiSub == iiIce)  then
                 write(LOGUNIT,*) "state value origin:",D2STATE_ICE(origin,grid_nr)
                 write(LOGUNIT,*) "state value destination:",D2STATE_ICE(destination,grid_nr)
#endif
#if defined INCLUDE_BEN
              elseif ( iiSub == iiBen)  then
                 write(LOGUNIT,*) "state value origin:",D2STATE_BEN(origin,grid_nr)
                 write(LOGUNIT,*) "state value destination:",D2STATE_BEN(destination,grid_nr)
#endif
              endif
              STDERR "Error in (scalar)flux function:negative flux !"
              call BFM_ERROR("flux","negative flux")
              if ( present(error)) error=1
            endif ! flow<0
#endif
            select case ( iiSub )
              case (iiPel)
                D3SINK(origin,destination,grid_nr)=flow*DAY_PER_SEC
                D3SOURCE(destination,origin,grid_nr)= flow*DAY_PER_SEC
#if defined INCLUDE_SEAICE
              case (iiIce)
                D2SINK_ICE(origin,destination,grid_nr)= flow*DAY_PER_SEC
                D2SOURCE_ICE(destination,origin,grid_nr)= flow*DAY_PER_SEC
#endif
#if defined INCLUDE_BEN
              case (iiBen)
                D2SINK_BEN(origin,destination,grid_nr)= flow*DAY_PER_SEC
                D2SOURCE_BEN(destination,origin,grid_nr)= flow*DAY_PER_SEC
#endif
            end select
          else
            select case ( iiSub )
              case (iiPel)
                if (flow > ZERO ) then
                  D3SOURCE(destination,origin,grid_nr)=  &
                           D3SOURCE(destination,origin,grid_nr)+flow*DAY_PER_SEC
                else
                  D3SINK(origin,destination,grid_nr)=    &
                         D3SINK(origin,destination,grid_nr)-flow*DAY_PER_SEC
                endif
#if defined INCLUDE_SEAICE
              case (iiIce)
                if (flow > ZERO ) then
                  D2SOURCE_ICE(destination,origin,grid_nr)=  &
                           D2SOURCE_ICE(destination,origin,grid_nr)+flow*DAY_PER_SEC
                else
                  D2SINK_ICE(origin,destination,grid_nr)=    &
                         D2SINK_ICE(origin,destination,grid_nr)-flow*DAY_PER_SEC
                endif
#endif
#if defined INCLUDE_BEN
              case (iiBen)
                if (flow > ZERO ) then
                  D2SOURCE_BEN(destination,origin,grid_nr)=  &
                           D2SOURCE_BEN(destination,origin,grid_nr)+flow*DAY_PER_SEC
                else
                  D2SINK_BEN(origin,destination,grid_nr)=    &
                         D2SINK_BEN(origin,destination,grid_nr)-flow*DAY_PER_SEC
                endif
#endif
            end select
          endif
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !END compute
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          return
          end subroutine flux
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        function Source_D3_vector(iistate)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the pelagic
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO, SEC_PER_DAY
          implicit none
          integer, intent(IN) ::iistate
          real(RLEN) ::Source_D3_vector(size(D3SOURCE,DIM=3))
          integer :: i
          Source_D3_vector = ZERO
          do i = 1, NO_D3_BOX_STATES 
               Source_D3_vector(:)=Source_D3_vector(:) + &
                                    D3SOURCE(iistate,i,:)
               Source_D3_vector(:)=Source_D3_vector(:) - &
                                     D3SINK(iistate,i,:)
          end do
          Source_D3_vector(:)=Source_D3_vector(:)*SEC_PER_DAY
        end function Source_D3_vector
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#if defined INCLUDE_SEAICE
        function Source_D2_vector_ice(iistate)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the seaice
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO, SEC_PER_DAY
          implicit none
          integer, intent(IN) ::iistate
          real(RLEN) ::Source_D2_vector_ice(size(D2SOURCE_ICE,DIM=3))
          integer :: i
          Source_D2_vector_ice = ZERO
          do i = 1,NO_D2_BOX_STATES_ICE 
               Source_D2_vector_ice(:)=Source_D2_vector_ice(:) + &
                                   D2SOURCE_ICE(iistate,i,:) 
               Source_D2_vector_ice(:)=Source_D2_vector_ice(:) - &
                                    D2SINK_ICE(iistate,i,:)
          end do
          Source_D2_vector_ice(:)=Source_D2_vector_ice(:)*SEC_PER_DAY
        end function Source_D2_vector_ice
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#if defined INCLUDE_BEN
        function Source_D2_vector_ben(iistate)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the benthic
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO, SEC_PER_DAY
          implicit none
          integer, intent(IN) ::iistate
          real(RLEN) ::Source_D2_vector_ben(size(D2SOURCE_BEN,DIM=3))
          integer :: i
          Source_D2_vector_ben = ZERO
          do i = 1,NO_D2_BOX_STATES_BEN
               Source_D2_vector_ben(:)=Source_D2_vector_ben(:) + &
                                   D2SOURCE_BEN(iistate,i,:) 
               Source_D2_vector_ben(:)=Source_D2_vector_ben(:) - &
                                    D2SINK_BEN(iistate,i,:)
          end do
          Source_D2_vector_ben(:)=Source_D2_vector_ben(:)*SEC_PER_DAY
        end function Source_D2_vector_ben
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        function source(iiSub,iibox,iistate)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! function to get actual rate of change
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO, SEC_PER_DAY
          implicit none
          real(RLEN)  ::Source
          integer, intent(IN)  ::iiSub
          integer, intent(IN)  ::iibox
          integer, intent(IN)  ::iistate
          integer :: i
          Source = ZERO
          if ( iiSub == iiPel )  then
             do i = 1, NO_D3_BOX_STATES
                Source=Source+D3SOURCE(iistate,i,iibox)
                Source=Source-D3SINK(iistate,i,iibox)
             end do
             Source=Source*SEC_PER_DAY
#if defined INCLUDE_SEAICE
          elseif ( iiSub == iiIce )  then
             do i = 1, NO_D2_BOX_STATES_ICE
                Source=Source+D2SOURCE_ICE(iistate,i,iibox)
                Source=Source-D2SINK_ICE(iistate,i,iibox)
             end do
             Source=Source*SEC_PER_DAY
#endif
#if defined INCLUDE_BEN
          elseif ( iiSub == iiBen )  then
             do i = 1, NO_D2_BOX_STATES_BEN
                Source=Source+D2SOURCE_BEN(iistate,i,iibox)
                Source=Source-D2SINK_BEN(iistate,i,iibox)
             end do
             Source=Source*SEC_PER_DAY
#endif
          endif
        end function source
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        subroutine unicflux(grid_nr,iiSub,origin,destination)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO
          use global_mem, only: LOGUNIT
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          implicit none
          integer,intent(IN)    :: grid_nr
          integer,intent(IN)    :: origin
          integer,intent(IN)    :: iiSub
          integer,intent(IN)    :: destination
          real(RLEN) :: tot
          character(len=20):: type
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !BEGIN compute
          select case ( iiSub )
            case (iiPel)
              type="D3"
              if ( grid_nr <=0  ) then
                tot=sum(D3SINK(origin,destination,:))
              else
                tot=D3SINK(origin,destination,grid_nr)
              endif
#if defined INCLUDE_SEAICE
            case (iiIce)
              type="D2"
              if ( grid_nr <=0  ) then
                tot=sum(D2SINK_ICE(origin,destination,:))
              else
                tot=D2SINK_ICE(origin,destination,grid_nr)
              endif
#endif
#if defined INCLUDE_BEN
            case (iiBen)
              type="D2"
              if ( grid_nr <=0  ) then
                tot=sum(D2SINK_BEN(origin,destination,:))
              else
                tot=D2SINK_BEN(origin,destination,grid_nr)
              endif
#endif
            case (iiReset)
              D3SINK(:,:,:)=0.0D+00
#if defined INCLUDE_SEAICE
              D2SINK_ICE(:,:,:)=0.0D+00
#endif
#if defined INCLUDE_BEN
              D2SINK_BEN(:,:,:)=0.0D+00
#endif
              return
          end select
          if ( tot > 0.0D+00  ) then
            write(LOGUNIT,'(''Double definition '',A2,''-flux'')')type
            write(LOGUNIT,'(''origin:'',I3,'' destination:'',I3)') origin, destination
            if ( origin /= destination ) then
              STDERR 'double definition of fluxes'
              stop 1006
            endif
          endif
        !END compute
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        return
      end subroutine unicflux
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      subroutine fixed_quota_flux_vector(mode,iiSub,which,origin, &
                                          destination,flux,collect)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      use global_mem, only: LOGUNIT
      use constants, only: RLEN
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Implicit typing is never allowed
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      implicit none
      integer,intent(IN) :: mode
      integer,intent(IN) :: iiSub
      integer,intent(IN) :: which
      integer,intent(IN) :: origin
      integer,intent(IN) :: destination
      real(RLEN),intent(IN),dimension(:) :: flux
      real(RLEN),intent(INOUT),dimension(:) :: collect
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !BEGIN compute
      if ( origin> 0 .and.destination >0) then
         call flux_vector(iiSub,origin, destination,flux)
      else if ( origin > 0 ) then
         call flux_vector(iiSub,origin, origin,-flux)
      elseif ( destination > 0 ) then
         call flux_vector(iiSub,destination, destination,flux)
      elseif (iiSub < 0 ) then
         if ( mode==0)  return
         if ( sum(abs(flux)/(1.0D-80+abs(collect))-1.0D+00)> 1.0D-6) then
              if ( iiSub==-iiN) then
                write(LOGUNIT,'(''Warning: N:C quotum not fixed'')')
              elseif (iiSub==-iiP) then
                write(LOGUNIT,'(''Warning: P:C quotum not fixed'')')
              endif
              return
         endif
      endif      
      if ( mode==1 ) then
        if ( (which == origin) .and.(origin.ne.destination)) then
           collect=collect-flux
        else
           collect=collect+flux
        endif
      endif
      !END compute
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      end subroutine fixed_quota_flux_vector
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! NaN-check routines
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#ifdef IFORT
! Important note: the intrinsic isnan is not allowed with strict F95
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          subroutine testnan_vector(array,iiSub,origin,destination)
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use global_mem, only: LOGUNIT
          implicit none
            real(RLEN),intent(IN)    :: array(:)
            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            integer:: i=0
            do i=1,size(array)
              if (isnan(array(i))) then
                write(LOGUNIT,'(''at level:'',I4)') i
                write(LOGUNIT,'(''origin='',i4,'' destination='',i4)') &
                  origin,destination
                if( iiSub == iiPel ) then
                    write(LOGUNIT,*) "state value origin:",D3STATE(origin,i)
                    write(LOGUNIT,*) "state value destination:",D3STATE(destination,i)
#if defined INCLUDE_SEAICE
                elseif ( iiSub == iiIce ) then
                    write(LOGUNIT,*) "state value origin:",D2STATE_ICE(origin,i)
                    write(LOGUNIT,*) "state value destination:",D2STATE_ICE(destination,i)
#endif
#if defined INCLUDE_BEN
                elseif ( iiSub == iiBen ) then
                    write(LOGUNIT,*) "state value origin:",D2STATE_BEN(origin,i)
                    write(LOGUNIT,*) "state value destination:",D2STATE_BEN(destination,i)
#endif
                endif
                STDERR 'Nan value in flux'
                stop 1002
              endif
            enddo
          end subroutine testnan_vector
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          subroutine testnan(scalar,grid_nr,iiSub,origin,destination)
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use global_mem, only: LOGUNIT
          implicit none
            real(RLEN),intent(IN)    :: scalar
            integer,intent(IN) :: grid_nr
            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            if (isnan(scalar)) then
               write(LOGUNIT,*) 'Nan value in scalar flux'
               write(LOGUNIT,'(''origin='',i4,'' destination='',i4)') origin,destination
                 if ( iiSus == iiPel ) then
                    write(LOGUNIT,*) "Pelagic state value origin:",D3STATE(origin,grid_nr)
                    write(LOGUNIT,*) "Pelagic state value destination:",D3STATE(destination,grid_nr)
#if defined INCLUDE_SEAICE
               elseif ( iiSub == iiIce )  then
                    write(LOGUNIT,*) "Benthic state value origin:",D2STATE_ICE(origin,grid_nr)
                    write(LOGUNIT,*) "Benthic state value destination:",D2STATE_ICE(destination,grid_nr)
#endif
#if defined INCLUDE_BEN
               elseif ( iiSub == iiBen )  then
                    write(LOGUNIT,*) "Benthic state value origin:",D2STATE_BEN(origin,grid_nr)
                    write(LOGUNIT,*) "Benthic state value destination:",D2STATE_BEN(destination,grid_nr)
#endif
               endif
               stop "subroutine TESTNAN forced STOP"
            endif
          end subroutine testnan
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#endif

!EOC
