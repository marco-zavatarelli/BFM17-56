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
!
! !REVISION HISTORY:
!  Author(s): Marcello Vichi
!
!EOP
!-----------------------------------------------------------------------
!BOC

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! flux functions
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          subroutine flux_vector(iiSub,iiorigin,iidestination,flux)
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            use constants, only: RLEN, ZERO,  SEC_PER_DAY, DAY_PER_SEC
            use global_mem, only: LOGUNIT

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Implicit typing is never allowed
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            implicit none
            integer,intent(IN)       :: iiSub, iiorigin, iidestination
            real(RLEN),dimension(:)  :: flux !is going to be restored at the end if modified

            integer          :: origin, destination
            logical          :: fluxsign
            character(len=8) :: D23
            integer          :: j
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !BEGIN compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            origin      = iiorigin
            destination = iidestination
            fluxsign    = .FALSE. 

            if ( destination == 0 ) then
               ! call (iiSub,origin,origin,-flux)
               destination = origin
               flux        = -flux
               fluxsign    = .TRUE.
            elseif ( origin == 0 ) then
              ! call (iiSub,destination,destination,flux)
               origin = destination
            endif

#ifdef IFORT
            TESTNANVECTOR(flux,iiSub,origin,destination)
            CHECKFLUX(-1,iiSub,origin,destination)
#endif
            if ( origin /= destination ) then
              select case ( iiSub )
                case (iiPel)
                  D3SOURCE(origin,:) = D3SOURCE(origin,:) - &
                      flux*DAY_PER_SEC
                  D3SOURCE(destination,:) = D3SOURCE(destination,:) + &
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
                  D2SOURCE_ICE(origin,:) = D2SOURCE_ICE(origin,:) - &
                      flux*DAY_PER_SEC
                  D2SOURCE_ICE(destination,:) = D2SOURCE_ICE(destination,:) + &
                      flux*DAY_PER_SEC
                  if( allocated( D2FLUX_MATRIX_ICE ) .AND. &
                       allocated( D2FLUX_MATRIX_ICE(origin,destination)%p ) ) then
                     do j=1, SIZE(D2FLUX_MATRIX_ICE(origin,destination)%p)
                        D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) =      &
                             D2FLUX_FUNC_ICE( ABS(D2FLUX_MATRIX_ICE(origin,destination)%p(j)) ) + &
                             (SIGN( 1, D2FLUX_MATRIX_ICE(origin,destination)%p(j) ) * flux (1) )
                     end do
                  endif
#endif
#if defined INCLUDE_BEN
                case (iiBen)
                  D2SOURCE_BEN(origin,:) = D2SOURCE_BEN(origin,:) - &
                      flux*DAY_PER_SEC
                  D2SOURCE_BEN(destination,:) = D2SOURCE_BEN(destination,:) + &
                      flux*DAY_PER_SEC
                  if( allocated( D2FLUX_MATRIX_BEN ) .AND. &
                       allocated( D2FLUX_MATRIX_BEN(origin,destination)%p ) ) then
                     do j=1, SIZE(D2FLUX_MATRIX_BEN(origin,destination)%p)
                        D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) =      &
                             D2FLUX_FUNC_BEN( ABS(D2FLUX_MATRIX_BEN(origin,destination)%p(j)) ) + &
                             (SIGN( 1, D2FLUX_MATRIX_BEN(origin,destination)%p(j) ) * flux (1) )
                     end do
                  endif
#endif
              end select
            else ! origin==destination
              ! In this case the flux carries the proper sign
              select case ( iiSub )
                case (iiPel)
                  D3SOURCE(destination,:) = D3SOURCE(destination,:) + &
                      flux*DAY_PER_SEC
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
                  D2SOURCE_ICE(destination,:) = D2SOURCE_ICE(destination,:) + &
                      flux*DAY_PER_SEC
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
                  D2SOURCE_BEN(destination,:) = D2SOURCE_BEN(destination,:) + &
                      flux*DAY_PER_SEC
                  if( allocated( D2FLUX_MATRIX_BEN ) .AND. &
                       allocated(D2FLUX_MATRIX_BEN(origin,destination)%p ) ) then
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
            end if !origin <> destination

            ! if destination == 0 => need to restore original input sign 
            if ( fluxsign ) flux = -flux 
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !END compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            return
          end subroutine flux_vector
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          subroutine flux(grid_nr,iiSub,origin,destination,flow)
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO, SEC_PER_DAY, DAY_PER_SEC
          implicit none
          integer,intent(IN)                 :: grid_nr
          integer,intent(IN)                 :: iiSub
          integer,intent(IN)                 :: origin
          integer,intent(IN)                 :: destination
          real(RLEN),intent(IN)              :: flow
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !BEGIN compute
#ifdef IFORT
          TESTNAN(flow,grid_nr,iiSub,origin,destination)
          CHECKFLUX(grid_nr,iiSub,origin,destination)
#endif
          if ( origin /= destination ) then
            select case ( iiSub )
              case (iiPel)
                D3SOURCE(origin,grid_nr) =    &
                         D3SOURCE(origin,grid_nr)-flow*DAY_PER_SEC
                D3SOURCE(destination,grid_nr) =    &
                         D3SOURCE(destination,grid_nr)+flow*DAY_PER_SEC
#if defined INCLUDE_SEAICE
              case (iiIce)
                D2SOURCE_ICE(origin,grid_nr) =    &
                         D2SOURCE_ICE(origin,grid_nr)-flow*DAY_PER_SEC
                D2SOURCE_ICE(destination,grid_nr) =    &
                         D2SOURCE_ICE(destination,grid_nr)+flow*DAY_PER_SEC
#endif
#if defined INCLUDE_BEN
              case (iiBen)
                D2SOURCE_BEN(origin,grid_nr) =    &
                         D2SOURCE_BEN(origin,grid_nr)-flow*DAY_PER_SEC
                D2SOURCE_BEN(destination,grid_nr) =    &
                         D2SOURCE_BEN(destination,grid_nr)+flow*DAY_PER_SEC
#endif
            end select
          else
            ! In this case the flux carries the proper sign
            select case ( iiSub )
              case (iiPel)
                 D3SOURCE(destination,grid_nr) = D3SOURCE(destination,grid_nr) + &
                      flow*DAY_PER_SEC
#if defined INCLUDE_SEAICE
              case (iiIce)
                 D2SOURCE_ICE(destination,grid_nr) = D2SOURCE_ICE(destination,grid_nr) + &
                      flow*DAY_PER_SEC
#endif
#if defined INCLUDE_BEN
              case (iiBen)
                 D2SOURCE_BEN(destination,grid_nr) = D2SOURCE_BEN(destination,grid_nr) + &
                      flow*DAY_PER_SEC
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
          real(RLEN) ::Source_D3_vector(size(D3SOURCE,DIM=2))

          Source_D3_vector = ZERO
          Source_D3_vector(:) = D3SOURCE(iistate,:)*SEC_PER_DAY
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
          real(RLEN) ::Source_D2_vector_ice(size(D2SOURCE_ICE,DIM=2))

          Source_D2_vector_ice = ZERO
          Source_D2_vector_ice = D2SOURCE_ICE(iistate,:)*SEC_PER_DAY 
        end function Source_D2_vector_ice
#endif

#if defined INCLUDE_BEN
        function Source_D2_vector_ben(iistate)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the benthic
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          use constants, only: RLEN, ZERO, SEC_PER_DAY
          implicit none
          integer, intent(IN) ::iistate
          real(RLEN) ::Source_D2_vector_ben(size(D2SOURCE_BEN,DIM=2))

          Source_D2_vector_ben = ZERO
          Source_D2_vector_ben = D2SOURCE_BEN(iistate,:)*SEC_PER_DAY 
        end function Source_D2_vector_ben
#endif
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

          Source = ZERO
          if ( iiSub == iiPel )  then
             Source=D3SOURCE(iistate,iibox)*SEC_PER_DAY
#if defined INCLUDE_SEAICE
          elseif ( iiSub == iiIce )  then
             Source=D2SOURCE_ICE(iistate,iibox)*SEC_PER_DAY
#endif
#if defined INCLUDE_BEN
          elseif ( iiSub == iiBen )  then
             Source=D2SOURCE_BEN(iistate,iibox)*SEC_PER_DAY
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
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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
      integer,intent(IN) :: mode                       ! check fix quota mode
      integer,intent(IN) :: iiSub
      integer,intent(IN) :: which                      ! pointer to Zoo carbon
      integer,intent(IN) :: origin
      integer,intent(IN) :: destination
      real(RLEN),intent(IN),dimension(:) :: flux
      real(RLEN),intent(INOUT),dimension(:) :: collect ! save fluxes  
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !BEGIN compute
      if (iiSub >= 0 ) then
          call flux_vector(iiSub,origin, destination,flux)
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
                 if( iiSub == iiPel ) then
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
