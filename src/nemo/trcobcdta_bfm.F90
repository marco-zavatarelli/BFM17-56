#include"cppdefs.h" 

#ifdef key_obc
   SUBROUTINE trcobcdta_bfm( kt )
      !!--------------------------------------------------------------------
      !!              ***  SUBROUTINE trcobc_dta  ***
      !!                   
      !! ** Purpose :   Find the climatological boundary arrays for the specified date,
      !!   The boundary arrays are netcdf files. Three possible cases:
      !!   - one time frame only in the file (time dimension = 1).
      !!     in that case the boundary data does not change in time.
      !!   - many time frames. In that case,  if we have 12 frames
      !!     we assume monthly fields.
      !!     Else, we assume that time_counter is in seconds
      !!     since the beginning of either the current year or a reference
      !!     year given in the namelist.
      !!     (no check is done so far but one would have to check the "unit"
      !!     attribute of variable time_counter).
      !!
      !! History :
      !!        !  98-05 (J.M. Molines) Original code
      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
      !!   9.0  !  04-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
      !!--------------------------------------------------------------------

! MODULES USED


   !!------------------------------------------------------------------------------
   !!   'key_obc'         :                                Open Boundary Conditions
   !!------------------------------------------------------------------------------
   !!   obc_dta           : read u, v, t, s data along each open boundary
   !!   obc_dta_psi       : read psi data along each open boundary (rigid lid only)
   !!------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE daymod          ! calendar
   USE in_out_manager  ! I/O logical units
   USE lib_mpp         ! distributed memory computing
   USE iom
!GELSOMINA
   USE trc             ! ocean passive tracers variables
   USE trcobc_oce      ! tracers open boundary conditions
   USE trcobcdta
   USE trcobcini
!GELSOMINA
      
   use mem, only: D3STATE,NO_D3_BOX_STATES, &
                  NO_BOXES, D3STATEOBC
   use global_mem, only:LOGUNIT,LOGUNITOBC
   use api_bfm
   ! NEMO
!   use oce_trc          ! ocean dynamics and active tracers variables

!

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
!      INTEGER, INTENT( out ) ::   m          ! ocean time-step index
      !! * Local declarations
      INTEGER ::   ji, jj, jk, ii, ij   ! dummy loop indices
      INTEGER ::   itimo, iman, imois
      INTEGER ::   i15
      INTEGER ::   n,m
      REAL(wp) ::   zxy
      !! * Ajouts FD
      INTEGER ::  isrel              ! number of seconds since 1/1/1992
      INTEGER, DIMENSION(1) ::  itobce, itobcw,  & ! number of time steps in OBC files
                                itobcs, itobcn     !    "       "       "       "
      INTEGER ::  istop        
      INTEGER ::  iprint                              ! frequency for printouts.
      INTEGER ::  idvar, id_e, id_w, id_n, id_s       ! file identifiers
      LOGICAL :: llnot_done
      CHARACTER(LEN=25) :: cl_vname
      !!--------------------------------------------------------------------
      INTEGER ::   &
      nlecto,   &  ! switch for the first read
      ntobc1,   &  ! first record used
      ntobc2,   &  ! second record used
      ntobc        ! number of time steps in OBC files 
!      REAL(wp), DIMENSION(:), ALLOCATABLE :: tcobc      ! time_counter variable of BCs
      INTEGER :: AllocStatus,DeallocStatus
      INTEGER,SAVE :: first=0

      ! 1.   First call: check time frames available in files.
      ! -------------------------------------------------------

        WRITE(LOGUNIT,*)  'prima di kt=nit000 trcobcdta_bfm : find boundary data'
        WRITE(LOGUNIT,*)  'kt=',kt
        call flush (LOGUNIT)

      IF ( first == 0 ) THEN
           first=1
         IF( lp_obc_north)   THEN

          allocate(trndta(jpind:jpinf,jpk,jptobc,NO_D3_BOX_STATES),stat=AllocStatus)
          if (AllocStatus  /= 0) stop "error allocating trndta"

            ! initialisation to zero
            trndta(:,:,:,:) = 0.e0
          ENDIF

          IF( lp_obc_south)   THEN
 
          allocate(trsdta(jpisd:jpisf,jpk,jptobc,NO_D3_BOX_STATES),stat=AllocStatus)
          if (AllocStatus  /= 0) stop "error allocating trsdta"
   
            ! initialisation to zero
            trsdta(:,:,:,:) = 0.e0
          ENDIF
  
          IF( lp_obc_east)   THEN

          allocate(tredta(jpjed:jpjef,jpk,jptobc,NO_D3_BOX_STATES),stat=AllocStatus)
          if (AllocStatus  /= 0) stop "error allocating tredta"
  
            ! initialisation to zero
            tredta(:,:,:,:) = 0.e0
          ENDIF 

          IF( lp_obc_west)   THEN  
 
          allocate(trwdta(jpjwd:jpjwf,jpk,jptobc,NO_D3_BOX_STATES),stat=AllocStatus)
          if (AllocStatus  /= 0) stop "error allocating trwdta"

            ! initialisation to zero
            trwdta(:,:,:,:) = 0.e0
          ENDIF

         ENDIF

         nlecto = 0
         itobce(1) = 0   ;    itobcw(1) = 0
         itobcn(1) = 0   ;    itobcs(1) = 0

        WRITE(LOGUNIT,*)  'trcobcdta_bfm : find boundary data'
        WRITE(LOGUNIT,*)      '~~~~~~~'       
        call flush (LOGUNIT)

         IF ( ntrcobc_dta == 0 ) THEN
           WRITE(LOGUNIT,*)  'Passive Tracers OBC data taken from initial conditions in BFM.'
            call flush (LOGUNIT)
            ntobc1 = 1
            ntobc2 = 1
         ELSE
            CALL ctl_stop( 'trcobcdta_bfm: Time-varying Biological OBC data not ready yet!')
         ENDIF

       ! 1.2  Data set to zero
       !                        or initial conditions if ntrcobc_dta == 0
       ! --------------------------------------------------------------


           
      IF ( kt == nit000 ) THEN

       DO m = 1,NO_D3_BOX_STATES

         IF (D3STATEOBC(m)/= 0.) THEN
!        WRITE(LOGUNIT,*)  'controllo di m per D3STATESOBC=', m 
!        call flush (LOGUNIT)
            ! remap the biological states to the OBC array
#ifdef USEPACK
            trn(:,:,:,1) = unpack(D3STATE(m,:),SEAmask,ZEROS)
#else
            DO n=1,NO_BOXES
               trn(iwet(n),jwet(n),kwet(n),1)=D3STATE(m,n)
            ENDDO
#endif
            
         IF( lp_obc_east )   THEN


            !                                    ! ================== !
            IF( ntrcobc_dta == 0 )   THEN        ! initial state used
                                                 ! ================== !
               !  Fills trcedta (global arrays)
               !  Remark: this works for njzoom = 1.
               !          Should the definition of ij include njzoom?
                DO ji = nie0, nie1
                  DO jk = 1, jpkm1
                     DO jj = nje0p1, nje1m1
                        ij = jj -1 + njmpp
                        tredta(ij,jk,1,m) = trn(ji,jj,jk,1)*tmask(ji,jj,jk)
                     END DO
                  END DO
                END DO
            ENDIF
         ENDIF

         IF(lp_obc_west )   THEN


            !                                    ! ================== !
            IF( ntrcobc_dta == 0 )   THEN           ! initial state used !
               !                                 ! ================== !
               !  Fills swdta, twdta, uwdta (global arrays)
               !  Remark: this works for njzoom = 1. 
               !          Should the definition of ij include njzoom?
               DO ji = niw0, niw1
                  DO jk = 1, jpkm1
                     DO jj = njw0p1, njw1m1
                        ij = jj -1 + njmpp
                        trwdta(ij,jk,1,m) = trn(ji,jj,jk,1)*tmask(ji,jj,jk)

                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_north)   THEN



            !                                    ! ================== !
            IF( ntrcobc_dta == 0 )   THEN           ! initial state used
               !                                 ! ================== !
               !  Fills sndta, tndta, vndta (global arrays)
               !  Remark: this works for njzoom = 1. 
               !          Should the definition of ij include njzoom?
               DO jj = njn0, njn1
                  DO jk = 1, jpkm1
                     DO ji = nin0p1, nin1m1
                        ii = ji -1 + nimpp
                        trndta(ii,jk,1,m) = trn(ji,jj,jk,1)*tmask(ji,jj,jk)
                      IF(m.eq.1.and. jk.eq.1) THEN
                            WRITE(LOGUNITOBC,*) 'VARIABLE m=',m
                            WRITE(LOGUNITOBC,*) 'ii=' ,ii,'jj=',jj,'ji=',ji
                            WRITE(LOGUNITOBC,*) 'trndta(ii,jk,1,m)=', trndta(ii,jk,1,m)
                            CALL flush(LOGUNITOBC)
                      ENDIF
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_south )   THEN
         

            !                                    ! ================== !
            IF( ntrcobc_dta == 0 )   THEN           ! initial state used
               !                                 ! ================== !
               !  Fills ssdta, tsdta, vsdta (global arrays)
               !  Remark: this works for njzoom = 1. 
               !          Should the definition of ij include njzoom?
               DO jj = njs0, njs1
                  DO jk = 1, jpkm1
                     DO ji = nis0p1, nis1m1
                        ii = ji -1 + nimpp
                        trsdta(ii,jk,1,m) = trn(ji,jj,jk,1)*tmask(ji,jj,jk)
                          IF (jk.eq.1.and.m.eq.1) then 
                            WRITE(LOGUNIT,*) 'trsdta######'
                            WRITE(LOGUNIT,*) 'm=',m,'ii=' ,ii,'jk=',jk
                            WRITE(LOGUNIT,*) 'trsdta(ii,jk,1,m)=', trsdta(ii,jk,1,m)
                            CALL flush(LOGUNIT)
                          ENDIF
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF
         ENDIF
      END DO   ! loop over NO_D3_STATES

      ENDIF        !       end if kt == nit000
 
   !  IF (kt==nit000+1) then DEALLOCATION
     
      
      ! 2.  Initialize the time we are at.
      !     Does this every time the routine is called,
      !     excepted when ntrobc_dta = 0
      !---------------------------------------------------------------------
!      IF( ntrcobc_dta == 0 )   THEN
!         itimo = 1
!         zxy   = 0
!      ELSE
!         IF( ntobc == 1 )   THEN
!            itimo = 1
!         ELSE IF( ntobc == 12 )   THEN      !   BC are monthly   
!!            ! we assume we have climatology in that case
!            iman  = 12
!            i15   = nday / 16
!            imois = nmonth + i15 - 1
!            IF( imois == 0 )   imois = iman
!            itimo = imois   
!         ELSE
!            WRITE(LOGUNIT,*) 'data other than constant or monthly', kt
!            call flush (LOGUNIT)
!            iman  = ntobc
!            itimo = FLOOR( kt*rdt / (tcobc(2)-tcobc(1)) )
!            isrel = kt*rdt
!         ENDIF
!      ENDIF
      
      ! 2.1 Read two records in the file if necessary
      ! ---------------------------------------------

      
      ! 3.  Call at every time step :
      !     Linear interpolation of BCs to current time step 
      ! ----------------------------------------------------

!      IF( ntobc == 1 .OR. ntrcobc_dta == 0 )   THEN 
!         zxy = 0.
!      ELSE IF( ntobc == 12 )   THEN         
!         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!      ELSE
!         zxy = (tcobc(ntobc1)-FLOAT(isrel))/(tcobc(ntobc1)-tcobc(ntobc2))
!      ENDIF
!Per evitare la duplicazione degli array trfoe-w-n-s l'interpolazione viene fatta in trcobc.F90      
      
        WRITE(LOGUNIT,*)  'exit trcobcdta_bfm'
        call flush (LOGUNIT)

   END SUBROUTINE trcobcdta_bfm
#else
   SUBROUTINE trcobcdta_bfm( kt )
   ! empty routine when OBC are not defined
   END SUBROUTINE trcobcdta_bfm
#endif
