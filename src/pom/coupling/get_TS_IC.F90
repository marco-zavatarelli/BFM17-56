!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: get_TS_IC
!
!DESCRIPTION    
!
! This subroutine opens  and read files containing the T&S initial conditions
! Files are read in direct access mode reading path
! specified in pom_input nml
!
! !INTERFACE
!
     subroutine get_TS_IC
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use global_mem, ONLY: error_msg_prn, NML_OPEN, NML_READ
!
      use Service, only: wind_input,surfaceS_input,radiance_input,&
                         ism_input,Sal_input,Temp_input,W_input,Weddy_input1,&
                         Weddy_input2,Sprofile_input,Tprofile_input,heat_input,&
                         surfNut_input,bottNut_input,ISM,read_restart
!    
      use pom, ONLY: KB,T,TB,S,SB

!-------------------------------------------------------------------------!
!
!BOC
!
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Implicit typing is never allowed
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      IMPLICIT NONE 
!
!     -----LOOP COUNTER-----
!
      INTEGER :: K
!
!      -----RECORD LENGTH-----
!
      INTEGER :: RLENGTH
!
!     -----NAMELIST READING UNIT-----
!
      integer,parameter  :: namlst=10
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Open nml with forcing data path specified and read data
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

       namelist /pom_input/ wind_input,surfaceS_input,radiance_input, &
                          ism_input,Sal_input,Temp_input,W_input,Weddy_input1,&
                          Weddy_input2,Sprofile_input,Tprofile_input,heat_input,&
                          surfNut_input,bottNut_input,read_restart
!
       open(namlst,file='pom_input.nml',status='old',action='read',err=100)
       read(namlst,nml=pom_input, err=102)
       rewind(namlst)
       close(namlst)
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !   OPEN SALINITY INITIAL CONDITION FILE 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
       inquire(IOLENGTH=rlength) SB(1)
       open(29,file=Sprofile_input,form='unformatted',access='direct',recl=rlength)
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !   OPEN TEMPERATURE INITIAL CONDITION FILE 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
       inquire(IOLENGTH=rlength) TB(1)
       open(10,file=Tprofile_input,form='unformatted',access='direct',recl=rlength)
!
!
!    -----READ T&S INITIAL CONDITIONS-----
!
     DO K = 1,KB
!
           READ (29,REC=K) SB(K)
           READ (10,REC=K) TB(K)
!
     END DO
!
     T(:)=TB(:)
     S(:)=SB(:)
!
     return
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Define error messages
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

100   call error_msg_prn(NML_OPEN,"get_TS_IC.F90","pom_input.nml")
102   call error_msg_prn(NML_READ,"get_TS_IC.F90","pom_input")
!
      end subroutine get_TS_IC
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
