!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  MODEL POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!  ROUTINE: get_rst
!
! DESCRIPTION
!
! This subroutine opens and reads the fort file containing the restart
! The file path is specified in pom_input.nml
! 
! INTERFACE
!
      subroutine get_rst
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
      use POM, ONLY: time0,u,ub,v,vb,t,tb,s,sb,q2,q2b,q2l,q2lb,kh,km,kq,l, &
                     wubot,wvbot,rho
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
!     -----NAMELIST READING UNIT-----
!
      integer,parameter  :: namlst=10
!
!      -----RECORD LENGTH-----
!
      INTEGER :: RLENGTH
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
     !   Open and read restart file 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
       open(70,file=read_restart,form='unformatted',status='old')
!
       read(70) time0,            &
                u,ub,v,vb,        &
                t,tb,s,sb,        &
                q2,q2b,q2l,q2lb,  &
                kh,km,kq,         &
                l,                &
                wubot,wvbot,      &
                rho
!
      return
!
100   call error_msg_prn(NML_OPEN,"get_rst.F90","pom_input.nml")
102   call error_msg_prn(NML_READ,"get_rst.F90","pom_input")
!
      end subroutine get_rst
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
