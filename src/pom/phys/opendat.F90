!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: Opendat
!
!DESCRIPTION    
!
! This subroutine opens forcing files in direct access mode reading path
! specified in pom_input nml
!
! !INTERFACE
    subroutine opendat(WSU1,WSV1,SSS1,SLUX1,ISM1,SCLIM1,TCLIM1,WCLIM1,&
         WEDDY1,WEDDY2,SB,TB,SWRAD1,WTSURF1,QCORR1,NO3_s1,NH4_s1,PO4_s1,&
         SIO4_s1,O2_b1,NO3_b1,PO4_b1,PON_b1)
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      use global_mem,ONLY:RLEN,LOGUNIT,error_msg_prn,NML_OPEN,NML_READ
      use Service, only: wind_input,surfaceS_input,radiance_input,&
                         ism_input,Sal_input,Temp_input,W_input,Weddy_input1,&
                         Weddy_input2,Sprofile_input,Tprofile_input,heat_input,&
                         surfNut_input,bottNut_input,ISM,read_restart
!    
      use pom, ONLY: KB
!-------------------------------------------------------------------------!
!BOC
!
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     ! Implicit typing is never allowed
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      IMPLICIT NONE 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Scalar Arguments
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      INTEGER :: RLENGTH
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Local Scalars
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      integer,parameter  :: namlst=10
      integer,parameter  :: N_COMP=KB-1
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Array Arguments 
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      real(RLEN) :: WSU1,WSV1,SSS1,SLUX1,SWRAD1,WTSURF1,QCORR1, &
                    NO3_s1,NH4_s1,PO4_s1,SIO4_s1,O2_b1,NO3_b1,PO4_b1,PON_b1
      REAL(RLEN),dimension(KB) :: SCLIM1,TCLIM1,WCLIM1,SB,TB,WEDDY1,WEDDY2
      REAL(RLEN),dimension(N_COMP) :: ISM1

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
       close(namlst)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Wind speed (u,v)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) WSU1,WSV1
       open(11,file=wind_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 11 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Surface salinity
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) SSS1
       open(13,file=surfaceS_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 13 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Radiance
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) SLUX1
       open(17,file=radiance_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 17 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Inorganic suspended matter
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) ISM1(1)
       open(19, file=ism_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 19 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Salinity climatology (diagnostic mode)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) SCLIM1(1)
       open(20, file=Sal_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 20 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Temperature climatology (diagnostic mode)
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) TCLIM1(1)
       open(15, file=Temp_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 15 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  General Circulation W Velocity climatology
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) WCLIM1(1)
       open(215, file=w_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 215 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Intermittant Eddy W Velocity 1
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) WEDDY1(1)
       open(216, file=weddy_input1, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 216 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Intermittant Eddy W Velocity 2
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) WEDDY2(1)
       open(217, file=weddy_input2, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 217 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Salinity initial profile
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) SB(1)
       open(29,file=Sprofile_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 29 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Temperature initial profile
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) TB(1)
       open(10,file=Tprofile_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 10 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Heat flux
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) SWRAD1,WTSURF1,QCORR1
       open(21, file=heat_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 21 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Surface nutrients
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) NO3_s1,NH4_s1,PO4_s1,SIO4_s1
       open(18, file=surfNut_input, form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 18 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Bottom nutrients
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) O2_b1,NO3_b1,PO4_b1,PON_b1
       open(300, file=bottNut_input, form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 300 done'

       write(6,*) 'open units done'
!
     return
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Define error messages
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

100   call error_msg_prn(NML_OPEN,"opendat.F90","pom_input.nml")
102   call error_msg_prn(NML_READ,"opendat.F90","pom_input")

      end subroutine opendat
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
