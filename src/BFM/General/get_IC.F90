!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM17 - 17 Eq. Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: get_IC
!
!DESCRIPTION    
!
! This subroutine opens and read files containing the initial conditions
! for P2c, Z4c, R1c, R6c, N1p, N3n, N4n, and O2o
! Files are read in direct access mode reading path
! specified in bfm_IC_input nml
!
! !INTERFACE
!
     subroutine get_IC
! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      use Service, only: phyto_input,zoop_input,poc_input,&
                         doc_input,phos_input,nit_input,am_input,&
                         oxy_input
      use mem, ONLY: R1c, R6c, O2o, N3n, N4n, N1p, P2c, Z5c, NO_BOXES
!    
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

       namelist /bfm_IC_nml/ phyto_input,zoop_input,poc_input, &
                          doc_input,phos_input,nit_input,am_input, &
                          oxy_input
       rlength = 1
!
       open(namlst,file='BFM_General.nml',status='old',action='read')!,err=100)
       read(namlst,nml=bfm_IC_nml)!, err=102)
       rewind(namlst)
       close(namlst)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Phytoplankton Carbon IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) P2c(1)
       open(1001,file=phyto_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 1001 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Zooplankton Carbon IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) Z5c(1)
       open(1003,file=zoop_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 1003 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Particulate Organic Carbon IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) R6c(1)
       open(1007,file=poc_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 1007 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Disolved Organic Carbon IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) R1c(1)
       open(1009, file=doc_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 1009 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Phosphate IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) N1p(1)
       open(2000, file=phos_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 2000 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Nitrate IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) N3n(1)
       open(1005, file=nit_input, form='unformatted', access='direct', recl=rlength)
       write(6,*) 'open 1005 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Ammonium IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) N4n(1)
       open(2009,file=am_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 2009 done'
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Oxygen IC
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       inquire(IOLENGTH=rlength) O2o(1)
       open(1000,file=oxy_input,form='unformatted',access='direct',recl=rlength)
       write(6,*) 'open 1000 done'

       write(6,*) 'open units done'
!
!
!    -----READ T&S INITIAL CONDITIONS-----
!
     do k = 1,NO_BOXES

           read (1001,rec=k) P2c(k)
           read (1003,rec=k) Z5c(k)
           read (1007,rec=k) R6c(k)
           read (1009,rec=k) R1c(k)
           read (2000,rec=k) N1p(k)
           read (1005,rec=k) N3n(k)
           !read (2009,rec=k) N4n(k)
           !N1p(k) = N3n(k)/16.0
           !N3n(k) = 0.04
           N4n(k) = 0.0
           read (1000,rec=k) O2o(k)

     enddo

     write(6,*) 'read units done'
!
     return
!
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     !  Define error messages
     ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!100   call error_msg_prn(NML_OPEN,"get_IC.F90","bfm_IC_input.nml")
!102   call error_msg_prn(NML_READ,"get_IC.F90","bfm_IC_input")
!
      end subroutine get_IC
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM17 - 17 Eq. Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
