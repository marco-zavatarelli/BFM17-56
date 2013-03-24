!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: FilterFeeder
!
! DESCRIPTION
!   Parameter value file for suspension feeders (Y3)
!   
!INTERFACE#
  module mem_FilterFeeder
!
! !USES:
#ifdef INCLUDE_BEN

  use global_mem

!  
!
! !AUTHORS
!   ERSEM group, HBB
!
!
!
! !REVISION_HISTORY
!   !
!
! COPYING
!   
!   Copyright (C) 2013 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! FilterFeeder PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer     :: sw_uptake ! 1=ERSEM 2=Modified ..Uptake 3= Mode-FIOed Uptake with filterlimitmit.
  real(RLEN)  :: p_dwat  ! Food layer in the water
  real(RLEN)  :: p_chu  ! Upper limit of total food
  real(RLEN)  :: p_clu  ! Lower limit of availability of a food source
  real(RLEN)  :: p_su  ! Growth rate
  real(RLEN)  :: p_q10  ! q10
  real(RLEN)  :: p_Rps  ! pseudefaxes production
  real(RLEN)  :: p_R6  ! Food matrix detritus on pelagic det. (R6)
  real(RLEN)  :: p_puePI  ! Excreted fraction of phytoplankton uptake
  real(RLEN)  :: p_pueZI  ! Excreted fraction of microzooplankton uptake
  real(RLEN)  :: p_pueQ6  ! Excreted fraction of detritus uptake
  real(RLEN)  :: p_srr  ! Relative respiration rate
  real(RLEN)  :: p_sra  ! respired Part of uptake  used for filtering
  real(RLEN)  :: p_pur  ! respired Part of uptake  used for digesting
  real(RLEN)  :: p_sd  ! Specific Mortality
  real(RLEN)  :: p_sd2  ! Specific Density Dependent Mortality (mort. o1 0.1 at 25000)
  real(RLEN)  :: p_qn  ! Fixed nutrient quotum N:C
  real(RLEN)  :: p_qp  ! Fixed nutrient quotum P:C
  real(RLEN)  :: p_clm  ! Upper depth of accessed sediment layer
  real(RLEN)  :: p_cm  ! Lower  depth of accessed sediment layer
  real(RLEN)  :: p_puQ6  ! Food matrix detritus on the sediment (Q6)
  real(RLEN)  :: p_PI  ! Food parameter Y3 on phytoplankton
  real(RLEN)  :: p_ZI  ! Food parameter Y3 on microzooplankton
  real(RLEN)  :: p_vum  ! Volume filtered by 1mgC Y3 
  real(RLEN)  :: p_clO2o  ! oxygen con. at which activity is lowered to half 
  real(RLEN)  :: p_height ! height of the layer from which is filtered
  real(RLEN)  :: p_max ! proportion of sedimentation entering gridlayer which can be used for food uptake.
  real(RLEN)  :: p_pePel =0.0 ! part of excretion and respiration which is coupled to pelagic 
  real(RLEN)  :: p_pR6Pel =0.0 ! part of produced R6 which is excreted to the pelagic 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitFilterFeeder
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitFilterFeeder()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /FilterFeeder_parameters/ sw_uptake, p_dwat, p_su, p_q10, p_Rps, p_R6, p_clu, p_chu, &
    p_puePI, p_pueZI, p_pueQ6, p_srr, p_sra, p_pur, p_sd, p_sd2, p_qn, p_qp, p_clm, p_cm, p_puQ6, &
    p_PI,p_ZI,p_vum,p_clO2o,p_height,p_max,p_pePel,p_pR6Pel
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
    write(LOGUNIT,*) "#  Reading FilterFeeder parameters.."
    open(NMLUNIT,file='FilterFeeder.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=FilterFeeder_parameters,err=101)
    close(NMLUNIT)
    select case (sw_uptake)
     case(1)
         write(LOGUNIT,*) "Uptake according original ERSEM"
         if  ( p_chu ==0.0 .or. p_dwat ==0.0 .or. p_Rps > 0.0 ) then
         write(LOGUNIT,*) "w_uptake==1 and p_chu==0.0 or p_dwat ==0 : wrong combination"
         write(LOGUNIT,*) "w_uptake==1 and p_Rps>0.0"
         goto 101
        endif
     case(2) 
        write(LOGUNIT,*) "Uptake according Holling modified response"
        if  ( p_vum ==0.0 ) then
           write(LOGUNIT,*) "w_uptake=2 and p_vum==0 : wrong combination"
           goto 101
        endif
     case(3) 
        write(LOGUNIT,*) "Uptake according Holling modified response"
        write(LOGUNIT,*) "Lower Threhold in food uptake is set by comparing energy gain with loss"
        if  ( p_vum ==0.0.or.p_sra==0.0.or.p_clu >0.0) then
          write(LOGUNIT,*) "w_uptake=3 and p_vum==0 or p_sra==0.0: wrong combination"
          if ( p_clu > 0.0 ) then
             write(LOGUNIT,*)"p_sra is indirectly used as a threshold instead p_chu"
             write(LOGUNIT,*)"Error: p_clu > 0.0 !"
          endif
          goto 101
        endif
     end select
         
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=FilterFeeder_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitFilterFeeder.f90","FilterFeeder.nml")
101 call error_msg_prn(NML_READ,"InitFilterFeeder.f90","FilterFeeder_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitFilterFeeder
#endif
  end module mem_FilterFeeder
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
