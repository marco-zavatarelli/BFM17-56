!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: set_initial_conditions
!
!DESCRIPTION
!
! This routine assigns initial conditioons of biochemical variables in POM
!
!********************************************************************************
!********************************************************************************
!**                                                                            **
!** IN THIS DEFAULT VERSION THE BFM STATE VAR'S ARE ASSIGNED THE VALUED SET IN **
!** THE NAMELIST "bfm_init_nml",  BUT OBVIOUSLY THIS ROUTINE COULD/SHOULD BE   **
!** TAILORED ACCCORDING TO THE IMPLEMENTATION CHARACTERISTIC AND TO THE DATA   **
!** AVAILABILITY. IF THE USER HAS A DATA FILE TO BE READ, THE READING          **
!** STATEMENT COULD BE PUT HERE IN SUBSTITUTION OF THE HARDWIRED I.C.'s        **
!**                                                                            **
!**                                              Marco.Zavatarelli@unibo.it    **
!**                                                                            **
!********************************************************************************
!********************************************************************************
!
! !INTERFACE
      subroutine set_initial_conditions_bfm17

!   Momme Butenschon, May 2004
!Marco Zavatarelli, Giulia Mussap, 2014
!   Dipartimento di Fisica
!   Universita di Bologna
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Modules (use of ONLY is strongly encouraged!)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      use global_mem, ONLY: ZERO,NML_OPEN,NML_READ,NMLUNIT,error_msg_prn
      use Mem
      use POM
      use Service
!-------------------------------------------------------------------------!
!BOC
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Implicit typing is never allowed
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      IMPLICIT NONE
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Local Variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      integer :: ib
!
      real(RLEN) :: dd, d1,d2,d1cc,d1cn,d1cp, d1cs, d1ci
      real(RLEN),dimension(KB) :: d2cc
!
      real(RLEN),parameter :: p_nRc=0.0126,p_pRc=0.7862e-3,p_sRc=0.0118,&
                              p_iRc=1./25.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Definition of Initial Pelagic (D3) state variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
     real(RLEN) :: O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, B1n0, B1p0, &
      P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, P3n0, P3p0, &
      P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, Z4p0, Z5c0, &
      Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, R6c0, R6n0, &
      R6p0, R6s0, O3c0, O3h0
!
     namelist /bfm_init_nml/ O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, &
      B1n0, B1p0, P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, &
      P3n0, P3p0, P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, &
      Z4p0, Z5c0, Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, &
      R6c0, R6n0, R6p0, R6s0, O3c0, O3h0
!
#ifdef INCLUDE_BEN
!
     real(RLEN) :: G3c0, G3h0, G13c0, G13h0, G23c0, G23h0, Y1c0, Y1n0, Y1p0, &
      Y2c0, Y2n0, Y2p0, Y3c0, Y3n0, Y3p0, Y4c0, Y4n0, Y4p0, Y5c0, Y5n0, Y5p0, &
      Q6c0, Q6n0, Q6p0, Q6s0, Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, H1c0, H1n0, &
      H1p0, H2c0, H2n0, H2p0, K1p0, K11p0, K21p0, K4n0, K14n0, K24n0, K6r0, &
      K16r0, K26r0, K3n0, K5s0, G2o0, G4n0, D1m0, D2m0, D6m0, D7m0, D8m0, D9m0
!
     namelist /bfm_init_nml_ben/ G3c0, G3h0, G13c0, G13h0, G23c0, G23h0, Y1c0, &
      Y1n0, Y1p0, Y2c0, Y2n0, Y2p0, Y3c0, Y3n0, Y3p0, Y4c0, Y4n0, Y4p0, Y5c0, &
      Y5n0, Y5p0, Q6c0, Q6n0, Q6p0, Q6s0, Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, &
      H1c0, H1n0, H1p0, H2c0, H2n0, H2p0, K1p0, K11p0, K21p0, K4n0, K14n0, K24n0, &
      K6r0, K16r0, K26r0, K3n0, K5s0, G2o0, G4n0, D1m0, D2m0, D6m0, D7m0, D8m0, &
      D9m0
!
#endif
!
     open(NMLUNIT,file='BFM_General.nml',status='old',action='read',err=100)
     read(NMLUNIT,nml=bfm_init_nml,err=101)
!
!
#ifdef INCLUDE_BEN
!
     read(NMLUNIT,nml=bfm_init_nml_ben,err=102)
!
#endif
!
     close(NMLUNIT)
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Definition of biogeochemical  global variables
! IrrOPT  in the equation of Steele and light
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
        eir(:) = ZERO
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Definition of general pelagic state variables:
! pelagic gases
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!        if(set_IC_from_file)then
!           write(6,*) 'set_IC_from_file .TRUE.'
           call get_IC
!        else
!            write(6,*) 'set_IC_from_file .FALSE.'
!           O2o(:) = O2o0
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! pelagic nutrients (mMol /m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!           n1p(:) = N1p0
!           n3n(:) = N3n0
!           n4n(:) = N4n0
!           R6c(:) = R6c0
!           R1c(:) = R1c0
!           Z4c(:) = Z4c0
!        endif
        n5s(:) = N5s0
        O4n(:) = O4n0
        N6r(:) = N6r0
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! pelagic detritus  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
        R6n(:) = r6c(:)*p_nRc
        R6p(:) = r6c(:)*p_pRc
        R6s(:) = r6c(:)*p_sRc
        R2c(:) = R2c0
        R3c(:) = R3c0
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! dissolved organic matter
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
         R1n(:) = R1c(:)*p_nRc*0.5
         R1p(:) = R1c(:)*p_pRc*0.5
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! State variables for phytoplankton model
! pelagic diatoms  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
        d1cc=P1c0
        d1cn=d1cc*p_nRc
        d1cp=d1cc*p_pRc
        d1cs=d1cc*p_sRc
        d1ci=d1cc*p_iRc
      if (CalcPhytoPlankton(1)) then
        P1c(:) = d1cc
        P1n(:) = d1cn
        P1p(:) = d1cp
        P1s(:) = d1cs
        P1l(:) = d1ci
      else
        P1c(:) = ZERO
        P1n(:) = ZERO
        P1p(:) = ZERO
        P1s(:) = ZERO
        P1l(:) = ZERO
      end if
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! pelagic flagellates  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      d1cc = P2c0
      d1cn = d1cc*p_nRc
      d1cp = d1cc*p_pRc
      d1ci = d1cc*p_iRc*.5
      if (CalcPhytoPlankton(2)) then
!         if(set_IC_from_file)then
            d2cc(:) = P2c(:)
            P2n(:)  = d2cc(:)*p_nRc
            P2p(:)  = d2cc(:)*p_pRc
            P2l(:)  = d2cc(:)*p_iRc*.5
!         else
!            P2c(:) = d1cc
!            P2n(:) = d1cn
!            P2p(:) = d1cp
!            P2l(:) = d1ci
!         endif
      else
        P2c(:) = ZERO
        P2n(:) = ZERO
        P2p(:) = ZERO
        P2l(:) = ZERO
      end if
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! picophytoplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
        d1cc = P3c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
      if (CalcPhytoPlankton(3)) then
        P3c(:) = d1cc
        P3n(:) = d1cn
        P3p(:) = d1cp
        P3l(:) = d1ci
      else
        P3c(:) = ZERO
        P3n(:) = ZERO
        P3p(:) = ZERO
        P3l(:) = ZERO
      end if
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Large phytoplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
        d1cc = P4c0
        d1cn = d1cc*p_nRc
        d1cp = d1cc*p_pRc
        d1ci = d1cc*p_iRc*.5
      if (CalcPhytoPlankton(4)) then
        P4c(:) = d1cc
        P4n(:) = d1cn
        P4p(:) = d1cp
        P4l(:) = d1ci
      else
        P4c(:) = ZERO
        P4n(:) = ZERO
        P4p(:) = ZERO
        P4l(:) = ZERO
      end if
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! State variables for mesozooplankton model
! carnivorous mesozooplankton ( mg C/m3 )
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      if (CalcMesoZooPlankton(1)) then
        Z3c(:) = Z3c0
        Z3n(:) = Z3c(:)*p_nRc
        Z3p(:) = Z3c(:)*p_pRc
      else
        Z3c(:) = ZERO
        Z3n(:) = ZERO
        Z3p(:) = ZERO
      end if
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! omnivorous mesozooplankton ( mg C/m3 )
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      if (CalcMesoZooPlankton(2)) then
        Z4c(:) = Z4c0
        Z4n(:) = Z4c(:)*p_nRc
        Z4p(:) = Z4c(:)*p_pRc
      else
        Z4c(:) = ZERO
        Z4n(:) = ZERO
        Z4p(:) = ZERO
       end if

! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! State variables for microzooplankton model
! pelagic microzooplankton  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      if (CalcMicroZooPlankton(1)) then
!         Z5c(:) = Z5c0
         Z5n(:) = Z5c(:)*p_nRc
         Z5p(:) = Z5c(:)*p_pRc
      else
         Z5c(:) = ZERO
         Z5n(:) = ZERO
         Z5p(:) = ZERO
      end if
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! heterotrophic flagellates (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      if (CalcMicroZooPlankton(2)) then
         Z6c(:) = Z6c0
         Z6n(:) = Z6c(:)*p_nRc
         Z6p(:) = Z6c(:)*p_pRc
      else
         Z6c(:) = ZERO
         Z6n(:) = ZERO
         Z6p(:) = ZERO
      end if
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! State variables for pelagic bacteria model B1
! pelagic bacteria  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
      if (CalcPelBacteria(1)) then
        B1c(:) = B1c0
        B1n(:) = B1c(:)*p_nRc
        B1p(:) = B1c(:)*p_pRc
      else
        B1c(:) = ZERO
        B1n(:) = ZERO
        B1p(:) = ZERO
      end if
!
!
#ifdef INCLUDE_BEN
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! State variables for the benthic modules
! zoobenthos
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      if (CalcBenOrganisms(1)) then
        Y1c(:) = Y1c0
      else
        Y1c(:) = ZERO
      end if

      if (CalcBenOrganisms(2)) then
        Y2c(:) = Y2c0
      else
        Y2c(:) = ZERO
      end if

      if (CalcBenOrganisms(3)) then
        Y3c(:) = Y3c0
      else
        Y3c(:) = ZERO
      end if

      if (CalcBenOrganisms(4)) then
        Y4c(:) = Y4c0
      else
        Y4c(:) = ZERO
      end if

      if (CalcBenOrganisms(5)) then
        Y5c(:) = Y5c0
      else
        Y5c(:) = ZERO
      end if

! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! bacteria
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      if (CalcBenBacteria(1)) then
        H1c(:) = H1c0
      else
        H1c(:) = ZERO
      end if

      if (CalcBenBacteria(2)) then
        H2c(:) = H2c0
      else
        H2c(:) = ZERO
      end if

! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! benthic nutrients
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
       K5s(:)  = 20.75
       K6r(:)  = K6r0
       k4n(:)  = K4n0
       K14n(:) = K14n0
       k24n(:) = K24n0
       k1p(:)  = K1p0
       K11p(:) = K11p0
       K21p(:) = K21p0
       K3n(:)  = K3n0
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! benthic detritus  (respectively mg C/m3 mMol N/m3 mMOL P/m3)
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
       Q1c(:)  = Q1c0
       Q11c(:) = Q11c0
!
#ifdef IMPFLUX
!
       Q6c(:) = 1.E9
       Q6n(:) = 1.E9
       Q6p(:) = 1.E9
       Q6s(:) = 1.E9
!
#else
!
       Q6c(:) = Q6c0
       Q6n(:) = Q6n0
       Q6p(:) = Q6p0
       Q6s(:) = Q6s0
!
#endif
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! gases
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
       G2o(:) = G2o0
       G3c(:) = G3c0
       G4n(:) = 37.0
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! layers
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
        D1m(:) = D1m0
        D2m(:) = D2m0
        D6m(:) = D6m0
        D7m(:) = D7m0
        D8m(:) = D8m0
        D9m(:) = D9m0
!
#endif
!
      return
 100 call error_msg_prn(NML_OPEN,"InitParam.f90","BFM_General.nml")
 101 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters")
 102 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters_ben")

end subroutine set_initial_conditions_bfm17
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
