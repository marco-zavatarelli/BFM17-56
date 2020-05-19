! !ROUTINE: PROFTS
!
! **************************************************************
! **************************************************************
! **                                                          **
! ** ONE-DIMENSIONAL BFM-POM  MODELING SYSTEM (BFM-POM1D)     **
! **                                                          **
! ** The modeling system originate from the direct on-line    **
! ** coupling of the 1D Version of the Princeton Ocean model  **
! ** "POM" and the Biological Flux Model "BFM".               **
! **                                                          **
! ** The whole modelling system and its documentation are     **
! ** available for download from the BFM web site:            **
! **                                                          **
! **                  bfm-community.eu                        **
! **                                                          **
! ** For questions and/or information please address to the   **
! ** BFM system team:                                         **
! **                                                          **
! **                 (bfm_st@lists.cmcc.it)                   **
! **                                                          **
! ** Version 1.0 2016                                         **
! **                                                          **
! ** This release has been finalised by Marco Zavatarelli,    **
! ** Giulia Mussap and Nadia Pinardi. However, previous       **
! ** significant contributions were provided also by          **
! ** Momme Butenschoen and Marcello Vichi.                    **
! ** Thanks are due to Prof. George L. Mellor that allowed us **
! ** to modify, use and distribute the one dimensional        **
! ** version of the Princeton Ocean Model.                    **
! **                                                          **
! **                            Marco.Zavatarelli@unibo.it    **
! **                                                          **
! ** This program is free software; you can redistribute it   **
! ** and/or modify it under the terms of the GNU General      **
! ** Public License as published by the Free Software         **
! ** Foundation.                                              **
! ** This program is distributed in the hope that it will be  **
! ** useful,but WITHOUT ANY WARRANTY; without even the        **
! ** implied warranty of  MERCHANTEABILITY or FITNESS FOR A   **
! ** PARTICULAR PURPOSE.  See the GNU General Public License  **
! ** for more details.                                        **
! ** A copy of the GNU General Public License is available at **
! ** http://www.gnu.org/copyleft/gpl.html or by writing to    **
! ** the Free Software Foundation, Inc. 59 Temple Place,      **
! ** Suite 330, Boston, MA 02111, USA.                        **
! **                                                          **
! **************************************************************
! **************************************************************

!DESCRIPTION
!
! !INTERFACE
!
      Subroutine PROFTS(FF,WFSURF,WFBOT,SWRAD,FSURF,NBC,DT2,NTP,UMOL)
!
! DESCRIPTION
!
! This subroutine solves for the conservative (Temperature and Salinity)
! and non-conservative (BFM state var's) scalars of BFM-POM1D.
! It handles the surface and bottom boundary condition.
! When used to compute temperature it handles also the solar radiation
! penetration along the water column, Based on:
!
! Paulson C. A., Simpson J.J. (1977)
! Irradiance measurements in the upper ocean.
! Journal of Physical Oceanography, 7, 952-956.
!
! Note that when the system is run in diagnostic mode (prescribed
! Temperature and salinity values), the soutine is used only to compute
! the vertical profiles of the non-conservative BFM scalars.
!
! The routine dummy arguments are:
! FF:     Property to be computed
! WFSURF: Property surface flux (for temperature it lacks the incoming
!         surface solar radiation).
! WFBOT:  Property bottom flux.
! SWRAD:  Incoming solar radiation
! FSURF:  Prescribed surface property value
! NBC:    Flag for definition of the surface boundary condition
! DT2:    Twice the Time step.
! NTP:    Flag to choose the Optical (Jerlov) Water type
! UMOL:   Background diffusivity.
!
!****************************************************************************
!
!     -----MODULES (USE OF ONLY IS STRONGLY ENCOURAGED)-----
!
      use global_mem,ONLY: RLEN, ZERO,ONE
!
      use POM,ONLY: H,KB,A,C,KH,DZ,DZZ,VH,VHP,Z,ilong
!
!     -----IMPLICIT TYPING IS NEVER ALLOWED-----
!
      IMPLICIT NONE
!
!     -----TWICE THE TIME STEP
!
      REAL(RLEN)                          :: DT2
!
!     -----SURFACE TEMPERATURE / SALINITY / TRACER-----
!
      REAL(RLEN)                          :: FSURF
!
!     -----SURFACE INCIDENT SHORT WAVE RADIATION-----
!
      REAL(RLEN)                          :: SWRAD
!
!      -----SURFACE/BOTTOM HEAT FLUX LOSS TERM OR SALINITY / TRACER FLUX-----
!
      REAL(RLEN)                          :: WFSURF, WFBOT
!
!     -----FLAG FOR BOUNDARY CONDITION DEFINITION-----
!
!     ************************************************
!     ************************************************
!     **                                            **
!     ** NBC=1: SURF. B.C. IS WFSURF+SWRAD. NO      **
!     **        RADIATIVE PENETRATION.              **
!     ** NBC=2; SURF. B.C. IS WFSURF. SWRAD         **
!     **        PENETRATES WATER COLUMN             **
!     ** NBC=3; SURF. B.C. IS TSURF. NO SWRAD       **
!     **        RADIATIVE PENETRATION               **
!     ** NBC=4; SURF. B.C. IS TSURF. SWRAD          **
!     **        PENETRATES WATER COLUMN             **
!     **                                            **
!     ** NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE   **
!     ** NEGATIVE VALUES WHEN FLUX IS "IN" THE      **
!     ** WATER COLUMN                               **
!     **                                            **
!     ************************************************
!     ************************************************
!
      INTEGER(ilong)                      :: NBC
!
!
!     -----FLAG FOR JERLOV WATER TYPE CHOICE-----
!
!     ************************************************
!     ************************************************
!     **                                            **
!     ** JERLOV WATER TYPE CHOICE IS RELEVANT ONLY  **
!     ** WHEN NBC = 2 OR NBC = 4                    **
!     **                                            **
!     ************************************************
!     ************************************************
!
      INTEGER(ilong)                      :: NTP
!
!     -----BACKGROUND DIFFUSIVITY------
!
      REAL(RLEN)                          :: UMOL
!
!     -----TEMPERATURE/SALINITY/TRACER-----
!
      REAL(RLEN)                          :: FF(KB)
!
!     -----LOOP COUNTERS-----
!
      INTEGER(ilong)                      :: K,KI
!
!     -----SW PROFILE----
!
      REAL(RLEN)                          :: RAD(KB)
!
!     -----IRRADIANCE PARAMETERS AFTER PAULSON & SIMPSON JPO 1977, 952-956-----
!
      ! REAL(RLEN)                          :: RCP(5), AD1(5), AD2(5)
      REAL(RLEN)                          :: RP(5), AD1(5), AD2(5)
!
!     -----JERLOV WATER TYPES-----
!
!     NTP        = 1           2            3           4          5
!    JERLOV TYPE = I           IA           IB          II         III
!
     DATA RP  /    0.58_RLEN,  0.62_RLEN,   0.67_RLEN,  0.77_RLEN, 0.78_RLEN/
     DATA AD1 /    0.35_RLEN,  0.60_RLEN,   1.00_RLEN,  1.50_RLEN, 1.40_RLEN/
     DATA AD2 /   23.00_RLEN, 20.00_RLEN,  17.00_RLEN, 14.00_RLEN, 7.90_RLEN/
!
!     -----INTRINSIC FUNCTION-----
!
      INTRINSIC EXP
!
!     -----START COMPUTATION OF VERTICAL PROFILE----
!
      DO K = 2,KB - 1
!
         A(K-1) = -DT2 * (KH(K)+UMOL)/ (DZ(K-1) * DZZ(K-1) * H * H)
         C(K)   = -DT2 * (KH(K)+UMOL)/ (DZ(K)   * DZZ(K-1) * H * H)
!
      END DO
!
      RAD(:)=ZERO
!
!     -----SURFACE BOUNDARY CONDITION-----
!
      select case (NBC)
!
         case (2, 4)
!
!        ***********************************************************
!        ***********************************************************
!        **                                                       **
!        ** PENETRATIVE RADIATION CALCULATION.                    **
!        ** AT THE BOTTOM ANY                                     **
!        ** UNATTENUATED IS DEPOSITED IN THE BOTTOM LAYER         **
!        **                                                       **
!        ***********************************************************
!        ***********************************************************
!
         RAD(:) = SWRAD * ((    RP(NTP) * EXP(Z(:)*H/AD1(NTP))               &
                        +  (ONE-RP(NTP) * EXP(Z(:)*H/AD2(NTP)))))
!
         RAD(KB)=ZERO
!
      end select

      select case (NBC)

         case (1)
!
              VH(1)  =  A(1)/(A(1)-ONE)
              VHP(1) = -DT2*(WFSURF+SWRAD)/(-DZ(1)*H) - FF(1)
              VHP(1) =  VHP(1)/(A(1)-ONE)
!
              RAD(:) = ZERO
!
         case (2)
!
              VH(1)  = A(1)/(A(1)-ONE)
              VHP(1) = DT2*(WFSURF+RAD(1)-RAD(2))/(DZ(1)*H) -FF(1)
              VHP(1) = VHP(1)/(A(1)-ONE)
!
         case (3, 4)

              VH(1)  = ZERO
              VHP(1) = FSURF
!
      end select
!
!    ********************************************************
!    ********************************************************
!    **                                                    **
!    ** The following section solves the equation          **
!    ** DT2*(KH*FF')' -FF = -FB                            **
!    **                                                    **
!    ********************************************************
!    ********************************************************
!
      DO K = 2,KB - 2
!
        VHP(K) = ONE/(A(K)+C(K)*(ONE-VH(K-1))-ONE)
        VH(K)  = A(K)*VHP(K)
        VHP(K) = (C(K)*VHP(K-1)-FF(K)+DT2*(RAD(K)-RAD(K+1))/(H*DZ(K)))*VHP(K)
        
!
      END DO
!
!     -----APPLY A NON ADIABATIC BOTTOM BOUNDARY CONDITION-----
!
      FF(KB-1) = (C(KB-1)*VHP(KB-2)-FF(KB-1)+(WFBOT*DT2/(DZ(KB-1)*H)))/ &
                 (C(KB-1)* (ONE-VH(KB-2))-ONE)
!
!     -----FINAL SCALAR COMPUTATION-----
!
      DO K = 2,KB - 1
!
           KI = KB - K
           FF(KI) = VH(KI)*FF(KI+1) + VHP(KI)
!
      END DO
!
!     -----ZEROING-----
!
      VH(:)  = ZERO
      VHP(:) = ZERO
      A(:)   = ZERO
      C(:)   = ZERO
!
      RETURN
!
     end subroutine PROFTS

!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  POM - Princeton Ocean Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
