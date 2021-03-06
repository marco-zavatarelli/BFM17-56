#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelChem
!
! DESCRIPTION
!       This process describes the additional dynamics of dissolved
!       compounds in the watercolumn. Parameterized processes are:
!       - nitrification
!       - denitrification
!       - reoxidation of reduction equivalents
!       - dissolution of biogenic silica
!       This function also calls the carbonate system dynamics
!       (INCLUDE_PELCO2) and iron dynamics (INCLUDE_PELFE)
!       if activated
!
! !INTERFACE
  subroutine PelChemDynamics
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_mem, ONLY:RLEN,ZERO
#ifdef NOPOINTERS
  use mem
#else
  use mem,  ONLY: N4n, N3n, O2o, O4n, N6r, R6s, N5s, P1s
  use mem, ONLY: ppN4n, ppN3n, ppO2o, ppO4n, ppN6r, ppR6s, ppN5s,    &
    flN3O4n, ETW, flPTN6r, NO_BOXES, iiBen, iiPel, flN4N3n, &
    flux_vector, ppO3c, ppR6c, ppR1c, ppR2c,ppR3c, ppR6p, ppR1p, ppR6n, &
    ppR1n, R6c, R6p, R1c, R1p, R1n, R2c, ppN1p
#ifdef INCLUDE_PELCO2
  use mem, ONLY: ppO3c,CalcPelBacteria
#endif
#endif
  use mem_Param,  ONLY: p_qon_nitri, p_qro, p_qon_dentri, p_small
  use mem_PelChem
  use mem_globalfun,   ONLY: MM_vector, eTq_vector, insw_vector
#ifdef INCLUDE_PELCO2
  use mem_CO2, ONLY: CalcBioAlkFlag
#endif
  use constants, ONLY: MW_C
  use POM, ONLY: H
!
!
! !AUTHORS
!   Original version by P. Ruardij and M. Vichi
!
!
!
! !REVISION_HISTORY
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, save :: first =0
  integer       :: AllocStatus, DeallocStatus
  real(RLEN),allocatable,save,dimension(:) :: fN4N3n,fN6O2r,eo,     &
                                              er,osat,rPAo,fR6N5s
#ifndef INCLUDE_PELCO2
  integer,parameter :: ppO3c = 0
#endif

  real(RLEN),allocatable,save,dimension(:) :: fRIO3c, fRIN1p, fRIN4n, bacprof
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  if (first==0) then
     first=1
     allocate(fN6O2r(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fN6O2r"
     allocate(eo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating eo"
     allocate(er(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating er"
     allocate(rPAo(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating rPAo"
     allocate(fR6N5s(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fR6N5s"
     allocate(osat(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating osat"
     allocate(fN4N3n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fN4N3n"
   if (.NOT.CalcPelBacteria(1)) then
     allocate(fRIO3c(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fRIO3c"
     allocate(fRIN1p(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fRIN1p"
     allocate(fRIN4n(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating fRIN4n"
     allocate(bacprof(NO_BOXES),stat=AllocStatus)
     if (AllocStatus  /= 0) stop "error allocating bacprof"
   end if
  end if

#ifdef INCLUDE_PELCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Carbonate chemistry
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  call PelCO2Dynamics( )
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Regulating factors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  eo  =   MM_vector(  max(p_small,O2o(:)),  p_clO2o)
  er  =   MM_vector(  N6r(:),  p_clN6r)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Nitrification in the water
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  flN4N3n(:) =  max(ZERO,p_sN4N3* N4n(:)* eTq_vector(  ETW(:),  p_q10N4N3)* eo)
  call flux_vector( iiPel, ppN4n,ppN3n, flN4N3n(:) )
  call flux_vector( iiPel, ppO2o,ppO2o,-( flN4N3n(:)* p_qon_nitri) )

#ifndef BFM17
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Denitrification in the water
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  rPAo  =   flPTN6r(:)/ p_qro
  flN3O4n(:) = max(ZERO,p_sN3O4n* eTq_vector( ETW(:), p_q10N4N3)* er* rPAo/ p_rPAo* &
               N3n(:))
  call flux_vector( iiPel, ppN3n,ppO4n, flN3O4n(:) )
  call flux_vector( iiPel, ppN6r,ppN6r,-( p_qro* flN3O4n(:)* p_qon_dentri* &
  insw_vector( -( O2o(:)- N6r(:)/ p_qro))) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Reoxidation of reduction equivalents
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  fN6O2r  =   p_rOS* N6r(:)* eo
  call flux_vector( iiPel, ppN6r,ppN6r,-( fN6O2r) )
  call flux_vector( iiPel, ppO2o,ppO2o,-( fN6O2r/ p_qro) )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Dissolution of biogenic silicate
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  fR6N5s  =   p_sR6N5* eTq_vector(  ETW(:),  p_q10R6N5)* R6s(:)
  call flux_vector( iiPel, ppR6s,ppN5s, fR6N5s )
#endif

#ifdef INCLUDE_PELCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Corrections of nitrogen cycle biogeochemistry on Total Alkalinity
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if ( ppO3c > 0 .and. CalcBioAlkFlag)  call AlkalinityDynamics( )
#endif

#ifdef INCLUDE_PELFE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Iron Chemistry (dissolution and scavenging)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call PelIronDynamics()
#endif

if (.NOT.CalcPelBacteria(1)) then
!!CASE for constant oranic matter remineralisation
     bacprof(:) = (1/log(H))*(log(H) - log(H - Depth(:)))
     !fRIO3c=p_sR6O3*bacprof(:)*R6c(:)
     fRIO3c=p_sR6O3*R6c(:)
     call flux_vector(iiPel,ppR6c,ppO3c,fRIO3c)
     call flux_vector(iipel,ppO2o,ppO2o,-fRIO3c/MW_C)
     !fRIO3c=p_sR1O3*bacprof(:)*R1c(:)
     fRIO3c=p_sR1O3*R1c(:)
     call flux_vector(iiPel,ppR1c,ppO3c,fRIO3c)
     call flux_vector(iipel,ppO2o,ppO2o,-fRIO3c/MW_C)
#ifndef BFM17
     fRIO3c=P_sR2O3*R2c(:)
     call flux_vector(iiPel,ppR2c,ppO3c,fRIO3c)
     call flux_vector(iipel,ppO2o,ppO2o,-fRIO3c/MW_C)
#endif
     !fRIN1p=p_sR6N1*bacprof(:)*R6p(:)
     fRIN1p=p_sR6N1*R6p(:)
     call flux_vector(iiPel,ppR6p,ppN1p,fRIN1p)
     fRIN1p=P_sR1N1*bacprof(:)*R1p(:)
     fRIN1p=P_sR1N1*R1p(:)
     call flux_vector(iiPel,ppR1p,ppN1p,fRIN1p)
     !fRIN4n=p_sR6N4*bacprof(:)*R6n(:)
     fRIN4n=p_sR6N4*R6n(:)
     call flux_vector(iiPel,ppR6n,ppN4n,fRIN4n)
     !fRIN4n=p_sR1N4*bacprof(:)*R1n(:)
     fRIN4n=p_sR1N4*R1n(:)
     call flux_vector(iiPel,ppR1n,ppN4n,fRIN4n)
#ifndef BFM17
     fRIO3c=p_sR3O3*R3c(:)
     call flux_vector(iiPel,ppR3c,ppO3c,fRIO3c)
     call flux_vector(iipel,ppO2o,ppO2o,-fRIO3c/MW_C)
#endif
end if

  end subroutine PelChemDynamics
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
