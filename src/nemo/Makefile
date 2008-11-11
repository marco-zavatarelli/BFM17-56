#
# Makefile to build the pelagic configuration of the BFM model coupled with NEMO
# Author: M. Vichi (vichi@bo.ingv.it)
# 
# $NEMOLIBDIR is obtained from the parent Makefile via export
#
BFMLIB     = $(NEMOLIBDIR)/oce/libopa.a

# BFMDIR path is an environmental variable
ifndef BFMDIR
  $(error The environment variable BFMDIR is not defined!)
endif
BFMSRC 	   = $(BFMDIR)/src/BFM
BFMNEMOSRC = $(BFMDIR)/src/nemo
BFMSHARE   = $(BFMDIR)/src/share

## -------------------------------------------------------------------------
## BFM CONFIGURATION 
## -------------------------------------------------------------------------
DEFINES += -DDUMMY
# Pelagic CO2 flags (activates compilation and macros, true by default)
INCLUDE_PELCO2=true
ifeq ($(INCLUDE_PELCO2),true)
  DEFINES += -DINCLUDE_PELCO2
endif

## NEMO modules (with absolute path)
F_O += -I$(NEMOLIBDIR)/oce
## BFM include files and pp macros
BFMINCDIR = $(BFMDIR)/src/BFM/include
F_O += -I$(BFMINCDIR) -I$(BFMDIR)/include -I$(NCDF_INC) -I$(NEMOWORKDIR) -I$(NEMOLIBDIR)/oce
P_P += -DBFM_NEMO  $(DEFINES)
## specific options for NEC-SX
ifeq ($(FORTRAN_COMPILER),NECSX6)
AR=${SXAR}
P_P += -DBFM_NOPOINTERS -DNECSX -DBFM_PARALLEL 
INLINE = -pi auto line=200 expin=$(BFMDIR)/src/BFM/General,$(BFMDIR)/src/share,$(BFMDIR)/src/BFM/PelB,$(BFMDIR)/src/BFM/CO2 exp=sourcesink_flux_vector,flux_vector,flux,fixed_quota_flux_vector,Source_D3_vector nest=10 rexp=PhytoDynamics,MicroZooDynamics,MesoZooDynamics,fixed_quota_flux_vector,flux_vector,Source_D3_vector,calcco2system
P_P += -R2 -EP -sx6 -Cvopt $(INLINE) -Wf,-P nh -Wf,-A idbl4 -Wf,-pvctl expand=5 noassume loopcnt=500000 fullmsg -Wf,-L fmtlist transform map summary noinclist 
endif

BFM_PELMOD = \
	${BFMLIB}(${BFMSRC}/General/ModuleGlobalMem.o)			\
	${BFMLIB}(${BFMSRC}/General/ModuleConstants.o)			\
	${BFMLIB}(${BFMSRC}/General/ModuleGlobFun.o)			\
	${BFMLIB}(${BFMSRC}/General/bfm_error_msg.o)			\
	${BFMLIB}(${BFMSRC}/General/ModuleMem.o)			\
	${BFMLIB}(${BFMSRC}/General/ModuleParam.o)			\
	${BFMLIB}(${BFMSRC}/General/ModuleInterface.o)			\
	${BFMLIB}(${BFMSRC}/Light/ModuleLightAdaptation.o)		\
	${BFMLIB}(${BFMSRC}/Light/ModulePhotoAvailableRadiation.o)	\
	${BFMLIB}(${BFMSRC}/Oxygen/ModuleWindOxReaeration_3.o)		\
	${BFMLIB}(${BFMSRC}/PelB/ModuleMesoZoo.o)			\
	${BFMLIB}(${BFMSRC}/PelB/ModuleMicroZoo.o)			\
	${BFMLIB}(${BFMSRC}/PelB/ModulePelBac.o)			\
	${BFMLIB}(${BFMSRC}/PelB/ModulePelChem.o)			\
	${BFMLIB}(${BFMSRC}/PelB/ModulePelGlobal.o)			\
	${BFMLIB}(${BFMSRC}/PelB/ModulePhyto.o)				\
	${BFMLIB}(${BFMSRC}/PelBen/ModuleSettling.o)


BFMNEMO_OBJ   = \
	${BFMLIB}(${BFMSHARE}/api_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/netcdf_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/ResetFluxes.o)		\
	${BFMLIB}(${BFMSHARE}/calcmean_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/string_functions.o)	\
	${BFMLIB}(${BFMSHARE}/init_cnps.o)		\
	${BFMLIB}(${BFMSHARE}/init_var_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/make_flux_output.o)	\
	${BFMLIB}(${BFMNEMOSRC}/trc_ini_bfm.o)		\
	${BFMLIB}(${BFMNEMOSRC}/trcbfm.o)		\
	${BFMLIB}(${BFMNEMOSRC}/envforcing_bfm.o)	\
	${BFMLIB}(${BFMNEMOSRC}/trcrad_bfm.o)		\
	${BFMLIB}(${BFMNEMOSRC}/trc_trp_bfm.o)		\
	${BFMLIB}(${BFMNEMOSRC}/trc_set_bfm.o)		\
	${BFMLIB}(${BFMNEMOSRC}/trc_sink_bfm.o)		\
	${BFMLIB}(${BFMNEMOSRC}/trc_sink_muscl_bfm.o)	\
	${BFMLIB}(${BFMNEMOSRC}/trc_dia_bfm.o)		\
	${BFMLIB}(${BFMNEMOSRC}/GetDelta.o)		\
	${BFMLIB}(${BFMNEMOSRC}/D3toD1.o)		\
	${BFMLIB}(${BFMNEMOSRC}/D2toD1.o)

BFM_PELOBJ = \
	${BFMLIB}(${BFMSRC}/General/set_var_info_bfm.o)			\
	${BFMLIB}(${BFMSRC}/General/AllocateMem.o)			\
	${BFMLIB}(${BFMSRC}/General/Ecology.o)				\
	${BFMLIB}(${BFMSRC}/General/InitBoxParams.o)			\
	${BFMLIB}(${BFMSRC}/General/Initialize.o)			\
	${BFMLIB}(${BFMSRC}/General/InitTransportStateTypes.o)		\
	${BFMLIB}(${BFMSRC}/General/eTq.o)				\
	${BFMLIB}(${BFMSRC}/General/CheckMassConservation.o)		\
	${BFMLIB}(${BFMSRC}/Light/LightAdaptation.o)			\
	${BFMLIB}(${BFMSRC}/Light/PhotoAvailableRadiation.o)		\
	${BFMLIB}(${BFMSRC}/Oxygen/WindOxReaeration_3.o)		\
	${BFMLIB}(${BFMSRC}/Oxygen/CalcOxygenSaturation_3.o)		\
	${BFMLIB}(${BFMSRC}/Oxygen/CalcSchmidtNumberOx.o)		\
	${BFMLIB}(${BFMSRC}/PelB/CalcChlorophylla.o)			\
	${BFMLIB}(${BFMSRC}/PelB/CalcVerticalExtinction.o)		\
	${BFMLIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${BFMLIB}(${BFMSRC}/PelB/MesoZoo.o)				\
	${BFMLIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${BFMLIB}(${BFMSRC}/PelB/PelBac.o)				\
	${BFMLIB}(${BFMSRC}/PelB/PelChem.o)				\
	${BFMLIB}(${BFMSRC}/PelB/PelGlobal.o)				\
	${BFMLIB}(${BFMSRC}/PelB/PelagicSystem.o)			\
	${BFMLIB}(${BFMSRC}/PelB/Phyto.o)				\
	${BFMLIB}(${BFMSRC}/PelBen/BentoPelCoup.o)			\
	${BFMLIB}(${BFMSRC}/PelBen/Settling.o)

ifeq ($(INCLUDE_PELCO2),true)
BFM_PELCO2MOD = \
	${BFMLIB}(${BFMSRC}/CO2/ModuleCO2.o)		\
	${BFMLIB}(${BFMSRC}/CO2/ModuleCO2System.o)
BFM_PELCO2OBJ = \
	${BFMLIB}(${BFMSRC}/CO2/PelCO2.o)		\
	${BFMLIB}(${BFMSRC}/CO2/Alkalinity.o)		\
	${BFMLIB}(${BFMSRC}/CO2/CO2Flux.o)		\
	${BFMLIB}(${BFMSRC}/CO2/CO2System_vector.o)	
endif

.PHONY: all depend clean

# Include the dependency-list created by makedepf90 below
# This part is still experimental and does not work yet
#include .depend

all: $(BFM_PELMOD) $(BFM_PELCO2MOD) $(BFMNEMO_OBJ) $(BFM_PELOBJ) \
	$(BFM_PELCO2OBJ) depend 

$(BFM_PELMOD) : $(BFMSRC)/General/ModuleMem.F90

${BFMSRC}/General/ModuleMem.F90 : ${BFMSRC}/General/GlobalDefsBFM.model \
				${BFMSRC}/General/FluxFunctions.F90
	$(BFMSRC)/scripts/GenerateGlobalBFMF90Code  $(DEFINES) \
		-read $(BFMSRC)/General/GlobalDefsBFM.model \
		-from $(BFMSRC)/proto -to $(BFMSRC)/General \
		-actions statemem allocmem netcdfmem \
                -to $(BFMSRC)/include -actions headermem
${BFMSRC}/General/AllocateMem.F90: $(BFMSRC)/General/ModuleMem.F90
${BFMSRC}/General/set_var_info_bfm.F90: $(BFMSRC)/General/ModuleMem.F90

depend .depend:
	makedepf90 ${BFMSRC}/General/*.?90 \
		${BFMSRC}/PelB/*.?90 \
		${BFMSRC}/CO2/*.?90 \
		${BFMSRC}/Light/*.?90 \
		${BFMSRC}/Oxygen/*.?90 \
		${BFMSRC}/Forcing/*.?90 \
		${BFMSRC}/PelBen/ModuleSettling.f90 \
		${BFMSRC}/PelBen/Settling.F90 \
		${BFMSRC}/PelBen/BentoPelCoup.F90 \
		${BFMSHARE}/ResetFluxes.F90 \
		${BFMSHARE}/api_bfm.F90 \
		${BFMSHARE}/calcmean_bfm.F90 \
		${BFMSHARE}/init_cnps.F90 \
		${BFMSHARE}/init_var_bfm.F90 \
		${BFMSHARE}/netcdf_bfm.F90 \
		${BFMSHARE}/string_functions.F90 \
		${BFMNEMOSRC}/*.?90 > .depend
clean:
	rm -rf ${BFMSRC}/*/*.o ${BFMNEMOSRC}/*.o ${BFMSHARE}/*.o
	rm -rf ${BFMNEMOSRC}/*.mod
	rm -rf ${BFMSRC}/General/ModuleMem.F90
	rm -rf ${BFMSRC}/General/AllocateMem.F90
	rm -rf ${BFMSRC}/General/set_var_info_bfm.F90
	rm -rf ${BFMSRC}/include/INCLUDE.h


.SUFFIXES: 
.SUFFIXES: .f90 .F90

# rules
.DEFAULTS:

%.o: %.f90
		$(F_C) $(P_P) $(F_O) -c $< -o $@\
		|| { if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi ; false ; exit ; }
ifeq ($(FORTRAN_COMPILER),NECSX6)
		$(A_X) $(BFMLIB) $*.o  > /dev/null
		$(RM) $*.[of]
		if [ -f $*.o ] ; then mv $*.o $(NEMOLIBDIR)/oce ; fi
endif
		if [ -f $*.mod ] ; then mv $*.mod $(NEMOLIBDIR)/oce ; fi
		if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi
%.o: %.F90
		$(F_C) $(P_P) $(F_O) -c $< -o $@\
		|| { if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi ; false ; exit ; }
ifeq ($(FORTRAN_COMPILER),NECSX6)
		$(A_X) $(BFMLIB) $*.o  > /dev/null
		$(RM) $*.[of]
		if [ -f $*.o ] ; then mv $*.o $(NEMOLIBDIR)/oce ; fi
endif
		if [ -f $*.mod ] ; then mv $*.mod $(NEMOLIBDIR)/oce ; fi
		if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi



#-----------------------------------------------------------------------
# Copyright (C) 2006 - the BFM-team
#-----------------------------------------------------------------------
