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

# compilation macros (compiler dependent, DO SOMETHING)
REAL_4B = real\(4\)

## NEMO modules (with absolute path)
F_O += -I$(NEMOLIBDIR)/oce
## BFM include files and pp macros
BFMINCDIR = $(BFMDIR)/src/BFM/include
F_O += -I$(BFMINCDIR) -I$(BFMDIR)/include -I$(NCDF_INC) -I$(NEMOWORKDIR) -I$(NEMOLIBDIR)/oce
P_P += -DBFM_NEMO  -DREAL_4B=$(REAL_4B) -fPIC -shared-intel -mcmodel=large


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
	${BFMLIB}(${BFMSRC}/PelBen/ModuleSettling.o)			\
	${BFMLIB}(${BFMSRC}/CO2/ModuleCO2.o)


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
	${BFMLIB}(${BFMSRC}/PelBen/Settling.o)				\
	${BFMLIB}(${BFMSRC}/CO2/CO2Calc.o)				\
	${BFMLIB}(${BFMSRC}/CO2/CO2Flux.o)				\
	${BFMLIB}(${BFMSRC}/CO2/Drtsafe.o)				\
	${BFMLIB}(${BFMSRC}/CO2/Ta_Iter_1.o)


.PHONY: all depend clean

# Include the dependency-list created by makedepf90 below
include .depend

all .depend: depend $(BFM_PELMOD) $(BFMNEMO_OBJ) $(BFM_PELOBJ) 
#all : $(BFM_PELMOD) $(BFMNEMO_OBJ) $(BFM_PELOBJ) 

$(BFM_MOD) : $(BFMSRC)/General/ModuleMem.F90

${BFMSRC}/General/ModuleMem.F90 : $(BFMSRC)/General/GlobalDefsBFM.model
	${BFMSRC}/scripts/GenerateGlobalBFMF90Code  -read ${BFMSRC}/General/GlobalDefsBFM.model \
                -from ${BFMSRC}/proto -to ${BFMSRC}/General -actions statemem allocmem netcdfmem \
                -to ${BFMSRC}/include -actions headermem
depend:
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

.SUFFIXES: 
.SUFFIXES: .f90 .F90

# rules
.DEFAULTS:

%.o: %.f90
		$(F_C) $(P_P) $(F_O) -c $< -o $@\
		|| { if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi ; false ; exit ; }
		-@#$(A_C) $(BFMLIB) $*.o  > /dev/null
		-@#$(RM) $*.[of]
		-@if [ -f $*.mod ] ; then mv $*.mod $(NEMOLIBDIR)/oce ; fi
		-@#if [ -f $*.o ] ; then mv $*.o $(NEMOLIBDIR)/oce ; fi
		-@if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi
%.o: %.F90
		$(F_C) $(P_P) $(F_O) -c $< -o $@\
		|| { if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi ; false ; exit ; }
		-@#$(A_C) $(BFMLIB) $*.o  > /dev/null
		-@#$(RM) $*.[of]
		-@if [ -f $*.mod ] ; then mv $*.mod $(NEMOLIBDIR)/oce ; fi
		-@#if [ -f $*.o ] ; then mv $*.o $(NEMOLIBDIR)/oce ; fi
		-@if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi



#-----------------------------------------------------------------------
# Copyright (C) 2006 - the BFM-team
#-----------------------------------------------------------------------
