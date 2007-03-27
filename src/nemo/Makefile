#
# Makefile to build the BFM model coupled with NEMO
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

## BFM include files and pp macros
BFMINCDIR = $(BFMDIR)/src/BFM/include
F_O += -I$(BFMINCDIR) -I$(BFMDIR)/include -I$(NCDF_INC) -I$(NEMOWORKDIR)
P_P += -DBFM_NEMO  -DREAL_4B=$(REAL_4B)
##P_P += -DNOT_STANDALONE -DBFM_NOPOINTERS 


BFM_MOD = \
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
	${BFMLIB}(${BFMSRC}/PelBen/ModuleControlBenPartNutrientBuffers.o)	\
	${BFMLIB}(${BFMSRC}/Ben/ModuleBenOrganism.o)			\
	${BFMLIB}(${BFMSRC}/Ben/ModuleBenBac.o)				\
	${BFMLIB}(${BFMSRC}/Ben/ModuleFilterFeeder.o)			\
	${BFMLIB}(${BFMSRC}/Ben/ModuleBioturbation.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenAmmonium.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenAnoxic.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenDenitriDepth.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenNitrate.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenNutConstants.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenNutType.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenNutInterface.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenNutVariables.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenOxygen.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenPhosphate.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenQ1Transport.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenSilica.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenthicReturn1.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/ModuleBenthicReturn2.o)			 

BFMNEMO_OBJ   = \
	${BFMLIB}(${BFMSHARE}/api_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/netcdf_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/ResetFluxes.o)		\
	${BFMLIB}(${BFMSHARE}/calcmean_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/string_functions.o)	\
	${BFMLIB}(${BFMSHARE}/init_cnps.o)		\
	${BFMLIB}(${BFMSHARE}/init_var_bfm.o)		\
	${BFMLIB}(${BFMSHARE}/init_benthic_bfm.o)	\
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

BFM_OBJ = \
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
	${BFMLIB}(${BFMSRC}/Ben/ResetTotMassVar.o)			\
	${BFMLIB}(${BFMSRC}/PelBen/BentoPelCoup.o)			\
	${BFMLIB}(${BFMSRC}/PelBen/PelForcingForBen.o)			\
	${BFMLIB}(${BFMSRC}/PelBen/Sedimentation.o)			\
	${BFMLIB}(${BFMSRC}/PelBen/Settling.o)				\
	${BFMLIB}(${BFMSRC}/PelBen/ControlBenPartNutrientBuffers.o)	\
	${BFMLIB}(${BFMSRC}/Ben/BenBac.o)				\
	${BFMLIB}(${BFMSRC}/Ben/BenOrganism.o)				\
	${BFMLIB}(${BFMSRC}/Ben/BenthicSystem.o)			\
	${BFMLIB}(${BFMSRC}/Ben/Bioturbation.o)				\
	${BFMLIB}(${BFMSRC}/Ben/FilterFeeder.o)				\
	${BFMLIB}(${BFMSRC}/Ben/ResetTotMassVar.o)			\
	${BFMLIB}(${BFMSRC}/Bennut/BenAmmonium.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenAnoxic.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenDenitriDepth.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenNitrate.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenNitrogenShifting.o) 		\
	${BFMLIB}(${BFMSRC}/Bennut/BenOxygen.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenPhosphate.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenQ1Transport.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenSilica.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenthicNutrient2.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenthicNutrient3.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenthicReturn1.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/BenthicReturn2.o) 			\
        ${BFMLIB}(${BFMSRC}/Bennut/InitBenthicNutrient.o)		\
        ${BFMLIB}(${BFMSRC}/Bennut/InitBenthicNutrient3.o)		\
	${BFMLIB}(${BFMSRC}/Bennut/CalculateFromSet.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/CalculateSet.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/CalculateShift.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/CalculateTau.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/CompleteSet.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/InitializeSet.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/DefineSet.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/GetInfoFromSet.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/PrintSet.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/bess_exp.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/bessi0.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/bessi1.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/bessk0.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/bessk1.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/calculate_equation.o) 		\
	${BFMLIB}(${BFMSRC}/Bennut/calculatelayer.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/funcalc.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/input_para.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/kfind.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/lubksb.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/ludcmp.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/manage_coeff.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/noutput.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/qgaus_exp.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/re_store.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/set_max_sing.o) 			\
	${BFMLIB}(${BFMSRC}/Bennut/svbksb.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/svdcmp.o) 				\
	${BFMLIB}(${BFMSRC}/Bennut/transfer.o) 

.PHONY: all depend clean

# Include the dependency-list created by makedepf90 below
include .depend

all .depend: depend $(BFM_MOD) $(BFMNEMO_OBJ) $(BFM_OBJ) 
#all : $(BFM_MOD) $(BFMNEMO_OBJ) $(BFM_OBJ) 

$(BFM_MOD) : $(BFMSRC)/General/ModuleMem.F90

${BFMSRC}/General/ModuleMem.F90 : $(BFMSRC)/General/GlobalDefsBFM.model
	${BFMSRC}/scripts/GenerateGlobalBFMF90Code  -read ${BFMSRC}/General/GlobalDefsBFM.model \
                -from ${BFMSRC}/proto -to ${BFMSRC}/General -actions statemem allocmem netcdfmem \
                -to ${BFMSRC}/include -actions headermem
depend:
	makedepf90 ${BFMSRC}/*/*.?90 ${BFMNEMOSRC}/*.?90 ${BFMSHARE}/*.?90 > .depend
clean:
	rm -rf ${BFMSRC}/*/*.o ${BFMNEMOSRC}/*.o ${BFMSHARE}/*.o

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
		-@if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi
%.o: %.F90
		$(F_C) $(P_P) $(F_O) -c $< -o $@\
		|| { if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi ; false ; exit ; }
		-@#$(A_C) $(BFMLIB) $*.o  > /dev/null
		-@#$(RM) $*.[of]
		-@if [ -f $*.mod ] ; then mv $*.mod $(NEMOLIBDIR)/oce ; fi
		-@if [ -f $(PREF)$*.L ] ; then mv $(PREF)$*.L $(NEMOTMPDIR) ; fi



#-----------------------------------------------------------------------
# Copyright (C) 2006 - the BFM-team
#-----------------------------------------------------------------------
