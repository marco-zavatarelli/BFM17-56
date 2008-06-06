## Rules.make for BFM
## Modified after the GOTM Makefile (Karsten Bolding & Hans Burchard)

SHELL   = /bin/sh

## The compilation mode is obtained from $COMPILATION_MODE -
#### default production - else debug or profiling
ifndef COMPILATION_MODE
  compilation=production
else
  compilation=$(COMPILATION_MODE)
endif

FEATURES	=
FEATURE_LIBS	=
EXTRA_LIBS	=
INCDIRS		=
LDFLAGS		=

## Netcdf specifications
ifndef NETCDFINC 
  $(error The environment variable NETCDFINC is not defined!)
endif
ifndef NETCDFLIBDIR 
  $(error The environment variable NETCDFLIBDIR is not defined!)
endif
INCDIRS		+= -I$(NETCDFINC)
## NETCDFLIB	= $(NETCDFLIBNAME)
NETCDFLIB	= -lnetcdf
LDFLAGS		+= -L$(NETCDFLIBDIR)

##
## phony targets
##
.PHONY: clean realclean distclean dummy

## Default root directory BFM
# BFMDIR path must be given with an environmental variable
ifndef BFMDIR
  $(error The environment variable BFMDIR is not defined!)
endif

CPP = /lib/cpp

EXTRA_LIBS += $(NETCDFLIB)

## Directory related settings.
ifndef BINDIR
BINDIR	= $(BFMDIR)/bin
endif
ifndef LIBDIR
LIBDIR	= $(BFMDIR)/lib/$(FORTRAN_COMPILER)
endif
ifndef MODDIR
MODDIR	= $(BFMDIR)/modules/$(FORTRAN_COMPILER)
endif

INCDIRS	+= -I/usr/local/include -I$(BFMDIR)/include -I$(MODDIR)

## -------------------------------------------------------------------------
## BFM CONFIGURATION 
## -------------------------------------------------------------------------
# BFM source directories
BFMSRCDIR =$(BFMDIR)/src/BFM/General 
BFMSRCDIR +=$(BFMDIR)/src/BFM/PelB
BFMSRCDIR +=$(BFMDIR)/src/BFM/Oxygen
BFMSRCDIR +=$(BFMDIR)/src/BFM/Light
BFMSRCDIR +=$(BFMDIR)/src/BFM/Forcing
BFMSRCDIR +=$(BFMDIR)/src/BFM/PelBen
BFMSRCDIR +=$(BFMDIR)/src/share

# BFM include files and the library
BFMINCDIR = $(BFMDIR)/src/BFM/include
INCDIRS		+= -I$(BFMINCDIR)

# Option to compile without pointers (default with NECSX)
ifeq ($(FORTRAN_COMPILER),NECSX6)
  DEFINES += -DBFM_NOPOINTERS 
endif
# Option to compile without allocatable memory (default with GOTM)
## DEFINES += -DNOT_STANDALONE

# Pelagic CO2 flags (activates compilation and macros, false by default)
INCLUDE_PELCO2=true
ifeq ($(INCLUDE_PELCO2),true)
  DEFINES += -DINCLUDE_PELCO2
  BFMSRCDIR += $(BFMDIR)/src/BFM/CO2
endif
# Benthic ecosystem flags (activates compilation and macros, true by default)
INCLUDE_BEN = false
INCLUDE_BENCO2=false
INCLUDE_BENPROFILES=false
ifeq ($(INCLUDE_BEN),true)
  DEFINES += -DINCLUDE_BEN
  ifeq ($(INCLUDE_BENCO2),true)
    DEFINES += -DINCLUDE_BENCO2
    DEFINES += -DINCLUDE_PELCO2
    BFMSRCDIR += $(BFMDIR)/src/BFM/CO2
  endif
  ifeq ($(INCLUDE_BENPROFILES),true)
    DEFINES += -DINCLUDE_BENPROFILES
  endif
endif
# Sea-ice ecosystem (activates compilation and macros, false by default)
INCLUDE_SEAICE = false
ifeq ($(INCLUDE_SEAICE),true)
  # the benthic system must be activated because the sea-ice model is 
  # defined as a benthic layer. If the user wants to remove the 
  # benthic variables for memory reasons, it is suggested to use
  # the GlobalDefsBFM.model.seaice file in BFM/General/Configurations
  DEFINES += -DINCLUDE_BEN
  DEFINES += -DINCLUDE_SEAICE
  BFMSRCDIR += $(BFMDIR)/src/BFM/Seaice
endif
# Silt submodel (only for GOTM)
INCLUDE_SILT = false
ifeq ($(INCLUDE_SILT),true)
  DEFINES += -DINCLUDE_SILT
  BFMSRCDIR += $(BFMDIR)/src/BFM/Silt
endif

## -------------------------------------------------------------------------
## Normally no changes below this line
## -------------------------------------------------------------------------

## The Fortran compiler is determined from the EV FORTRAN_COMPILER - options 
## so far: Intel ifort, gfortran, NEC SXF90
## Default is none
ifndef FORTRAN_COMPILER
  $(error The environment variable FORTRAN_COMPILER is not defined!)
endif

include $(BFMDIR)/compilers/compiler.$(FORTRAN_COMPILER)

DEFINES += -DREAL_4B=$(REAL_4B)

## Sets options for debug compilation
ifeq ($(compilation),debug)
buildtype = _debug
DEFINES += -DDEBUG $(STATIC)
FLAGS   = $(DEBUG_FLAGS) 
endif

## Sets options for profiling compilation
ifeq ($(compilation),profiling)
buildtype = _prof
DEFINES += -DPROFILING $(STATIC)
FLAGS   = $(PROF_FLAGS) 
endif

## Sets options for production compilation
ifeq ($(compilation),production)
buildtype = _prod
DEFINES += -DPRODUCTION $(STATIC)
FLAGS   = $(PROD_FLAGS) 
endif

## For making the source code documentation.
PROTEX	= protex -b -n -s

.SUFFIXES:
.SUFFIXES: .F90 .f90

LINKDIR	= -L$(LIBDIR)

CPPFLAGS	= $(DEFINES) $(INCDIRS)
FFLAGS  	= $(DEFINES) $(FLAGS) $(MODULES) $(INCDIRS) $(EXTRAS)
F90FLAGS  	= $(FFLAGS)
LDFLAGS		+= $(FFLAGS) $(LINKDIR)

##
## Common rules
##
ifeq  ($(can_do_F90),true)
%.o: %.F90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
%.o: %.f90 
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
else
%.f90: %.F90
##	$(CPP) $(CPPFLAGS) $< -o $@
	$(F90_to_f90)
%.o: %.f90
	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@
endif
## Special option for GFORTRAN
ifeq ($(FORTRAN_COMPILER),GFORTRAN)
%.o: %.mod
endif

