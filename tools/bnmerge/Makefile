ifeq ($(DEBUG),yes)
	OPTDEBUG=-DDEBUG
endif

ifeq ($(PARALLEL),yes)
	MODE=openmp
	OPTPAR=-DPARALLEL
	include ../../compilers/ifort_athena_openmp.inc
else
	MODE=serial
	include ../../compilers/gfortran.inc
endif



.DEFAULT:
	-touch $@
all: bnmerge
bnmerge.o: src/bnmerge.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/bnmerge.F90 $(OPTPAR) $(OPTDEBUG)
create_outputfile_merged.o: src/$(MODE)/create_outputfile_merged.F90 mod_bnmerge.o
	$(FC) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/$(MODE)/create_outputfile_merged.F90 $(OPTPAR) $(OPTDEBUG)
merge_vars_merged.o: src/$(MODE)/merge_vars_merged.F90 mod_bnmerge.o
	$(FC) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/$(MODE)/merge_vars_merged.F90 $(OPTPAR) $(OPTDEBUG)
mod_bnmerge.o: src/mod_bnmerge.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/mod_bnmerge.F90 $(OPTPAR) $(OPTDEBUG)
read_input.o: src/read_input.F90 mod_bnmerge.o
	$(FC) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/read_input.F90 $(OPTPAR) $(OPTDEBUG)

SRC = src/${MODE}/create_outputfile_merged.F90 src/mod_bnmerge.F90 src/${MODE}/merge_vars_merged.F90 src/bnmerge.F90 src/read_input.F90
OBJ = create_outputfile_merged.o mod_bnmerge.o merge_vars_merged.o bnmerge.o read_input.o
clean: neat
	-rm -f runscript .cppdefs *.mod $(OBJ) bnmerge.x LinkerMap
neat:
	-rm -f $(TMPFILES)
bnmerge: $(OBJ) 
	$(LD) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(LDFLAGS) -o bnmerge.x $(OBJ) $(OPTPAR) $(OPTDEBUG)
