ifeq ($(PARALLEL),yes)
	MODE=openmp
	include ../../compilers/ifort_athena_openmp.inc
else
	MODE=concurrent
	include ../../compilers/ifort_athena.inc
endif


.DEFAULT:
	-touch $@
all: bnmerge
bnmerge.o: src/$(MODE)/bnmerge.F90
	$(FC) $(DEBUG) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/$(MODE)/bnmerge.F90
create_outputfile.o: src/$(MODE)/create_outputfile.F90 mod_bnmerge.o
	$(FC) $(DEBUG) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/$(MODE)/create_outputfile.F90
merge_vars.o: src/$(MODE)/merge_vars.F90 mod_bnmerge.o
	$(FC) $(DEBUG) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/$(MODE)/merge_vars.F90
mod_bnmerge.o: src/$(MODE)/mod_bnmerge.F90
	$(FC) $(DEBUG) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/$(MODE)/mod_bnmerge.F90
read_input.o: src/read_input.F90 mod_bnmerge.o
	$(FC) $(DEBUG) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/read_input.F90

SRC = src/${MODE}/create_outputfile.F90 src/${MODE}/mod_bnmerge.F90 src/${MODE}/merge_vars.F90 src/${MODE}/bnmerge.F90 src/read_input.F90
OBJ = create_outputfile.o mod_bnmerge.o merge_vars.o bnmerge.o read_input.o
clean: neat
	-rm -f runscript .cppdefs *.mod $(OBJ) bnmerge LinkerMap
neat:
	-rm -f $(TMPFILES)
bnmerge: $(OBJ) 
	$(LD) $(OBJ) -o bnmerge  $(LDFLAGS)
