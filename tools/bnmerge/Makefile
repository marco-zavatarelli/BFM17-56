# Makefile bnmerge

ifeq ($(DEBUG),yes)
	OPTDEBUG=-DDEBUG
endif

ifeq ($(PARALLEL),yes)
	COMPILER=ifort_openmp.inc
	OPTPAR= -DPARAL
else
	COMPILER=ifort.inc
endif

include ${BFMDIR}/compilers/$(COMPILER)

.DEFAULT:
	-touch $@
all: bnmerge
bnmerge.o: src/bnmerge.F90
	$(FC) $(FFLAGS) $(OPTPAR) $(OPTDEBUG) -c src/bnmerge.F90
create_output.o: src/create_output.F90 mod_bnmerge.o
	$(FC) $(FFLAGS) $(OPTPAR) $(OPTDEBUG) -c src/create_output.F90
merge_vars.o: src/merge_vars.F90 mod_bnmerge.o
	$(FC) $(FFLAGS) $(OPTPAR) $(OPTDEBUG) -c src/merge_vars.F90
mod_bnmerge.o: src/mod_bnmerge.F90
	$(FC) $(FFLAGS) $(OPTPAR) $(OPTDEBUG) -c src/mod_bnmerge.F90
read_input.o: src/read_input.F90 mod_bnmerge.o
	$(FC) $(FFLAGS) $(OPTPAR) $(OPTDEBUG) -c src/read_input.F90

OBJ = create_output.o mod_bnmerge.o merge_vars.o bnmerge.o read_input.o
clean: neat
	-rm -f runscript .cppdefs *.mod $(OBJ) bnmerge.x LinkerMap
neat:
	-rm -f $(TMPFILES)
bnmerge: $(OBJ) 
	$(LD) $(FFLAGS) $(OPTPAR) $(OPTDEBUG) $(OBJ) -o bnmerge.x $(LDFLAGS)
