# Makefile bnremap

include ${BFMDIR}/compilers/ifort.inc

# CPPDEFS=-DDEBUG

.DEFAULT:
	-touch $@
all: bnremap
bnremap.o: src/bnremap.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/bnremap.F90
create_outputfile.o: src/create_outputfile.F90 mod_bnremap.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/create_outputfile.F90
remap_vars.o: src/remap_vars.F90 mod_bnremap.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/remap_vars.F90
mod_bnremap.o: src/mod_bnremap.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/mod_bnremap.F90
read_input.o: src/read_input.F90 mod_bnremap.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/read_input.F90
SRC = src/create_outputfile.F90 src/mod_bnremap.F90 src/remap_vars.F90 src/bnremap.F90 src/read_input.F90
OBJ = create_outputfile.o mod_bnremap.o remap_vars.o bnremap.o read_input.o
clean: neat
	-rm -f .cppdefs *.mod $(OBJ) bnremap.x
neat:
	-rm -f $(TMPFILES)
bnremap: $(OBJ) 
	$(LD) $(OBJ) -o bnremap.x  $(LDFLAGS)
