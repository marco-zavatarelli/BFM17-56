# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

# compile for serial queue
include ../../compilers/xlf90_calypso_serial.inc

# compile for parallel queue
# include ../../compilers/mpxlf90_calypso.inc

.DEFAULT:
	-touch $@
all: bnremap
bnremap.o: ./bnremap.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./bnremap.F90
create_outputfile.o: ./create_outputfile.F90 mod_bnremap.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./create_outputfile.F90
remap_vars.o: ./remap_vars.F90 mod_bnremap.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./remap_vars.F90
mod_bnremap.o: ./mod_bnremap.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./mod_bnremap.F90
read_input.o: ./read_input.F90 mod_bnremap.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./read_input.F90
SRC = ./create_outputfile.F90 ./mod_bnremap.F90 ./remap_vars.F90 ./bnremap.F90 ./read_input.F90
OBJ = create_outputfile.o mod_bnremap.o remap_vars.o bnremap.o read_input.o
clean: neat
	-rm -f .cppdefs *.mod $(OBJ) bnremap
neat:
	-rm -f $(TMPFILES)
bnremap: $(OBJ) 
	$(LD) $(OBJ) -o bnremap  $(LDFLAGS)
