# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include ../../compilers/xlf90.inc


.DEFAULT:
	-touch $@
all: bnmerge
bnmerge.o: ./bnmerge.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./bnmerge.F90
create_outputfile.o: ./create_outputfile.F90 mod_bnmerge.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./create_outputfile.F90
merge_vars.o: ./merge_vars.F90 mod_bnmerge.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./merge_vars.F90
mod_bnmerge.o: ./mod_bnmerge.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./mod_bnmerge.F90
read_input.o: ./read_input.F90 mod_bnmerge.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./read_input.F90
SRC = ./create_outputfile.F90 ./mod_bnmerge.F90 ./merge_vars.F90 ./bnmerge.F90 ./read_input.F90
OBJ = create_outputfile.o mod_bnmerge.o merge_vars.o bnmerge.o read_input.o
clean: neat
	-rm -f .cppdefs *.mod $(OBJ) bnmerge
neat:
	-rm -f $(TMPFILES)
bnmerge: $(OBJ) 
	$(LD) $(OBJ) -o bnmerge  $(LDFLAGS)
