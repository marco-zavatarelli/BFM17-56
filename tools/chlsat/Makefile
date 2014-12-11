# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include ${BFMDIR}/compilers/gfortran.inc

# CPPDEFS=-DDEBUG

.DEFAULT:
	-touch $@
all: chlsat
chlsat.o: ./chlsat.F90 mod_chlsat.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./chlsat.F90
mod_chlsat.o: ./mod_chlsat.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./mod_chlsat.F90
read_create.o: ./read_create.F90 mod_chlsat.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	./read_create.F90
SRC = ./read_create.F90 ./mod_chlsat.F90 ./chlsat.F90
OBJ = read_create.o mod_chlsat.o chlsat.o
clean: neat
	-rm -f .cppdefs $(OBJ) chlsat
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
chlsat: $(OBJ) 
	$(LD) $(OBJ) -o chlsat  $(LDFLAGS)
