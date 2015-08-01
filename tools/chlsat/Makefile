# Makefile chlsat 

include ${BFMDIR}/compilers/ifort.inc

# CPPDEFS=-DDEBUG

.DEFAULT:
	-touch $@
all: chlsat
chlsat.o: src/chlsat.F90 mod_chlsat.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/chlsat.F90
mod_chlsat.o: src/mod_chlsat.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/mod_chlsat.F90
read_create.o: src/read_create.F90 mod_chlsat.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	src/read_create.F90
SRC = src/read_create.F90 src/mod_chlsat.F90 src/chlsat.F90
OBJ = read_create.o mod_chlsat.o chlsat.o
clean: neat
	-rm -f .cppdefs $(OBJ) *.mod chlsat.x
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
chlsat: $(OBJ) 
	$(LD) $(OBJ) -o chlsat.x  $(LDFLAGS)
