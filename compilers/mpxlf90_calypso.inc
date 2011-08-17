#NCDF_INC     netcdf include file 
#NCDF_LIB     netcdf library (also hdf5 with nc4)
#FC           Fortran compiler command
#FCFLAGS      Fortran compiler flags
#FFLAGS       Fortran 77 flags
#LD 	      Linker
#LDFLAGS      Linker flags 
#AR	      Assembler
#ARFLAGS      Assembler flags
#OTHERFLAGS   User defined additional flags



NCDF_INC=/usr/local/netcdf-3.6.3/include 
NCDF_LIB=-L /usr/local/netcdf-3.6.3/lib -lnetcdf 
FC=mpxlf90_r -cpp
FCFLAGS=
FFLAGS=-O3 -qsuffix=f=f90  -qstrict -qalias=intptr -qsuppress=1501-245 -qarch=pwr6 -q64 -qrealsize=8 -qtune=pwr6 -qextname -qlargepage -qmaxmem=-1 -I${NCDF_INC}
LD=mpxlf90_r
LDFLAGS=-I${NCDF_INC} ${NCDF_LIB} ${OTHERLIBS} -bloadmap:LinkerMap  
AR=ar
ARFLAGS=-rs
OTHERFLAGS=
