#NCDF_DIR     netcdf path
#NCDF_INC     netcdf include file 
#NCDF_LIB     netcdf library (also hdf5 with nc4)
#FC           Fortran compiler command
#FCFLAGS      Fortran compiler flags
#FFLAGS       Fortran 77 flags
#LD 	        Linker
#LDFLAGS      Linker flags 
#AR	          Assembler
#ARFLAGS      Assembler flags
#OTHERFLAGS   User defined additional flags

  NCDF_DIR=${NETCDF}
  NCDF_LIB=${NCDF_DIR}/lib -lnetcdf 
  NCDF_INC=${NCDF_DIR}/include
  FC=mpxlf90_r -cpp
  FCFLAGS=
  LD=mpxlf90_r
  AR=ar
  ARFLAGS=-rs
  OTHERFLAGS=
  FFLAGS=-O3 -qsuffix=f=f90  -qstrict -qalias=intptr -qsuppress=1501-245 -qarch=pwr6 -q64 -qrealsize=8 -qtune=pwr6 -qextname -qlargepage -qmaxmem=-1 -I${NCDF_INC}
  LDFLAGS=-I${NCDF_INC} -L${NCDF_LIB} -bloadmap:LinkerMap
