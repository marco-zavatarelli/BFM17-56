#NCDF_DIR     netcdf path
#NCDF_INC     netcdf include file
#NCDF_LIB     netcdf library (also hdf5 with nc4)
#FC           Fortran compiler command
#FCFLAGS      Fortran compiler flags
#FFLAGS       Fortran 77 flags
#LD           Linker
#LDFLAGS      Linker flags
#AR           Assembler
#ARFLAGS      Assembler flags
#OTHERFLAGS   User defined additional flags

  NCDF_DIR=${NETCDF}
  NCDF_LIB=${NCDF_DIR}/lib -lnetcdff -lnetcdf
  NCDF_INC=${NCDF_DIR}/include
  FC=ifort
  LD=ifort
  AR=ar
  ARFLAGS=-r
  CPPFLAGS=
  FFLAGS= -O3 -I${NCDF_INC}
  LDFLAGS= -I${NCDF_INC} -L${NCDF_LIB} 
