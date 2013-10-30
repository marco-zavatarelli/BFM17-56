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
  NCDF_LIB= -L${NCDF_DIR}/lib -lnetcdf
  NCDF_INC= -I${NCDF_DIR}/include

  FC=xlf90_r
  LD=xlf90_r 

  AR=ar
  ARFLAGS=-r

  CPPFLAGS=
  FCFLAGS -cpp
  FFLAGS= -cpp -O3 -qstrict -qsuffix=f=f90  -qarch=pwr6 -q64 -qrealsize=8 -qtune=pwr6 -qextname=flush ${NCDF_INC}
  LDFLAGS= ${NCDF_INC} ${NCDF_LIB} -bloadmap:LinkerMap
