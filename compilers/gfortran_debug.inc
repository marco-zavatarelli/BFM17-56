# Include file for compilation with gfortran 4.2-4.6
# Note that the debugging option -g is default and can be removed
# you may also need the following netcdf libraries 
# -L/usr/local/lib -lhdf5_hl -lhdf5 -lz

  NCDF_DIR=${NETCDF}
  NCDF_LIB= -L${NCDF_DIR}/lib -lnetcdff -lnetcdf 
  NCDF_INC= -I${NCDF_DIR}/include

  FC=gfortran
  LD=gfortran

  AR=ar
  ARFLAGS=-r

  CPPFLAGS=
  FFLAGS= -ffpe-trap=invalid,zero,overflow -fbounds-check -Wall -Wuninitialized -g -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize ${NCDF_INC}
  LDFLAGS= ${NCDF_INC} ${NCDF_LIB} 
