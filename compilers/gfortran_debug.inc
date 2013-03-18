# Include file for compilation with gfortran 4.2-4.6
# Note that the debugging option -g is default and can be removed
# you may also need the following netcdf libraries 
# -L/usr/local/lib -lhdf5_hl -lhdf5 -lz
	FC=gfortran
	LD=gfortran
	AR=ar
	ARFLAGS=-r
        NCDF_LIB= -L/usr/local/lib -lnetcdff -lnetcdf 
        NCDF_INC= /usr/local/include 
	CPPFLAGS=
	FFLAGS= -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -Wall -Wuninitialized -g -fdefault-real-8 -funroll-all-loops -fcray-pointer -ftree-vectorize -I${NCDF_INC}
	LDFLAGS= -I${NCDF_INC} ${NCDF_LIB} 
