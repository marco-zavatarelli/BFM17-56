	FC=gfortran
	LD=gfortran
	AR=ar
	ARFLAGS=-r
        NCDF_LIB= -L /opt/netcdf/lib -lnetcdf
        NCDF_INC= /opt/netcdf/include
	CPPFLAGS=
	FFLAGS= -ffree-form -Wall -std=f95 -pedantic -O3 -ftree-vectorize -ftree-vectorizer-verbose=5
	LDFLAGS= -I${NCDF_INC} ${NCDF_LIB} 
