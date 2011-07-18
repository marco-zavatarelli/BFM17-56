	FC=gfortran
	LD=gfortran
	AR=ar
	ARFLAGS=-r
        NCDF_LIB= -L /usr/lib -lnetcdff -lnetcdf
        NCDF_INC= /usr/include
	CPPFLAGS=
	FFLAGS= -ffree-form -Wall -std=f95 -pedantic -O3 -ftree-vectorize -ftree-vectorizer-verbose=5 -I${NCDF_INC}
	LDFLAGS= -I${NCDF_INC} ${NCDF_LIB} 
