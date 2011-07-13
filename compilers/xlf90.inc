	FC=xlf90_r -cpp
	LD=xlf90_r 
	AR=ar
	ARFLAGS=-r
        NCDF_LIB= -L /usr/local/netcdf-3.6.3/lib -lnetcdf
        NCDF_INC= /usr/local/netcdf-3.6.3/include
	CPPFLAGS=
	FFLAGS= -O3 -qstrict -qalias=intptr -qsuffix=f=f90  -qarch=pwr6 -q64 -qrealsize=8 -qtune=pwr6 -qextname=flush -I${NCDF_INC}
	LDFLAGS= -I${NCDF_INC} ${NCDF_LIB} 
