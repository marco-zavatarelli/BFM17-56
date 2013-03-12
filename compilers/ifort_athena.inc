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

	FC=mpiifort
	LD=mpiifort
	AR=ar
	ARFLAGS=-r
	NCDF_LIB= -L /users/home/opt/netcdf/netcdf-4.2.1.1_netcdf3/lib -lnetcdf -lnetcdff
	NCDF_INC= /users/home/opt/netcdf/netcdf-4.2.1.1_netcdf3/include
	CPPFLAGS=
	FFLAGS= -r8 -O3 -no-prec-div -xHost -fp-model fast -traceback -I${NCDF_INC}
	LDFLAGS= -I${NCDF_INC} ${NCDF_LIB} 
