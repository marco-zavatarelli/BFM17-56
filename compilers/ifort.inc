#NCDF_INC     netcdf include file
#NCDF_LIB     netcdf library (also hdf5 with nc4)
#MPIINC       mpi include files
#MPILIB       mpi library
#FC           Fortran compiler command
#FCFLAGS      Fortran compiler flags
#FFLAGS       Fortran 77 flags
#LD           Linker
#LDFLAGS      Linker flags
#AR           Assembler
#ARFLAGS      Assembler flags
#OTHERFLAGS   User defined additional flags

	FC       = ifort
	LD       = ifort

	NCDF_LIB  = -L${NETCDF}/lib -lnetcdff -lnetcdf
	NCDF_INC  = -I${NETCDF}/include
  OTHER_LIB = 

  FFLAGS   = -r8 -O3 -xHost -fp-model source -traceback ${NCDF_INC}
	LDFLAGS  = ${NCDF_LIB} ${OTHER_LIB}

#NOTE: you have to load modules: 
#      - NETCDF/netcdf-4.3