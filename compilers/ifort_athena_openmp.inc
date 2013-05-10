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

	FC = ifort
	LD = ifort
	AR      = ar
	ARFLAGS = -r
	NCDF_LIB = -L${PNETCDF}/lib -lpnetcdf
	NCDF_INC = -I${PNETCDF}/include
  MPI_LIB =  -L/opt/intel/impi/4.1.0/lib64 -lmpi -lmpich -lmpiif
  MPI_INC =  -I/opt/intel/impi/4.1.0/include64
	CPPFLAGS=
  FFLAGS  =  -r8 -O3 -no-prec-div -xHost -fp-model fast -traceback ${NCDF_INC} ${MPI_INC} -openmp
	LDFLAGS = ${NCDF_INC} ${MPI_INC} ${NCDF_LIB} ${MPI_LIB} -openmp

#NOTE: you have to load module: NETCDF/parallel-netcdf-1.3.1