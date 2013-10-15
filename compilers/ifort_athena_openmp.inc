#NCDF_DIR     netcdf path
#NCDF_INC     netcdf include file
#NCDF_LIB     netcdf library (also hdf5 with nc4)
#MPI_INC      mpi include files
#MPI_LIB      mpi library
#FC           Fortran compiler command
#FCFLAGS      Fortran compiler flags
#FFLAGS       Fortran 77 flags
#LD           Linker
#LDFLAGS      Linker flags
#AR           Assembler
#ARFLAGS      Assembler flags
#OTHERFLAGS   User defined additional flags

  NCDF_DIR=${PNETCDF}
  NCDF_LIB=${NCDF_DIR}/lib -lpnetcdf
  NCDF_INC=${NCDF_DIR}/include
  FC = ifort
  LD = ifort
  AR = ar
  ARFLAGS = -r
  MPI_LIB = /opt/intel/impi/4.1.0/lib64 -lmpi -lmpich -lmpiif
  MPI_INC = /opt/intel/impi/4.1.0/include64
  CPPFLAGS=
  FFLAGS  =  -r8 -O3 -xHost -fp-model source -traceback -I${NCDF_INC} -I${MPI_INC} -openmp
  LDFLAGS = -I${NCDF_INC} -I${MPI_INC} -L${NCDF_LIB} -L${MPI_LIB} -openmp

#NOTE: you have to load module: NETCDF/parallel-netcdf-1.3.1
