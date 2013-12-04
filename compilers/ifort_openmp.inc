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
#NOTE: in athena you have to load modules: -HDF5/hdf5-1.8.11_parallel -NETCDF/parallel-netcdf-1.3.1 -NETCDF/netcdf-4.3_parallel

  NCDF_DIR=${NETCDF}
  NCDF_LIB= -L${NETCDF}/lib -lnetcdff -lnetcdf -L${PNETCDF}/lib -lpnetcdf
  NCDF_INC= -I${PNETCDF}/include -I${NETCDF}/include
  HDF5=/users/home/opt/hdf5/hdf5-1.8.11_parallel
  HDF5_LIB= -L${HDF5}/lib -lhdf5_hl -lhdf5
  HDF5_INC= -I${HDF5}/include

  FC=mpiifort
  LD=mpiifort

  CPPFLAGS=
  OTHER_LIB= -lz -lcurl -lgpfs
  FFLAGS= -r8 -O3 -xHost -fp-model source -traceback -openmp ${NCDF_INC} ${HDF5_INC}
	LDFLAGS= ${NCDF_LIB} ${HDF5_LIB} ${OTHER_LIB}
