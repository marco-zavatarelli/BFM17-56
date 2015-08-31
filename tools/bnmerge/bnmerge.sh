#!/bin/bash
#set -ex 

# set environment
BNMERGE_LST=bnmerge.nml
LOG_DIR=.
BNMERGE_EXE=${BFMDIR}/tools/bnmerge/bnmerge.x
if [ "${PARALLEL}" == 'yes' ]; then
    LIBS="INTEL/intel_xe_2013 HDF5/hdf5-1.8.11_parallel NETCDF/netcdf-4.3_parallel NETCDF/parallel-netcdf-1.3.1"
    export OMP_NUM_THREADS=16
    QUEUE='poe_short'
else
    LIBS="INTEL/intel_xe_2013 NETCDF/netcdf-4.3"
    export OMP_NUM_THREADS=1
    QUEUE='serial_30min'
fi

# compile
module load ${LIBS}
gmake clean
gmake


# create runscript
cat > runscript <<EOF
#! /bin/sh 
#BSUB -J bnmerge                   # name of the job.
#BSUB -o ${LOG_DIR}/bnmerge%J.out  # appends std output to file %J.out.
#BSUB -e ${LOG_DIR}/bnmerge%J.err  # appends std error to file %J.out.
#BSUB -P bnmerge                   # project name
#BSUB -q ${QUEUE}                  # queue
#BSUB -n 1                         # Number of CPUs
##BSUB -x                           # exclusive host mode
##BSUB -R "span[ptile=16]"          # use only nodes with 16 cores

if [ ${DEBUG} ]; then set -exv; fi

export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export MP_TASK_AFFINITY=core

# Launch the model
${BNMERGE_EXE} -f ${BNMERGE_LST}

echo " bnmerge DONEEE!!!"
EOF

# execute runscript and wait for output
echo "--------------"
echo "Be sure you have these libs loaded in your experiment folder: "
echo "    ${LIBS}"
echo "Move \"runscript\" to your experiment folder and execute: "
echo "    bsub < runscript"
