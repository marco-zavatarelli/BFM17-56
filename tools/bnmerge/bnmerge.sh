#!/bin/bash
#set -ex 

# set environment

BNMERGE_LST=bnmerge.nml
LOG_DIR=.
BNMERGE_EXE=${BFMDIR}/tools/bnmerge/bnmerge.x

if [ "${PARALLEL}" == 'yes' ]; then
    export OMP_NUM_THREADS=128
    module load INTEL/intel_xe_2013 HDF5/hdf5-1.8.11_parallel NETCDF/netcdf-4.3_parallel NETCDF/parallel-netcdf-1.3.1
    QUEUE='poe_short'
else
    export OMP_NUM_THREADS=1
    module load INTEL/intel_xe_2013 NETCDF/netcdf-4.3
    QUEUE='serial_30min'
fi


# create runscript

cat > runscript <<EOF
    #! /bin/sh 

    #BSUB -J bnmerge         # Name of the job.
    #BSUB -o ${LOG_DIR}/bnmerge%J.out  # Appends std output to file %J.out.
    #BSUB -e ${LOG_DIR}/bnmerge%J.err  # Appends std error to file %J.out.
    #BSUB -P bnmerge
    #BSUB -q ${QUEUE}    # queue
    #BSUB -n ${OMP_NUM_THREADS}            # Number of CPUs
    #BSUB -x 
    #BSUB -R "span[ptile=32]"

    if [ ${DEBUG} ]; then set -exv; fi

    export MP_TASK_AFFINITY=core

    # Launch the model
    ${BNMERGE_EXE} -f ${BNMERGE_LST}

    echo " bnmerge DONEEE!!!"

EOF

# execute runscript and wait for output
echo "Move \"runscript\" to your experiment folder and execute: \"bsub < runscript\""
