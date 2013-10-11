#!/bin/bash

set -ex 

# set environment
module load NETCDF/parallel-netcdf-1.3.1

export OMP_NUM_THREADS=16
LIST=bnmerge.nml
LOG_DIR=./


gmake clean
gmake


if [ "${PARALLEL}" == 'yes' ]; then
    #QUEUE='poe_short_smt'
    QUEUE='serial_30min'
else
    QUEUE='serial_30min'
fi


cat > runscript <<EOF
    #! /bin/sh 

    #BSUB -J bnmerge         # Name of the job.
    #BSUB -o ${LOG_DIR}/bnmerge%J.out  # Appends std output to file %J.out.
    #BSUB -e ${LOG_DIR}/bnmerge%J.err  # Appends std error to file %J.out.
    #BSUB -P bnmerge
    #BSUB -q ${QUEUE}    # queue
    #BSUB -n ${OMP_NUM_THREADS}            # Number of CPUs
    #BSUB -x 
    #BSUB -R "span[ptile=16]"

    VERBOSE=_VERBOSE_
    if [ $VERBOSE ]; then set -exv; fi

    export MP_TASK_AFFINITY=core

    # Launch the model

    time ./bnmerge.x -f ${LIST}

    echo " bnmerge DONEEE!!!"

EOF


bsub < runscript
while `true` ; do
    echo "Waiting..."
    sleep 10
    njobs=`bjobs | grep bnmerge | wc -l`
    if [ ${njobs} -eq 0 ]; then
        break
    fi
done
