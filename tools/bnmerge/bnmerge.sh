#!/bin/bash

export OMP_NUM_THREADS=16
LIST=bnmerge.nml
LOG_DIR=./
MODE=openmp
#DEBUG=-DDEBUG

gmake clean
gmake DEBUG=${DEBUG} MODE=${MODE}

# cat > runscript <<EOF
#     #! /bin/sh 

#     #BSUB -a poe
#     #BSUB -J bnmerge         # Name of the job.
#     #BSUB -o ${LOG_DIR}/bnmerge%J.out  # Appends std output to file %J.out.
#     #BSUB -e ${LOG_DIR}/bnmerge%J.err  # Appends std error to file %J.out.
#     #BSUB -P bnmerge
#     #BSUB -q poe_short    # queue
#     #BSUB -n ${OMP_NUM_THREADS}            # Number of CPUs
#     ##BSUB -x 
#     ###BSUB -R "span[ptile=32]"

#     VERBOSE=_VERBOSE_
#     if [ $VERBOSE ]; then set -exv; fi

#     export MP_WAIT_MODE=poll
#     export MP_POLLING_INTERVAL=30000000
#     export MP_SHARED_MEMORY=yes
#     export MP_EUILIB=us
#     export MP_EUIDEVICE=sn_all
#     export LDR_CNTRL=TEXTPSIZE=64K@STACKPSIZE=64K@DATAPSIZE=64K
#     export MP_TASK_AFFINITY=core

#     # Launch the model

#     time mpirun.lsf ./bnmerge -f ${LIST}

#     echo " bnmerge DONEEE!!!"

# EOF

# if [ "${MODE}" == 'openmp' ]; then
#     bsub < runscript
#     while `true` ; do
#         echo "Waiting..."
#         sleep 10
#         njobs=`bjobs | grep bnmerge | wc -l`
#         if [ ${njobs} -eq 0 ]; then
#             break
#         fi
#     done
# else
    ./bnmerge -f ${LIST}
# fi
