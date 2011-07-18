#!/bin/sh
## Configuration file for Pelagic BFM STANDALONE
#  This script creates a directory BLD_STANDALONE with the appropriate 
#  makefile
#  Current available macros (cppdefs) are:
#  INCLUDE_PELCO2
#  INCLUDE_SEAICE
#  INCLUDE_DIAG3D
#  Warning: still not working for benthic BFM

#----------------- User configuration -----------------
archfile="${BFMDIR}/compilers/xlf90.inc"
cppdefs="-DBFM_STANDALONE -DINCLUDE_PELCO2 -DINCLUDE_DIAG3D"
exe=${BFMDIR}/bin/bfm_standalone.x
#----------------- User configuration -----------------
#set -xv

if [ "${BFMDIR}" = "" ]; then
   echo "Environmental variable BFMDIR not defined!"
   exit 0
fi
BLDDIR="./BLD_STANDALONE"
MKMF="${BFMDIR}/bin/mkmf"
oflags="-I${BFMDIR}/include -I${BFMDIR}/src/BFM/include"

if [ ! -d ${BLDDIR} ]; then
  mkdir ${BLDDIR}
fi

cd ${BLDDIR}

# generate BFM files
${BFMDIR}/src/BFM/scripts/GenerateGlobalBFMF90Code   ${cppdefs} \
        -read ${BFMDIR}/src/BFM/General/GlobalDefsBFM.model \
        -from ${BFMDIR}/src/BFM/proto -to ${BFMDIR}/src/BFM/General \
        -actions statemem allocmem netcdfmem \
        -to ${BFMDIR}/src/BFM/include -actions headermem

find ${BFMDIR}/src/BFM/General -name "*.?90" -print > BFM.lst
find ${BFMDIR}/src/standalone -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/share -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/PelB -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/PelBen -name "*.?90" -print >> BFM.lst
#find ${BFMDIR}/src/BFM/Ben -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/Light -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/Oxygen -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/CO2 -name "*.?90" -print >> BFM.lst

${MKMF} -c "${cppdefs}" -o "${oflags}" -t ${archfile} -p ${exe} BFM.lst
