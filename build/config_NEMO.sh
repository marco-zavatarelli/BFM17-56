#!/bin/sh
## Configuration file for Pelagic BFM STANDALONE
#  This script creates a directory BLD_STANDALONE with the appropriate 
#  makefile

#  Currently available macros (cppdefs) are:
#  INCLUDE_SILT
#  INCLUDE_PELCO2, INCLUDE_BENCO2
#  INCLUDE_BEN, INCLUDE_BENPROFILES
#  INCLUDE_SEAICE
#  INCLUDE_DIAG3D, INCLUDE_DIAG2D

#  Warning: still not working for benthic BFM
#           don't use DIAG with D1SOURCE and ONESOURCE

#----------------- User configuration -----------------
myGlobalDef="GlobalDefsBFM.model.standard"
cppdefs="-DINCLUDE_PELCO2"

RELEASE="no"

BLDDIR="STD_MemLayout"
#----------------- User configuration -----------------
# set -xv
cp="cp"
mv="mv"

if [ "${BFMDIR}" = "" ]; then
   echo "Environmental variable BFMDIR not defined!"
   exit 0
fi

if [ ! -d ${BFMDIR}/build/${BLDDIR} ]; then
  mkdir ${BFMDIR}/build/${BLDDIR}
fi

cd ${BLDDIR}
  # Link to standard file
ln -sf ${BFMDIR}/build/Configurations/${myGlobalDef} ${BFMDIR}/src/BFM/General/GlobalDefsBFM.model

# generate BFM Memory Layout files
${BFMDIR}/build/scripts/GenerateGlobalBFMF90Code ${cppdefs} \
          -read ${BFMDIR}/build/Configurations/${myGlobalDef} \
          -from ${BFMDIR}/src/BFM/proto \
          -to ${BFMDIR}/build/${BLDDIR} -actions statemem allocmem netcdfmem \
          -to ${BFMDIR}/build/${BLDDIR} -actions headermem \
          -to ${BFMDIR}/build/${BLDDIR} -actions initmem

echo "Memory Layout generated in local folder: ${BLDDIR}."

# If a RELEASE Configuration copy files to final destination
if [ ${RELEASE} = "yes" ]; then
  # Link to standard file 
  ln -sf ${BFMDIR}/build/Configurations/${myGlobalDef} ${BFMDIR}/src/BFM/General/GlobalDefsBFM.model

  ${cp} ${BFMDIR}/build/${BLDDIR}/*.F90 ${BFMDIR}/src/BFM/General
  ${mv} ${BFMDIR}/src/BFM/General/init_var_bfm.F90 ${BFMDIR}/src/share
  ${cp} ${BFMDIR}/build/${BLDDIR}/INCLUDE.h ${BFMDIR}/src/BFM/include

echo "Files copied to target folders."
 
fi

echo "DONEEEE !!!!"
