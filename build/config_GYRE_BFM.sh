#!/bin/sh
#
## Configuration for the BFM-NEMO coupling in the PELAGOS configuration. 
#
#  This script creates a directory $BLDDIR with the Memory Layout files and the FCM include for the coupling. 
#  The COMPILE flag starts the compilation using the makenemo tool.
#  Requires the environmental variables BFMDIR and NEMODIR to be set and pointing to the
#  root directories of BFM and NEMO 
# 

#  Currently available macros (cppdefs) are:
#  INCLUDE_PELFE                          : use Iron component to the pelagic system
#  INCLUDE_PELCO2, INCLUDE_BENCO2         : activate Carbonate System 
#  INCLUDE_BEN, INCLUDE_BENPROFILES       : Add Benthic compartment
#  INCLUDE_SILT                           : use Silt component
#  INCLUDE_SEAICE                         : activate SeaIce Ecology 
#  INCLUDE_DIAG3D, INCLUDE_DIAG2D         : additional diagnostics available for output
#  BFM_PARALLEL                           : used to run BFM with MPP

#  Warnings
# 1. Still not working for benthic BFM don't use DIAG with D1SOURCE and ONESOURCE
# 2. Using the key DEBUG will add more output information

# Author: Tomas Lovato (CMCC)

#----------------- BEGIN User configuration -----------------
# myGlobalDef   : file used by tcl script (GenerateGlobalBFMF90Code) 
#                 to design the BFM Memory Layout
# cppdefs       : keys used to configure the model
# CONFIG        : NEMO Configuration to be built
# BLDDIR        : Local folder containing the generated BFM Memory Layout files
# COMPILE       : set "yes" to compile NEMO with makenemo
# ARCH          : NEMO specific architecture file (get them with makenemo -h all)
# NPROC         : number of procs for compilation
# -----------------------------------------------------
myGlobalDef="GlobalDefsBFM.model.standard"
cppdefs="-DBFM_PARALLEL -DINCLUDE_PELCO2 -DINCLUDE_DIAG3D" 
CONFIG="GYRE_BFM"
BLDDIR="STD_${CONFIG}"
COMPILE="yes"
ARCH="PW6_calypso"
NPROC=4
MODE=${1}
#----------------- END User configuration -----------------
# set -xv
cp="cp"
mv="mv"

# Control if BFMDIR and NEMODIR are defined among environment variables
if [ "${BFMDIR}" = "" ]; then
   echo "Environmental variable BFMDIR not defined!"
   exit 0
fi
if [ "${NEMODIR}" = "" ]; then
   echo "Environmental variable NEMODIR not defined!"
   exit 0
fi
# Control if TCLSH is defined among environment variables
if [ "${TCLSH}" = "" ]; then
   echo "Environmental variable TCLSH not defined!"
   echo "You cannot generate the BFM code."
   echo "tclsh 8.4 is REQUIRED."
   exit 0
fi

if [ ! -d ${BFMDIR}/build/${BLDDIR} ]; then
  mkdir ${BFMDIR}/build/${BLDDIR}
fi

cd ${BLDDIR}

# generate BFM Memory Layout files
${BFMDIR}/build/scripts/GenerateGlobalBFMF90Code ${cppdefs} \
          -read ${BFMDIR}/build/Configurations/${myGlobalDef} \
          -from ${BFMDIR}/src/BFM/proto \
          -to ${BFMDIR}/build/${BLDDIR} -actions statemem allocmem netcdfmem \
          -to ${BFMDIR}/build/${BLDDIR} -actions headermem \
          -to ${BFMDIR}/build/${BLDDIR} -actions initmem

# Generate the specific bfm.fcm include file for makenemo
cppdefs=`echo ${cppdefs} | sed -e "s/"-D"//g"` 
# some macros are default with NEMO
FCMMacros="BFM_NEMO USEPACK BFM_NOPOINTERS ${cppdefs}"
sed -e "s/_place_keys_/${FCMMacros}/" -e "s/_place_def_/${myGlobalDef}/" \
       ${BFMDIR}/build/Configurations/Default_bfm.fcm > ${BFMDIR}/build/${BLDDIR}/bfm.fcm

echo "Memory Layout generated in local folder: ${BLDDIR}."

# If COMPILE, copy files to final destination and launch makenemo
if [ ${COMPILE} = "yes" ]; then

  # Move BFM Layout files to target folders 
  ${cp} ${BFMDIR}/build/${BLDDIR}/*.F90 ${BFMDIR}/src/BFM/General
  ${mv} ${BFMDIR}/src/BFM/General/init_var_bfm.F90 ${BFMDIR}/src/share
  ${cp} ${BFMDIR}/build/${BLDDIR}/INCLUDE.h ${BFMDIR}/src/BFM/include
  ${cp} ${BFMDIR}/build/${BLDDIR}/bfm.fcm ${BFMDIR}/src/nemo
  echo "Files copied to target folders."

  back=${PWD}
  echo "Starting NEMO compilation..."
  cd ${NEMODIR}/NEMOGCM/CONFIG
  if [ ${MODE} = "clean" ]; then
    echo "Cleaning up ${CONFIG}..."
    ./makenemo -n ${CONFIG} -m ${ARCH} clean
  fi
  ./makenemo -n ${CONFIG} -m ${ARCH} -e ${BFMDIR}/src/nemo -j ${NPROC}
  cd ${back}
  echo "${CONFIG} compilation done!"
fi
