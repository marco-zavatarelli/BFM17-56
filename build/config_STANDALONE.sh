#!/bin/sh
#
## Configuration file for Pelagic BFM STANDALONE.
#
#  This script creates a directory BLD_STANDALONE with the appropriate makefile. 
#  The STANDALONE configuration is created in the default folder BLD_STANDALONE. 
#  Note: use the gmake to compile the Makefile. 

#  Currently available macros (cppdefs) are:
#  INCLUDE_PELFE                          : use Iron component to the pelagic system
#  INCLUDE_PELCO2, INCLUDE_BENCO2         : activate Carbonate System 
#  INCLUDE_BEN, INCLUDE_BENPROFILES       : Add Benthic compartment
#  INCLUDE_SILT                           : use Silt component
#  INCLUDE_SEAICE                         : activate SeaIce Ecology
#  INCLUDE_DIAG3D, INCLUDE_DIAG2D         : additional diagnostics available for output
#  USEPACK                                : use pack/unpack fortran intrinsic function

#  Warnings
# 1. Still not working for benthic BFM don't use DIAG with D1SOURCE and ONESOURCE
# 2. Adding the BFM_NOPOINTERS key to compile with gfortran 4.5 and older
# 3. Using the key DEBUG will add more output information

#----------------- User configuration -----------------
# myGlobalDef   : file used by tcl script (GenerateGlobalBFMF90Code) to design the BFM Memory Layout
# cppdefs       : keys used to configure the model
# archfile      : compiler architecture file
# exe           : default location of the BFM standalone executable
# -----------------------------------------------------
myGlobalDef="GlobalDefsBFM.model.standard"
cppdefs="-DBFM_STANDALONE"

archfile="${BFMDIR}/compilers/gfortran.inc"
exe=${BFMDIR}/bin/bfm_standalone.x
#----------------- User configuration -----------------
#set -xv

# Control if BFMDIR is defined among environment variables
if [ "${BFMDIR}" = "" ]; then
   echo "Environmental variable BFMDIR not defined!"
   exit 0
fi

# set makefile options and destination
BLDDIR="./BLD_STANDALONE"
MKMF="${BFMDIR}/bin/mkmf"
oflags="-I${BFMDIR}/include -I${BFMDIR}/src/BFM/include"

if [ ! -d ${BLDDIR} ]; then
  mkdir ${BLDDIR}
fi

cd ${BLDDIR}
# Link to the configuration file
 ln -sf ${BFMDIR}/build/Configurations/${myGlobalDef} ${BFMDIR}/src/BFM/General/GlobalDefsBFM.model

# generate BFM files
${BFMDIR}/build/scripts/GenerateGlobalBFMF90Code  ${cppdefs} \
          -read ${BFMDIR}/build/Configurations/${myGlobalDef} \
          -from ${BFMDIR}/src/BFM/proto \
          -to ${BFMDIR}/src/BFM/General -actions statemem allocmem netcdfmem \
          -to ${BFMDIR}/src/BFM/include -actions headermem \
          -to ${BFMDIR}/src/share -actions initmem

# list files
find ${BFMDIR}/src/BFM/General -name "*.?90" -print > BFM.lst
find ${BFMDIR}/src/standalone -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/share -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/PelB -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/PelBen -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/Ben -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/Light -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/Oxygen -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/Forcing -name "*.?90" -print >> BFM.lst
find ${BFMDIR}/src/BFM/CO2 -name "*.?90" -print >> BFM.lst

# Make makefile
${MKMF} -c "${cppdefs}" -o "${oflags}" -t ${archfile} -p ${exe} BFM.lst
