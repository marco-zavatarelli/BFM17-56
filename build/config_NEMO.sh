#!/bin/sh
#
## Configuration for the BFM-NEMO coupling. 
#
#  This script creates a directory $BLDDIR with the Memory Layout files and the FCM include for the coupling. 
#  The RELEASE flag is used to copy the configuration files to the final destination folders (before compilation).
# 

#  Currently available macros (cppdefs) are:
#  INCLUDE_PELFE                          : use Iron component to the pelagic system
#  INCLUDE_PELCO2, INCLUDE_BENCO2         : activate Carbonate System 
#  INCLUDE_BEN, INCLUDE_BENPROFILES       : Add Benthic compartment
#  INCLUDE_SILT                           : use Silt component
#  INCLUDE_SEAICE                         : activate SeaIce Ecology 
#  INCLUDE_DIAG3D, INCLUDE_DIAG2D         : additional diagnostics available for output
#  USEPACK                                : use pack/unpack fortran intrinsic function
#  BFM_PARALLEL                           : used to generate bfm Log file in parallel job 

#  Warnings
# 1. Still not working for benthic BFM don't use DIAG with D1SOURCE and ONESOURCE
# 2. Adding the BFM_NOPOINTERS key to compile with gfortran 4.5 and older
# 3. Using the key DEBUG will add more output information

#----------------- User configuration -----------------
# myGlobalDef   : file used by tcl script (GenerateGlobalBFMF90Code) to design the BFM Memory Layout
# cppdefs       : keys used to configure the model
# RELEASE       : set "yes" to copy the Memory Layout files to final destination folders
# BLDDIR        : Local folder containing the generated BFM Memory Layout files
# -----------------------------------------------------
myGlobalDef="GlobalDefsBFM.model.standard"
cppdefs=" -DBFM_PARALLEL -DONESOURCE -DINCLUDE_PELCO2" 

RELEASE="yes"

BLDDIR="STD_MemLayout"
#----------------- User configuration -----------------
# set -xv
cp="cp"
mv="mv"

# Control if BFMDIR is defined among environment variables
if [ "${BFMDIR}" = "" ]; then
   echo "Environmental variable BFMDIR not defined!"
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
FCMMacros="BFM_NEMO USEPACK BFM_NOPOINTERS ${cppdefs}"
sed -e "s/_place_keys_/${FCMMacros}/" -e "s/_place_def_/${myGlobalDef}/" \
       ${BFMDIR}/build/Configurations/Default_bfm.fcm > ${BFMDIR}/build/${BLDDIR}/bfm.fcm

echo "Memory Layout generated in local folder: ${BLDDIR}."

# If a RELEASE Configuration copy files to final destination
if [ ${RELEASE} = "yes" ]; then

  # Move BFM Layout files to target folders
  ${cp} ${BFMDIR}/build/${BLDDIR}/*.F90 ${BFMDIR}/src/BFM/General
  ${mv} ${BFMDIR}/src/BFM/General/init_var_bfm.F90 ${BFMDIR}/src/share
  ${cp} ${BFMDIR}/build/${BLDDIR}/INCLUDE.h ${BFMDIR}/src/BFM/include
  ${cp} ${BFMDIR}/build/${BLDDIR}/bfm.fcm ${BFMDIR}/src/nemo

echo "Files copied to target folders."
 
fi

echo "DONEEEE !!!!"
