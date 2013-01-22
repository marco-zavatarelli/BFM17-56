#!/bin/bash
#
## Configuration for the BFM-NEMO coupling in the PELAGOS configuration. 
#
#  This script creates a directory $blddir with the Memory Layout files and the FCM include for the coupling. 
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

# Author: Esteban Gutierrez (CMCC) esteban.gutierrez@cmcc.it
# -----------------------------------------------------


PRESET="STANDALONE"
ARCH_OTH="PW6_calypso"
ARCH_STD="gfortran"
NPROC=4
MODE="clean"
MODEL="standard"
GEN=0;
EXE=0;
CMP=0;
VERBOSE=0;

usage(){
    cat << EOF
usage: $0 -h
usage: $0 {-g -c -e} [options]

This script compile and/or execute the BFM model.

must specify at least one these options:
   -h      shows this help
   -g      generate
   -c      compile with makenemo (include generation)
   -e      execute

alternative options are:
   -p PRESET
              Preset to generate the configuration. Available presets are: STANDALONE, GYRE_BFM, NEMO and PELAGOS (Default: "STANDALONE")
   -v
              Verbose mode to print all messages (Deactivated by default)
   -b BFMDIR
              the environmental variable BFMDIR pointing to the root directory of BFM (Default: "${BFMDIR}")
   -n NEMODIR
              the environmental variable NEMODIR pointing to the root directory of NEMO. (Default: "${NEMODIR}")
   -m MODEL
              option used to design the BFM Memory Layout. Available options are: standard, iron, seaice (Default: "standard")
   -a ARCH
              NEMO specific architecture file (Default: "gfortran" for STANDALONE or "PW6_calypso" for others)
              - For STANDALONE, list dir : ${BFMDIR}/compilers
              - For other presets, execute command: ${NEMODIR}/NEMOGCM/CONFIG/makenemo -h all
   -r PROC
              number of procs used for compilation. Default: 4
EOF
}



while getopts "hecp:m:vb:n:a:r:" opt; do
    case $opt in
      h ) usage; exit                                     ;;
      g ) echo "generation activated"  ; GEN=1            ;;
      e ) echo "execution activated"   ; EXE=1            ;;
      c ) echo "compilation activated" ; CMP=1; GEN=1     ;;
      p ) echo "preset $OPTARG"        ; PRESET=$OPTARG   ;;
      m ) echo "model $OPTARG"         ; MODEL=$OPTARG    ;;
      v ) echo "verbose mode"          ; VERBOSE=1        ;;
      b ) echo "BFMDIR=$OPTARG"        ; BFMDIR=$OPTARG   ;;
      n ) echo "NEMODIR=$OPTARG"       ; NEMODIR=$OPTARG  ;;
      a ) echo "architecture $OPTARG"  ; ARCH_OPT=$OPTARG ;;
      r ) echo "n. procs $OPTARG"      ; PROC=$OPTARG     ;;
    esac
done

if [[ ${EXE} == 0 && ${CMP} == 0 && ${GEN} == 0 ]] ; then
    usage
    exit
fi
if [[ ${PRESET} && ${PRESET} != "PELAGOS" && ${PRESET} != "GYRE_BFM" && ${PRESET} != "STANDALONE" && ${PRESET} != "NEMO" ]]; then
    usage
    exit
fi
if [[ ${MODEL} && ${MODEL} != "standard" && ${MODEL} != "iron" && ${MODEL} != "ice" ]]; then
    usage
    exit
fi
if [[ ! $BFMDIR || ! $NEMODIR ]]; then
    echo "BFMDIR and/or NEMODIR not specified"
    exit
fi
if [[ ${PROC} ]] && ! [[ "$PROC" =~ ^[0-9]+$ ]] ; then
    echo "PROC must be a number"
    exit
fi





if [[ ${GEN} -eq 1 ]]; then
    myGlobalDef="GlobalDefsBFM.model.${MODEL}"
    blddir="${BFMDIR}/build/STD_${PRESET}"

    if [[ ${PRESET} == "STANDALONE"   ]]; then
        cppdefs="-DBFM_STANDALONE -DINCLUDE_PELCO2 -DINCLUDE_DIAG3D"
        if [[ ${ARCH_OPT} ]]; then ARCH=${ARCH_OPT}; else ARCH=${ARCH_STD}; fi
    else
        cppdefs="-DBFM_PARALLEL -DINCLUDE_PELCO2 -DINCLUDE_DIAG3D"
        if [[ ${ARCH_OPT} ]]; then ARCH=${ARCH_OPT}; else ARCH=${ARCH_OTH}; fi
    fi

    if [ ! -d ${blddir} ]; then
        mkdir ${blddir}
    fi

    cd ${blddir}
    rm -rf *

    # generate BFM Memory Layout files
    perl -I${BFMDIR}/build/scripts/conf/ ${BFMDIR}/build/scripts/conf/read_conf.pl ${cppdefs} \
        -r ${BFMDIR}/build/Configurations/${myGlobalDef}  \
        -n "asd"   \
        -f ${BFMDIR}/src/BFM/proto \
        -t ${blddir} || exit

    if [[ ${PRESET} == "STANDALONE" ]]; then
        # list files
        find ${BFMDIR}/src/BFM/General -name "*.?90" -print > BFM.lst
        find ${BFMDIR}/src/standalone -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/share -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Pel -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/PelBen -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Ben -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Light -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Oxygen -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Forcing -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/CO2 -name "*.?90" -print >> BFM.lst

        #change netcdf path in file
        sed -e "s/\/usr\/local/\/opt\/local/" ${BFMDIR}/compilers/${ARCH}.inc > ${blddir}/${ARCH}.inc

        # Make makefile
        "${BFMDIR}/bin/mkmf" \
            -c "${cppdefs}" \
            -o "-I${BFMDIR}/include -I${BFMDIR}/src/BFM/include" \
            -t "${blddir}/${ARCH}.inc" \
            -p "${BFMDIR}/bin/bfm_standalone.x" \
            BFM.lst

        # Link to the configuration file
        ln -sf ${BFMDIR}/build/Configurations/${myGlobalDef} GlobalDefsBFM.model
        echo "${PRESET} compilation done!"

        if [ ${CMP} == 1 ]; then
            if [ ${MODE} == "clean" ]; then
                gmake clean
            fi
            gmake
        fi
    else

        # Generate the specific bfm.fcm include file for makenemo
        cppdefs=`echo ${cppdefs} | sed -e "s/"-D"//g"` 
        # some macros are default with NEMO
        FCMMacros="BFM_NEMO USEPACK BFM_NOPOINTERS ${cppdefs}"
        sed -e "s/_place_keys_/${FCMMacros}/" -e "s/_place_def_/${myGlobalDef}/" \
            ${BFMDIR}/build/Configurations/Default_bfm.fcm > ${blddir}/bfm.fcm

        echo "Memory Layout generated in local folder: ${blddir}."

        # If COMPILE, copy files to final destination and launch makenemo
        if [ ${CMP} == 1 ]; then

            # Move BFM Layout files to target folders 
            cp ${blddir}/*.F90 ${BFMDIR}/src/BFM/General
            mv ${BFMDIR}/src/BFM/General/init_var_bfm.F90 ${BFMDIR}/src/share
            cp ${blddir}/INCLUDE.h ${BFMDIR}/src/BFM/include
            cp ${blddir}/bfm.fcm ${BFMDIR}/src/nemo
            echo "Files copied to target folders."

            cd ${NEMODIR}/NEMOGCM/CONFIG
            if [ ${MODE} == "clean" ]; then
                echo "Cleaning up ${PRESET}..."
                ./makenemo -n ${PRESET} -m ${ARCH} clean
            fi
            echo "Starting NEMO compilation..."
            ./makenemo -n ${PRESET} -m ${ARCH} -e ${BFMDIR}/src/nemo -j ${NPROC}
            echo "${PRESET} compilation done!"
        fi
    fi
fi
