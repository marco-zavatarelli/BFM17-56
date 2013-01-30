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

# Author: Esteban Gutierrez (CMCC) based on Tomas lobato (CMCC) scripts to config BFM
# -----------------------------------------------------


MODE="STANDALONE"
PRESET="STANDALONE"
ARCH_NEMO="PW6_calypso"
ARCH_STD="gfortran"
NPROC=4
CLEAN="clean"
CPPDEFS="-DINCLUDE_PELCO2 -DINCLUDE_DIAG3D"

EXP="EXP00"
PROC=4
QUEUE="poe_short"

MKMF="mkmf"
GMAKE="gmake"
PERL="perl"
GENCONF="generate_conf"
MKNEMO="makenemo"
BFMSTD="bfm_standalone.x"
NEMOEXE="opa"

usage(){
    cat << EOF
usage: $0 -h
usage: $0 {-g -c -e} [options]

This script compile and/or execute the BFM model.

MUST specify at least one these OPTIONS:
   -h      shows this help
   -g      generate
   -c      compile with makenemo (include generation)
   -e      execute

alternative COMPILATION OPTIONS are:
   -m MODE
              Use with or without NEMO. Valid values are STANDALONE or NEMO. (Default: "STANDALONE")
   -p PRESET
              Preset to generate the configuration. (Default: "STANDALONE")
              - For other presets, list files *.conf in: BFMDIR/build/Configurations
   -v
              Verbose mode to print all messages (Deactivated by default)
   -b BFMDIR
              the environmental variable BFMDIR pointing to the root directory of BFM (Default: "${BFMDIR}")
   -n NEMODIR
              the environmental variable NEMODIR pointing to the root directory of NEMO. (Default: "${NEMODIR}")
   -a ARCH
              NEMO specific architecture file (Default: "gfortran" for STANDALONE or "PW6_calypso" for NEMO)
              - For STANDALONE available archs, list dir : BFMDIR/compilers
              - For NEMO available archs, execute command: NEMODIR/NEMOGCM/CONFIG/makenemo -h all
   -r PROC
              number of procs used for compilation. Default: 4
alternative EXECUTION OPTIONS are:
   -x EXP
              Name of the experiment for generation of the output (Default: "EXP00")
   -l NMLDIR
              input dir where are the namelists to run the experiment (Default: "BFMDIR/build/PRESET")
   -r PROC
              number of procs used for running. Default: 4
   -q QUEUE
              name of the queue number of procs used for running. Default
EOF
}



while getopts "hgcem:p:vb:n:a:r:x:l:q:" opt; do
    case $opt in
      h )                   usage; exit                                     ;;
      g ) [ ${VERBOSE} ] && echo "generation activated"  ; GEN=1            ;;
      c ) [ ${VERBOSE} ] && echo "compilation activated" ; CMP=1; GEN=1     ;;
      e ) [ ${VERBOSE} ] && echo "execution activated"   ; EXE=1            ;;
      m ) [ ${VERBOSE} ] && echo "mode $OPTARG"          ; MODE=$OPTARG     ;;
      p ) [ ${VERBOSE} ] && echo "preset $OPTARG"        ; PRESET=$OPTARG   ;;
      v )                   echo "verbose mode"          ; VERBOSE=1        ;;
      b ) [ ${VERBOSE} ] && echo "BFMDIR=$OPTARG"        ; BFMDIR=$OPTARG   ;;
      n ) [ ${VERBOSE} ] && echo "NEMODIR=$OPTARG"       ; NEMODIR=$OPTARG  ;;
      a ) [ ${VERBOSE} ] && echo "architecture $OPTARG"  ; ARCH_OPT=$OPTARG ;;
      r ) [ ${VERBOSE} ] && echo "n. procs $OPTARG"      ; PROC=$OPTARG     ;;
      x ) [ ${VERBOSE} ] && echo "experiment $OPTARG"    ; EXP=$OPTARG      ;;
      l ) [ ${VERBOSE} ] && echo "namelist dir $OPTARG"  ; NMLDIR=$OPTARG   ;;
      q ) [ ${VERBOSE} ] && echo "queue name $OPTARG"    ; QUEUE=$OPTARG    ;;
    esac
done

if [[ ! ${EXE} && ! ${CMP} && ! ${GEN} ]]; then
    echo "YOU MUST specify one of the \"must\" arguments";
    usage;
    exit;
fi
if [[ ${MODE} != "STANDALONE" && ${MODE} != "NEMO" ]]; then
    echo "MODE is not correct\n";
    usage;
    exit;
fi
if [[ ! $BFMDIR || ! $NEMODIR ]]; then 
    echo "BFMDIR and/or NEMODIR not specified"; 
    exit; 
fi
if [[ ${PROC} ]] && ! [[ "$PROC" =~ ^[0-9]+$ ]] ; then 
    echo "PROC must be a number"; 
    exit; 
fi


if [ $VERBOSE ]; then
    cmd_mkmf="${MKMF} -v"
    cmd_gmake="${GMAKE}"
    cmd_gen="${GENCONF}.pl -v"
    cmd_mknemo="${MKNEMO}"
else
    cmd_mkmf="${MKMF}"
    cmd_gmake="${GMAKE} -s"
    cmd_gen="${GENCONF}.pl"
    cmd_mknemo="${MKNEMO} -v0"
fi

blddir="${BFMDIR}/build/${PRESET}"
if [ ${GEN} ]; then
    myGlobalDef="${PRESET}.conf"

    if [ ! -f ${BFMDIR}/build/Configurations/${myGlobalDef} ]; then
         echo "ERROR: ${BFMDIR}/build/Configurations/${myGlobalDef} not exsits"
         exit
    fi

    if [ ! -d ${blddir} ]; then mkdir ${blddir}; fi
    cd ${blddir}
    rm -rf *
    
    if [[ ${MODE} == "STANDALONE" ]]; then
        if [[ ${ARCH_OPT} ]]; then ARCH=${ARCH_OPT}; else ARCH=${ARCH_STD}; fi
    else
        if [[ ${ARCH_OPT} ]]; then ARCH=${ARCH_OPT}; else ARCH=${ARCH_NEMO}; fi
    fi

    # generate BFM Memory Layout files
    ${PERL} -I${BFMDIR}/build/scripts/conf/ ${BFMDIR}/build/scripts/conf/${cmd_gen} \
        ${CPPDEFS} \
        -r ${BFMDIR}/build/Configurations/${myGlobalDef}  \
        -f ${BFMDIR}/src/BFM/proto \
        -t ${blddir} || exit

    if [[ ${MODE} == "STANDALONE" ]]; then
        cppdefs="-DBFM_STANDALONE ${CPPDEFS}"
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

        #change netcdf path in file if Mac
        if [ `uname -s` == "Darwin" ]; then
            [ ${VERBOSE} ] && echo "changing netcd path for Mac!"
            sed -e "s/\/usr\/local/\/opt\/local/" ${BFMDIR}/compilers/${ARCH}.inc > ${blddir}/${ARCH}.inc
        fi

        # Make makefile
        ${BFMDIR}/bin/${cmd_mkmf} \
            -c "${cppdefs}" \
            -o "-I${BFMDIR}/include -I${BFMDIR}/src/BFM/include" \
            -t "${blddir}/${ARCH}.inc" \
            -p "${BFMDIR}/bin/bfm_standalone.x" \
            BFM.lst && echo ""

        # Link to the configuration file
        ln -sf ${BFMDIR}/build/Configurations/${myGlobalDef} GlobalDefsBFM.model
        [ ${VERBOSE} ] && echo "${PRESET} compilation done!"

        # If COMPILE, launch gmake
        if [ ${CMP} ]; then
            if [ ${CLEAN} == "clean" ]; then
                [ ${VERBOSE} ] && echo "Cleaning up ${PRESET}..."
                ${cmd_gmake} clean
            fi
            [ ${VERBOSE} ] && echo "Starting ${PRESET} compilation..."
            ${cmd_gmake}
            [ ${VERBOSE} ] && echo "${PRESET} compilation done!"
        fi
    else
        cppdefs="-DBFM_PARALLEL ${CPPDEFS}"
        # Generate the specific bfm.fcm include file for makenemo
        cppdefs=`echo ${cppdefs} | sed -e "s/"-D"//g"` 
        # some macros are default with NEMO
        FCMMacros="BFM_NEMO USEPACK BFM_NOPOINTERS ${cppdefs}"
        sed -e "s/_place_keys_/${FCMMacros}/" -e "s/_place_def_/${myGlobalDef}/" \
            ${BFMDIR}/build/Configurations/Default_bfm.fcm > ${blddir}/bfm.fcm
        [ ${VERBOSE} ] && echo "Memory Layout generated in local folder: ${blddir}."

        # Move BFM Layout files to target folders 
        cp ${blddir}/*.F90 ${BFMDIR}/src/BFM/General
        cp ${blddir}/init_var_bfm.F90 ${BFMDIR}/src/share
        cp ${blddir}/INCLUDE.h ${BFMDIR}/src/BFM/include
        cp ${blddir}/bfm.fcm ${BFMDIR}/src/nemo
        [ ${VERBOSE} ] && echo "Files copied to target folders."

        # If COMPILE, launch makenemo
        if [ ${CMP} ]; then
            cd ${NEMODIR}/NEMOGCM/CONFIG/

            if [ ${CLEAN} == "clean" ]; then
                [ ${VERBOSE} ] && echo "Cleaning up ${PRESET}..."
                ./${cmd_mknemo} -n ${PRESET} -m ${ARCH} clean
            fi
            [ ${VERBOSE} ] && echo "Starting NEMO compilation..."
            ./${cmd_mknemo} -n ${PRESET} -m ${ARCH} -e ${BFMDIR}/src/nemo -j ${NPROC}
            [ ${VERBOSE} ] && echo "${PRESET} compilation done!"
        fi
    fi
fi

if [ ${EXE} ]; then
    [ ${VERBOSE} ] && echo "Executing ${PRESET}"

    if [ ! -d ${blddir} ]; then
        echo "ERROR: directory ${blddir} not exists"
    fi
 
    exedir="${BFMDIR}/run/${EXP}"
    if [ ! -d ${exedir} ]; then 
        mkdir ${exedir}; 
    fi
    cd ${exedir}
    rm -rf *

    if [[ ${MODE} == "STANDALONE" ]]; then
        # copy and link necessary files
        if [ ${NMLDIR} ]; then cp ${NMLDIR}/* .; else cp ${blddir}/*.nml .; fi
        ln -sf ${BFMDIR}/bin/${BFMSTD} ${BFMSTD}

        #execute
        ./${BFMSTD}
    else
        #change values in runscript
        sed -e "s,_EXP_,${EXP},g"       \
            -e "s,_PRESET_,${PRESET},g" \
            -e "s,_QUEUE_,${QUEUE},g"   \
            -e "s,_PROC_,${PROC},g"     ${BFMDIR}/build/scripts/conf/runscript > ./runscript_${EXP}

        # copy and link necessary files
        cp ${BFMDIR}/build/scripts/conf/xmlio_server.def ./
        if [ ${NMLDIR} ]; then cp ${NMLDIR}/* .; else cp ${blddir}/*.nml .; fi
        ln -sf ${NEMODIR}/NEMOGCM/CONFIG/${PRESET}/BLD/bin/${NEMOEXE} ${NEMOEXE}.x

        #execute
        bsub < ./runscript_${EXP}
        [ ${VERBOSE} ] && echo "Execution logs will be generated in ${exedir}"
    fi
    [ ${VERBOSE} ] && echo "Output generated in ${exedir}"
fi
