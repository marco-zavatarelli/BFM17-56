# DESCRIPTION
#   BFM Configuration manager
#
# AUTHORS
#   Esteban Gutierrez esteban.gutierrez@cmcc.it
#   Tomas Lovato toma.lovato@cmcc.it
#
# COPYING
#  
#   Copyright (C) 2013 BFM System Team ( bfm_st@lists.cmcc.it )
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation;
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# -----------------------------------------------------

#!/bin/bash -e

#logging configuration
LOGFILE=logfile_$$.log
LOGDIR="Logs"

#local paths
TEMPDIR="build/tmp"
CONFDIR="build/configurations"
SCRIPTSDIR="build/scripts"
SCRIPTS_CONFDIR="${SCRIPTSDIR}/conf"

#programs
MKMF="mkmf"
GMAKE="gmake"
PERL="perl"
GENCONF="generate_conf"
MKNEMO="makenemo"
BFMSTD="bfm_standalone.x"
NEMOEXE="nemo.exe"

#nemo files
NEMO_FILES="iodef.xml namelist namelist_top xmlio_server.def";

#options
OPTS="hvgcdPp:m:k:b:n:a:r:ft:x:l:q:"
OPTIONS=(     MODE     CPPDEFS     BFMDIR     NEMODIR     ARCH     CLEAN     PROC     NETCDF     EXP     NMLDIR     PROC     QUEUE     )
OPTIONS_USR=( mode_usr cppdefs_usr bfmdir_usr nemodir_usr arch_usr clean_usr proc_usr netcdf_usr exp_usr nmldir_usr proc_usr queue_usr )

#error message
ERROR_MSG="Execute $0 -h for help if you don't know what the hell is going wrong. PLEASE read CAREFULLY before bother someone else"

#----------------- USER CONFIGURATION DEFAULT VALUES -----------------
MODE="STANDALONE"
CPPDEFS="BFM_STANDALONE INCLUDE_PELCO2 INCLUDE_DIAG3D"
PRESET="STANDALONE_PELAGIC"
ARCH="gfortran.inc"
PROC=4
EXP="EXP00"
QUEUE="poe_short"
CLEAN=1
# --------------------------------------------------------------------


#print usage message 
usage(){
    more << EOF
NAME
    This script compile and/or deploy the BFM model.

SYNOPSIS
    usage: $0 -h
    usage: $0 -P
    usage: $0 {-g -c -d} [options]

DESCRIPTION
    MUST specify at least one these OPTIONS:
       -h         Shows this help
       -P         List all available presets
       -g         Generate ".H" ".F90" and ".NML" files
       -c         Compile
       -d         Deploy the execution environment (experiment folder creation)

    alternative COMPILATION OPTIONS are:
       -v
                  Verbose mode to print all messages (Deactivated by default)
       -p PRESET
                  Preset to generate the configuration. (Default: "${PRESET}")
                  - For other presets, execute $0 -P 
       -m MODE
                  Mode for compilation and deployment. Available models are: (Default: "STANDALONE")
                  - STANDALONE (without NEMO. Compile and run in local machine)
                  - NEMO (with NEMO. Compile and run ONLY in LSF platform)
       -k CPPDEFS
                  Key options to configure the model. (Default: "BFM_STANDALONE INCLUDE_PELCO2 INCLUDE_DIAG3D")                 
       -b BFMDIR
                  The environmental variable BFMDIR pointing to the root directory of BFM (Default: "${BFMDIR}")
       -n NEMODIR
                  The environmental variable NEMODIR pointing to the root directory of NEMO. (Default: "${NEMODIR}")
       -a ARCH
                  Specify compilation Architecture file (Default: "gfortran.inc")
                  - For STANDALONE mode available archs, list dir : BFMDIR/compilers
                  - For NEMO mode available archs, execute command: NEMODIR/NEMOGCM/CONFIG/makenemo -h all
       -r PROC
                  Number of procs used for compilation. Default: 4
       -f
                  Fast mode. Dont execute "clean" command in compilation (clean is activated by default)
       -t NETCDF
                  Path to netcdf library and header files. (Default: /usr/local)
    alternative DEPLOYMENT OPTIONS are:
       -x EXP
                  Name of the experiment for generation of the output folder (Default: "EXP00")
       -l NMLDIR
                  Input dir where are the namelists to run the experiment (Default: "${TEMPDIR}/${PRESET}")
       -r PROC
                  Number of procs used for running (same as compilation). Default: 4
       -q QUEUE
                  Name of the queue number of procs used for running. Default
EOF
}



# ------------------------- PROGRAM STARTS HERE ----------------------------

#print in log file
if [ ! -d ${LOGDIR} ]; then mkdir ${LOGDIR}; fi
mkfifo ${LOGDIR}/${LOGFILE}.pipe
tee < ${LOGDIR}/${LOGFILE}.pipe ${LOGDIR}/${LOGFILE} &
exec &> ${LOGDIR}/${LOGFILE}.pipe
rm ${LOGDIR}/${LOGFILE}.pipe


#get user options from commandline
while getopts "${OPTS}" opt; do
    case $opt in
      h ) usage;            rm ${LOGDIR}/${LOGFILE}      ; exit             ;;
      v )                   echo "verbose mode"          ; VERBOSE=1        ;;
      g ) [ ${VERBOSE} ] && echo "generation activated"  ; GEN=1            ;;
      c ) [ ${VERBOSE} ] && echo "compilation activated" ; CMP=1            ;;
      d ) [ ${VERBOSE} ] && echo "deployment activated"  ; DEP=1            ;;
      P ) [ ${VERBOSE} ] && echo "list presets"          ; LIST=1           ;;
      p ) [ ${VERBOSE} ] && echo "preset $OPTARG"        ; PRESET=$OPTARG      ;;
      m ) [ ${VERBOSE} ] && echo "mode $OPTARG"          ; mode_usr=$OPTARG    ;;
      k ) [ ${VERBOSE} ] && echo "key options $OPTARG"   ; cppdefs_usr=$OPTARG ;;
      b ) [ ${VERBOSE} ] && echo "BFMDIR=$OPTARG"        ; bfmdir_usr=$OPTARG  ;;
      n ) [ ${VERBOSE} ] && echo "NEMODIR=$OPTARG"       ; nemodir_usr=$OPTARG ;;
      a ) [ ${VERBOSE} ] && echo "architecture $OPTARG"  ; arch_usr=$OPTARG    ;;
      f ) [ ${VERBOSE} ] && echo "fast mode activated"   ; clean_usr=0         ;;
      t ) [ ${VERBOSE} ] && echo "netcdf path $OPTARG"   ; netcdf_usr=$OPTARG  ;;
      x ) [ ${VERBOSE} ] && echo "experiment $OPTARG"    ; exp_usr=$OPTARG     ;;
      l ) [ ${VERBOSE} ] && echo "namelist dir $OPTARG"  ; nmldir_usr=$OPTARG  ;;
      r ) [ ${VERBOSE} ] && echo "n. procs $OPTARG"      ; proc_usr=$OPTARG    ;;
      q ) [ ${VERBOSE} ] && echo "queue name $OPTARG"    ; queue_usr=$OPTARG   ;;
      * ) echo "option not recognized"                   ; exit             ;;
    esac
done

#check must parameters
if [[ $LIST ]]; then
    for pre in `ls ${BFMDIR}/${CONFDIR}`; do printf " - $pre\n"; done
    exit
fi
if [[ ! ${DEP} && ! ${CMP} && ! ${GEN} ]]; then
    echo "ERROR: YOU MUST specify one of the \"must\" arguments"
    echo ${ERROR_MSG}
    exit
fi

#activate/deactivate verbose mode
if [ $VERBOSE ]; then
    #set -xv
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

myGlobalConf="${PRESET}/configuration"
myGlobalMem="${PRESET}/layout"
myGlobalNml="${PRESET}/namelists"

##### Overwrite options specified in configuration file
[ ${VERBOSE} ] && echo "Conf. file replace:"
bfmconf=`perl -ne "/BFM_conf/ .. /^\// and print" ../${CONFDIR}/${myGlobalConf}`
for option in "${OPTIONS[@]}"; do
    value=`perl -e "print ( \"${bfmconf}\" =~ m/\${option}\ *=\ *[\"\']*([^\"\'\,]+)[\"\']*[\,\/]*/ );"`
    if [ "${value}" ]; then 
        [ ${VERBOSE} ] && echo "- replacing ${option}=${value}"
        eval ${option}=\"\${value}\"
    fi
done

##### Overwrite options specified in command line by user
[ ${VERBOSE} ] && echo "Command line replace:"
for option in "${OPTIONS_USR[@]}"; do
    opt_name=`echo $option | sed -e 's/_usr//' | awk '{print toupper($0)}'`
    eval [ \"\${$option}\" ] && eval "${opt_name}"=\"\${$option}\"
    [ ${VERBOSE} ] && eval [ \"\${$option}\" ]  && eval echo "- replacing ${opt_name}="\"\${$option}\"
done

#specify build dir of BFM
blddir="${BFMDIR}/${TEMPDIR}/${PRESET}"

#Check some optional parameter values
if [[ ! $BFMDIR ]]; then 
    echo "ERROR: BFMDIR not specified"
    echo ${ERROR_MSG}
    exit
fi
if [[ ! $NEMODIR && "$MODE" == "NEMO" ]]; then
    echo "ERROR: NEMODIR not specified in NEMO mode"
    echo ${ERROR_MSG}
    exit
fi
if [[ "$MODE" != "STANDALONE" && "$MODE" != "NEMO" ]]; then 
    echo "ERROR: MODE value not valid ($MODE). Available values are: STANDALONE or NEMO."
    echo ${ERROR_MSG}
    exit
fi
if [[ ${PROC} ]] && ! [[ "$PROC" =~ ^[0-9]+$ ]] ; then 
    echo "ERROR: PROC must be a number"
    echo ${ERROR_MSG}
    exit
fi




# -----------------------------------------------------
# Memory and namelist files GENERATION
# -----------------------------------------------------


if [ ${GEN} ]; then

    if [ ! -f ${BFMDIR}/${CONFDIR}/${myGlobalMem} ]; then
         echo "ERROR: ${BFMDIR}/${CONFDIR}/${myGlobalMem} not exsits"
         echo ${ERROR_MSG}
         exit
    fi
    if [ ! -f ${BFMDIR}/${CONFDIR}/${myGlobalNml} ]; then
         echo "ERROR: ${BFMDIR}/${CONFDIR}/${myGlobalNml} not exsits"
         echo ${ERROR_MSG}
         exit
    fi


    if [ ! -d ${blddir} ]; then mkdir ${blddir}; fi
    cd ${blddir}
    rm -rf *
    
    #add -D to cppdefs
    cppdefs=`echo ${CPPDEFS} | sed -e 's/\([a-zA-Z_0-9]*\)/-D\1/g'` 

    # generate BFM Memory Layout files and namelists
    ${PERL} -I${BFMDIR}/${SCRIPTS_CONFDIR}/ ${BFMDIR}/${SCRIPTS_CONFDIR}/${cmd_gen} \
        ${cppdefs} \
        -r ${BFMDIR}/${CONFDIR}/${myGlobalMem}  \
        -n ${BFMDIR}/${CONFDIR}/${myGlobalNml}  \
        -f ${BFMDIR}/src/BFM/proto \
        -t ${blddir} || exit

    if [[ ${MODE} == "STANDALONE" ]]; then
        # list files
        find ${BFMDIR}/src/BFM/General -name "*.?90" -print > BFM.lst
        find ${BFMDIR}/src/standalone -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/share -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Pel -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/PelBen -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Light -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Oxygen -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Forcing -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/CO2 -name "*.?90" -print >> BFM.lst

        #change netcdf path in compiler file
        if [ ${NETCDF} ]; then
            [ ${VERBOSE} ] && echo "changing netcd path!"
            sed -e "s,/usr/local,${NETCDF}," ${BFMDIR}/compilers/${ARCH} > ${blddir}/${ARCH}
        else
            cp ${BFMDIR}/compilers/${ARCH} ${blddir}/${ARCH}
        fi

        # Make makefile
        ${BFMDIR}/${SCRIPTSDIR}/${cmd_mkmf} \
            -c "${cppdefs}" \
            -o "-I${BFMDIR}/include -I${BFMDIR}/src/BFM/include" \
            -t "${blddir}/${ARCH}" \
            -p "${BFMDIR}/bin/bfm_standalone.x" \
            BFM.lst && echo ""

        # Link to the configuration file
        #ln -sf ${BFMDIR}/${CONFDIR}/${myGlobalDef} GlobalDefsBFM.model
    else
        # Generate the specific bfm.fcm include file for makenemo
        # some macros are default with NEMO
        FCMMacros="BFM_NEMO USEPACK BFM_NOPOINTERS ${CPPDEFS}"
        sed -e "s/_place_keys_/${FCMMacros}/" -e "s;_place_def_;${myGlobalConf};" \
            ${BFMDIR}/${SCRIPTS_CONFDIR}/Default_bfm.fcm > ${blddir}/bfm.fcm
        [ ${VERBOSE} ] && echo "Memory Layout generated in local folder: ${blddir}."

        # Move BFM Layout files to target folders 
        cp ${blddir}/*.F90 ${BFMDIR}/src/BFM/General
        mv ${BFMDIR}/src/BFM/General/init_var_bfm.F90 ${BFMDIR}/src/share
        cp ${blddir}/init_var_bfm.F90 ${BFMDIR}/src/share
        cp ${blddir}/INCLUDE.h ${BFMDIR}/src/BFM/include
        cp ${blddir}/bfm.fcm ${BFMDIR}/src/nemo
    fi
    echo "${PRESET} generation done!"
fi


# -----------------------------------------------------
# COMPILATION of executable
# -----------------------------------------------------


if [ ${CMP} ]; then
    if [ ! -d ${blddir} ]; then
        echo "ERROR: directory ${blddir} not exists"
        echo ${ERROR_MSG}
    fi
    cd ${blddir}

    if [[ ${MODE} == "STANDALONE" ]]; then
        if [ ${CLEAN} == 1 ]; then
            [ ${VERBOSE} ] && echo "Cleaning up ${PRESET}..."
            ${cmd_gmake} clean
        fi
        echo " "
        echo "Starting ${PRESET} compilation..."
        rm -rf ${BFMDIR}/bin/${BFMSTD}
        ${cmd_gmake}
        if [ ! -f ${BFMDIR}/bin/${BFMSTD} ]; then 
            echo "ERROR in ${PRESET} compilation!" ; 
            exit 1; 
        else
            echo " "
            echo "${PRESET} compilation done!"
            echo " "
        fi
    else
        cd ${NEMODIR}/NEMOGCM/CONFIG/

        if [ ${CLEAN} == 1 ]; then
            [ ${VERBOSE} ] && echo "Cleaning up ${PRESET}..."
            ./${cmd_mknemo} -n ${PRESET} -m ${ARCH} clean
        fi
        [ ${VERBOSE} ] && echo "Starting ${PRESET} compilation..."
        rm -rf ${NEMODIR}/NEMOGCM/CONFIG/${PRESET}/BLD/bin/${NEMOEXE}
        ./${cmd_mknemo} -n ${PRESET} -m ${ARCH} -e ${BFMDIR}/src/nemo -j ${PROC}
        if [ ! -f ${NEMODIR}/NEMOGCM/CONFIG/${PRESET}/BLD/bin/${NEMOEXE} ]; then 
            echo "ERROR in ${PRESET} compilation!" ; 
            exit 1; 
        else
            echo "${PRESET} compilation done!"
        fi
    fi
fi


# -----------------------------------------------------
# EXPERIMENT folder creation
# -----------------------------------------------------


if [ ${DEP} ]; then
    [ ${VERBOSE} ] && echo "creating Experiment ${PRESET}"

    if [ ! -d ${blddir} ]; then
        echo "ERROR: directory ${blddir} not exists"
        echo ${ERROR_MSG}
    fi
 

    #Copy Namelists
    exedir="${BFMDIR}/run/${EXP}"
    if [ ! -d ${exedir} ]; then 
        mkdir -p ${exedir};
        # copy and link namelist files
        if [ ${NMLDIR} ]; then cp ${NMLDIR}/*.nml ${exedir}/; 
        else cp ${blddir}/*.nml ${exedir}/; fi
    else
        echo "WARNING: directory ${exedir} exists (not copying namelist files)"
    fi

    #Copy nemo files and executable 
    if [[ ${MODE} == "STANDALONE" ]]; then
        ln -sf ${BFMDIR}/bin/${BFMSTD} ${exedir}/${BFMSTD}
        printf "Go to ${exedir} and execute command:\n\t./${BFMSTD}\n"
    else
        # copy and link necessary files
        cd "${BFMDIR}/${CONFDIR}/${PRESET}"
        cp ${NEMO_FILES} ${exedir}/
        ln -sf ${NEMODIR}/NEMOGCM/CONFIG/${PRESET}/BLD/bin/${NEMOEXE} ${exedir}/${NEMOEXE}
        #change values in runscript
        sed -e "s,_EXP_,${EXP},g"       \
            -e "s,_EXE_,${NEMOEXE},g" \
            -e "s,_VERBOSE_,${VERBOSE},g" \
            -e "s,_PRESET_,${PRESET},g" \
            -e "s,_QUEUE_,${QUEUE},g"   \
            -e "s,_PROC_,${PROC},g"     ${BFMDIR}/${SCRIPTS_CONFDIR}/runscript > ${exedir}/runscript_${EXP}
        printf "Go to ${exedir} and execute command:\n\tbsub < ./runscript_${EXP}\n"
    fi
fi
