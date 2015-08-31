#!/bin/bash -e

# DESCRIPTION
#   BFM Configuration manager
#
# AUTHORS
#   Esteban Gutierrez esteban.gutierrez@cmcc.it
#   Tomas Lovato tomas.lovato@cmcc.it
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

#logging configuration
LOGFILE=logfile_$$.log
LOGDIR="Logs"

#local paths
VER=""
TEMPDIR="build${VER}/tmp"
CONFDIR="build${VER}/configurations"
SCRIPTSDIR="build${VER}/scripts"
SCRIPTS_BIN="${SCRIPTSDIR}/bin"
SCRIPTS_PROTO="${SCRIPTSDIR}/proto"

#programs
MKMF="mkmf"
GMAKE="gmake"
PERL="perl"
GENCONF="generate_conf"
GENNML="generate_namelist"
MKNEMO="makenemo"       
NEMOEXE="nemo.exe"
#optional
NMLLIST="";

#nemo optional
NEMOSUB="";

#options
OPTIONS=(     MODE     CPPDEFS     ARCH     CLEAN     EXP     EXPDIR     EXECMD     VALGRIND     PROC     QUEUE     BFMEXE     NEMOSUB     EXPFILES     FORCING     NMLLIST     RUNPROTO     )
OPTIONS_USR=( mode_usr cppdefs_usr arch_usr clean_usr exp_usr expdir_usr execmd_usr valgrind_usr proc_usr queue_usr bfmexe_usr nemosub_usr expfiles_usr forcing_usr nmllist_usr runproto_usr )

#error message
ERROR_MSG="Execute $0 -h for help if you don't know what is going wrong. PLEASE read CAREFULLY before seeking help."

#----------------- USER CONFIGURATION DEFAULT VALUES -----------------
DEBUG="";
BFMDIR_DEFAULT="${PWD}/.."
MODE="STANDALONE"
CPPDEFS="BFM_STANDALONE INCLUDE_PELCO2 INCLUDE_DIAG"
PRESET="STANDALONE_PELAGIC"
ARCH="gfortran.inc"
PROC=8
PROC_CMP=8
EXP="EXP00"
QUEUE="poe_short"
BFMEXE="bfm_standalone.x"
CLEAN=1
NETCDF_DEFAULT="/usr"
MPICMD="mpirun.lsf"
RUNPROTO="runscript"
# --------------------------------------------------------------------


#print usage message 
usage(){
    more << EOF
NAME
    This script configure, compile and/or deploy the BFM model.

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
                  - POM1D
                  - NEMO (with NEMO)
                  - NEMO_3DVAR (with NEMO and 3DVAR)
       -k CPPDEFS
                  Key options to configure the model. (Default: "BFM_STANDALONE INCLUDE_PELCO2 INCLUDE_DIAG")                 
       -S NEMOSUB
                  NEMO sub-directories to add to the NEMO configuration
       -N NMLLIST
                  namelist file list containing namelists to copy to executable directory
       -a ARCH
                  Specify compilation Architecture file (Default: "gfortran.inc")
                  - For STANDALONE mode available archs, list dir : BFMDIR/compilers
                  - For NEMO mode available archs, execute command: NEMODIR/NEMOGCM/CONFIG/makenemo -h all
       -f
                  Fast mode. Dont execute "clean" command in compilation (clean is activated by default)
       -o BFMEXE OUTPUT
                  BFM executable name output. (Default: bfm_standalone.x)
    alternative DEPLOYMENT OPTIONS are:
       -x EXP
                  Name of the experiment for generation of the output folder (Default: "${EXP}")
       -F EXPFILES 
                  files (generated or not)  to copy to experiment directory
       -i FORCING 
                  Necessary forcings to execute the model. To specify several files, separate them by colon ';' and surround all by quotes '"'.
                  (If the path is relative, the ROOT will be ${BFMDIR})
       -D EXPDIR
                  Input dir where are files to copy to experiment directory (Default: "${TEMPDIR}/${PRESET}")
       -e EXECMD
                  Executable command to insert in runscript (Default for NEMO: "${MPICMD}", empty for others)
       -V VALGRIND
                  Executable valgrind command to insert in runscript
       -r PROC
                  Number of procs used for running. Default: 8
       -q QUEUE
                  Name of the queue number of procs used for running. Default
       -R RUNPROTO
                  Specify prototype for generating the experiment running script (Default: runscript)

ENVIRONMENT VARIABLES
       BFMDIR
                  Path to the root directory of BFM. (Default: "${BFMDIR_DEFAULT}" if environment variable is not defined)
       NETCDF
                  Path to netcdf library and header files. (Default: "${NETCDF_DEFAULT}" if environment variable is not defined)
       BFMDIR_RUN
                  Path to folder where experiments will be created. (Default: "${BFMDIR}/run" if environment variable is not defined)
       NEMODIR
                  Path to the root directory of NEMO. Must define this environment variable if you use NEMO
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
while getopts "hvgcdPp:m:k:N:S:a:fx:F:i:D:e:V:r:q:o:R:" opt; do
    case $opt in
      h ) usage;            rm ${LOGDIR}/${LOGFILE}         ; exit             ;;
      v )                   echo "verbose mode"             ; VERBOSE=1        ;;
      g ) [ ${VERBOSE} ] && echo "generation activated"     ; GEN=1            ;;
      c ) [ ${VERBOSE} ] && echo "compilation activated"    ; CMP=1            ;;
      d ) [ ${VERBOSE} ] && echo "deployment activated"     ; DEP=1            ;;
      P ) [ ${VERBOSE} ] && echo "list presets"             ; LIST=1           ;;
      p ) [ ${VERBOSE} ] && echo "preset $OPTARG"           ; PRESET=$OPTARG       ;;
      m ) [ ${VERBOSE} ] && echo "mode $OPTARG"             ; mode_usr=$OPTARG     ;;
      k ) [ ${VERBOSE} ] && echo "key options $OPTARG"      ; cppdefs_usr=$OPTARG  ;;
      N ) [ ${VERBOSE} ] && echo "NMLLIST=$OPTARG"          ; nmllist_usr=$OPTARG  ;;
      S ) [ ${VERBOSE} ] && echo "NEMOSUB=$OPTARG"          ; nemosub_usr=$OPTARG  ;;
      a ) [ ${VERBOSE} ] && echo "architecture $OPTARG"     ; arch_usr=$OPTARG     ;;
      f ) [ ${VERBOSE} ] && echo "fast mode activated"      ; clean_usr=0          ;;
      x ) [ ${VERBOSE} ] && echo "experiment $OPTARG"       ; exp_usr=$OPTARG      ;;
      F ) [ ${VERBOSE} ] && echo "experiment files $OPTARG" ; expfiles_usr=$OPTARG ;;
      i ) [ ${VERBOSE} ] && echo "forcing files $OPTARG"    ; forcing_usr=$OPTARG  ;;
      D ) [ ${VERBOSE} ] && echo "namelist dir $OPTARG"     ; expdir_usr=$OPTARG   ;;
      e ) [ ${VERBOSE} ] && echo "exe command $OPTARG"      ; execmd_usr=$OPTARG   ;;
      V ) [ ${VERBOSE} ] && echo "valgrind command $OPTARG" ; valgrind_usr=$OPTARG ;;
      r ) [ ${VERBOSE} ] && echo "n. procs $OPTARG"         ; proc_usr=$OPTARG     ;;
      q ) [ ${VERBOSE} ] && echo "queue name $OPTARG"       ; queue_usr=$OPTARG    ;;
      o ) [ ${VERBOSE} ] && echo "executable name $OPTARG"  ; bfmexe_usr=$OPTARG   ;;
      R ) [ ${VERBOSE} ] && echo "run proto  name $OPTARG"  ; runproto_usr=$OPTARG ;;
      * ) echo "option not recognized"                      ; exit                 ;;
    esac
done

#set bfm path
if [[ ! ${BFMDIR} ]]; then
    BFMDIR="${BFMDIR_DEFAULT}"
    [ ${VERBOSE} ] && echo "setting BFM path with default: ${BFMDIR}"
fi

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
    cmd_gennml="${GENNML}.pl -v"
    cmd_mknemo="${MKNEMO}"
else
    cmd_mkmf="${MKMF}"
    cmd_gmake="${GMAKE} -s"
    cmd_gen="${GENCONF}.pl"
    cmd_gennml="${GENNML}.pl"
    cmd_mknemo="${MKNEMO} -v0"
fi
if [ $DEBUG ]; then
    cmd_gen=${cmd_gen}" -d"
    cmd_gennml=${cmd_gennml}" -d"
fi

myGlobalConf="${PRESET}/configuration"
myGlobalMem="${PRESET}/layout"
myGlobalNml="${PRESET}/namelists_bfm"

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
#specify preset dir
presetdir="${BFMDIR}/${CONFDIR}/${PRESET}"


#Check some optional parameter values
if [[ ! $NEMODIR && ( "$MODE" == "NEMO" || "$MODE" == "NEMO_3DVAR" ) ]]; then
    echo "ERROR: NEMODIR not specified in NEMO mode"
    echo ${ERROR_MSG}
    exit
fi
if [[ "$MODE" != "STANDALONE" && "$MODE" != "POM1D" && "$MODE" != "NEMO" && "$MODE" != "NEMO_3DVAR" ]]; then 
    echo "ERROR: MODE value not valid ($MODE). Available values are: STANDALONE or POM1D or NEMO or NEMO_3DVAR."
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
         echo "ERROR: ${BFMDIR}/${CONFDIR}/${myGlobalMem} not exists"
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
    cp ${presetdir}/* ${blddir}
    
    #add -D to cppdefs
    cppdefs=`echo ${CPPDEFS} | sed -e 's/\([a-zA-Z_0-9]*\)/-D\1/g'`

    # generate BFM Memory Layout files and namelists
    ${PERL} -I${BFMDIR}/${SCRIPTS_BIN}/ ${BFMDIR}/${SCRIPTS_BIN}/${cmd_gen} \
        ${cppdefs} \
        -r ${BFMDIR}/${CONFDIR}/${myGlobalMem}  \
        -n ${BFMDIR}/${CONFDIR}/${myGlobalNml}  \
        -f ${BFMDIR}/${SCRIPTS_PROTO} \
        -t ${blddir} || exit

    # generate other namelists specified in the argument NMLLIST
    if [ "${NMLLIST}" ]; then
        ${PERL} -I${BFMDIR}/${SCRIPTS_BIN}/ ${BFMDIR}/${SCRIPTS_BIN}/${cmd_gennml} \
            -i ${presetdir} \
            -n "${NMLLIST}"  \
            -o ${blddir} || exit
    fi


    if [[ ${MODE} == "STANDALONE" || "$MODE" == "POM1D" ]]; then
        # list files
        find ${BFMDIR}/src/BFM/General -name "*.?90" -print > BFM.lst
        if [[ ${MODE} == "STANDALONE" ]]; then
            find ${BFMDIR}/src/standalone -name "*.?90" -print >> BFM.lst
        else
            find ${BFMDIR}/src/pom -name "*.?90"  -print >> BFM.lst
            find ${BFMDIR}/src/pom -name "*.[fF]" -print >> BFM.lst
        fi
        find ${BFMDIR}/src/share -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Pel -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Light -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Oxygen -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/Forcing -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/CO2 -name "*.?90" -print >> BFM.lst
        find ${BFMDIR}/src/BFM/PelBen -name "*.?90" -print >> BFM.lst
        if echo "$cppdefs" | grep -q "\-DINCLUDE_BEN" ; then
            [ $VERBOSE ] && echo "include BENTHIC in BFM.lst"        
            find ${BFMDIR}/src/BFM/Ben -name "*.?90" -print >> BFM.lst
            find ${BFMDIR}/src/BFM/Bennut -name "*.?90" -print >> BFM.lst
        fi
        if echo "$cppdefs" | grep -q "\-DINCLUDE_SEAICE" ; then
            [ $VERBOSE ] && echo "include SEAICE in BFM.lst"
            find ${BFMDIR}/src/BFM/Seaice -name "*.?90" -print >> BFM.lst
        fi

        #set netcdf path in compiler file
        if [ ${NETCDF} ]; then
            [ ${VERBOSE} ] && echo "setting netcd path with environment variable: ${NETCDF}"
            sed -e "s,\${NETCDF},${NETCDF}," ${BFMDIR}/compilers/${ARCH} > ${blddir}/${ARCH}
        else
            [ ${VERBOSE} ] && echo "setting netcd path with default: ${NETCDF_DEFAULT}"
            sed -e "s,\${NETCDF},${NETCDF_DEFAULT}," ${BFMDIR}/compilers/${ARCH} > ${blddir}/${ARCH}
        fi

        if [ ${VERBOSE} ]; then
            echo "Executing: "
            echo "${BFMDIR}/${SCRIPTS_BIN}/${cmd_mkmf}"
            echo "    -c \"${cppdefs}\" "
            echo "    -o \"-I${BFMDIR}/include -I${BFMDIR}/src/BFM/include\" "
            echo "    -t \"${blddir}/${ARCH}\" "
            echo "    -p \"${BFMDIR}/bin/${BFMEXE}\" "
            echo "    BFM.lst && echo \" "
        fi

        # Make makefile
        ${BFMDIR}/${SCRIPTS_BIN}/${cmd_mkmf} \
            -c "${cppdefs}" \
            -o "-I${BFMDIR}/include -I${BFMDIR}/src/BFM/include" \
            -t "${blddir}/${ARCH}" \
            -p "${BFMDIR}/bin/${BFMEXE}" \
            BFM.lst && echo ""


        # Move BFM Layout files to target folders 
        cp ${blddir}/*.F90 ${BFMDIR}/src/BFM/General
        mv ${BFMDIR}/src/BFM/General/init_var_bfm.F90 ${BFMDIR}/src/share
        cp ${blddir}/init_var_bfm.F90 ${BFMDIR}/src/share
        cp ${blddir}/INCLUDE.h ${BFMDIR}/src/BFM/include
    elif [[ "$MODE" == "NEMO" || "$MODE" == "NEMO_3DVAR" ]]; then 
        #Generate NEMO configuration with subdirs
        if [ ! -d ${NEMODIR}/NEMOGCM/CONFIG/${PRESET} ]; then
            if [ "$NEMOSUB" ] ; then
                ${NEMODIR}/NEMOGCM/CONFIG/${cmd_mknemo} -j0 -n ${PRESET} -m ${ARCH} -d "${NEMOSUB}"
                if [ ! -f "${presetdir}/cpp_${PRESET}.fcm" ]; then 
                    echo "ERROR: cpp_${PRESET}.fcm must exist in ${presetdir}" 
                    exit
                fi
            else
                echo "ERROR: NEMO configuration not exists. If you want to create a NEMO configuration, you have to specify NEMOSUB option"
                echo ${ERROR_MSG}
                exit
            fi
        fi

        #if cpp_PRESET exists in configuration folder => copy to nemo preset folder
        if [ -f "${presetdir}/cpp_${PRESET}.fcm" ]; then 
            cp "${presetdir}/cpp_${PRESET}.fcm" "${NEMODIR}/NEMOGCM/CONFIG/${PRESET}"
        fi
        #substitute BFM/src/nemo path in cpp according to the version
        [ ${VERBOSE} ] && echo "Changing nemo version to: \"$VER\""
        sed -ie "s;inc.*;inc \$BFMDIR/src/nemo${VER}/bfm\.fcm;" ${NEMODIR}/NEMOGCM/CONFIG/${PRESET}/cpp_${PRESET}.fcm


        # Generate the specific bfm.fcm include file for makenemo
        # some macros are default with NEMO
        if echo "$cppdefs" | grep -q "\-DINCLUDE_BEN" ; then
            [ $VERBOSE ] && echo "include BENTHIC in bfm.fcm"
            FCMBen="src::bfm::ben             ${BFMDIR}/src/BFM/Ben\&src::bfm::bennut          $BFMDIR/src/BFM/Bennut"
        fi
        if echo "$cppdefs" | grep -q "\-DINCLUDE_SEAICE" ; then
            [ $VERBOSE ] && echo "include SEAICE in bfm.fcm"
            FCMIce="src::bfm::seaice          ${BFMDIR}/src/BFM/Seaice"
        fi
        FCMMacros="BFM_NEMO BFM_NOPOINTERS ${CPPDEFS}"
        sed -e "s/_place_keys_/${FCMMacros}/" -e "s;_place_def_;${myGlobalConf};" \
            -e "s;_place_ben_;${FCMBen};g"     -e "s;_place_ice_;${FCMIce};"       \
            ${BFMDIR}/${SCRIPTS_PROTO}/Default_bfm.fcm | tr "\&" "\n" > ${blddir}/bfm.fcm
        [ ${VERBOSE} ] && echo "Memory Layout generated in local folder: ${blddir}"

        # Move BFM Layout files to target folders 
        cp ${blddir}/*.F90 ${BFMDIR}/src/BFM/General
        mv ${BFMDIR}/src/BFM/General/init_var_bfm.F90 ${BFMDIR}/src/share
        cp ${blddir}/init_var_bfm.F90 ${BFMDIR}/src/share
        cp ${blddir}/INCLUDE.h ${BFMDIR}/src/BFM/include
        cp ${blddir}/bfm.fcm ${BFMDIR}/src/nemo${VER}
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

    if [[ ${MODE} == "STANDALONE" || "$MODE" == "POM1D" ]]; then
        if [ ${CLEAN} == 1 ]; then
            [ ${VERBOSE} ] && echo "Cleaning up ${PRESET}..."
            [ ${VERBOSE} ] && echo "Command: ${cmd_gmake} clean"
            ${cmd_gmake} clean
        fi
        echo " "
        echo "Starting ${PRESET} compilation..."
        rm -rf ${BFMDIR}/bin/${BFMEXE}
        [ ${VERBOSE} ] && echo "Command: ${cmd_gmake}"
        ${cmd_gmake}
        if [ ! -f ${BFMDIR}/bin/${BFMEXE} ]; then 
            echo "ERROR in ${PRESET} compilation!" ; 
            exit 1; 
        else
            echo " "
            echo "${PRESET} compilation done!"
            echo " "
        fi
    elif [[ ${MODE} == "NEMO" || "$MODE" == "NEMO_3DVAR" ]]; then

        if [ ${CLEAN} == 1 ]; then
            [ ${VERBOSE} ] && echo "Cleaning up ${PRESET}..."
            [ ${VERBOSE} ] && echo "Command: ${NEMODIR}/NEMOGCM/CONFIG/${cmd_mknemo} -n ${PRESET} -m ${ARCH} clean"
            ${NEMODIR}/NEMOGCM/CONFIG/${cmd_mknemo} -n ${PRESET} -m ${ARCH} clean
        fi
        [ ${VERBOSE} ] && echo "Starting ${PRESET} compilation..."
        rm -rf ${NEMODIR}/NEMOGCM/CONFIG/${PRESET}/BLD/bin/${NEMOEXE}
        if [[ "$MODE" == "NEMO" ]]; then 
            [ ${VERBOSE} ] && echo "Command: ${NEMODIR}/NEMOGCM/CONFIG/${cmd_mknemo} -n ${PRESET} -m ${ARCH} -e ${BFMDIR}/src/nemo${VER} -j ${PROC_CMP}"
            ${NEMODIR}/NEMOGCM/CONFIG/${cmd_mknemo} -n ${PRESET} -m ${ARCH} -e ${BFMDIR}/src/nemo${VER} -j ${PROC_CMP}
        else
            [ ${VERBOSE} ] && echo "Command: ${NEMODIR}/NEMOGCM/CONFIG/${cmd_mknemo} -n ${PRESET} -m ${ARCH} -e \"${BFMDIR}/src/nemo${VER};${NEMODIR}/3DVAR\" -j ${PROC_CMP}"
            ${NEMODIR}/NEMOGCM/CONFIG/${cmd_mknemo} -n ${PRESET} -m ${ARCH} -e "${BFMDIR}/src/nemo${VER};${NEMODIR}/3DVAR" -j ${PROC_CMP}
        fi
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

    #set run dir path
    if [ ${BFMDIR_RUN} ]; then
        exedir="${BFMDIR_RUN}/${EXP}"
        [ ${VERBOSE} ] && echo "setting run dir path with environment variable: ${exedir}"
    else
        exedir="${BFMDIR}/run/${EXP}"
        [ ${VERBOSE} ] && echo "setting run dir path with default: ${exedir}"
    fi

    #Copy Namelists and other files
    if [ ! -d ${exedir} ]; then 
        if [ ! -d ${blddir} ]; then
            echo "ERROR: directory ${blddir} not exists"
            echo ${ERROR_MSG}
            exit
        fi

        mkdir -p ${exedir};

        #copy files from blddir to exedir
        cd ${blddir}
        cp *.nml ${exedir}/
        if [ "${EXPFILES}" ]; then cp ${EXPFILES} ${exedir}/; fi
        if [ "${EXPDIR}"   ]; then cp ${EXPDIR}/* ${exedir}/; fi

        #link forcing files to exedir
        if [ "${FORCING}" ]; then
            [ ${VERBOSE} ] && echo "linking Forcings: ${FORCING}"
            forcing_list=(${FORCING//;/ })
            for idx_for in "${!forcing_list[@]}"; do
                #if is relative path, start with BFMDIR
                if [[ ${forcing_list[$idx_for]} =~ ^\s*\/ ]]; then
                    ln -sf ${forcing_list[$idx_for]} ${exedir}/
                else
                    ln -sf ${BFMDIR}/${forcing_list[$idx_for]} ${exedir}/
                fi
            done
        fi

        #link reference nemo files from the shared directory
        if [[ "$MODE" == "NEMO" || "$MODE" == "NEMO_3DVAR" ]] && [ -d ${NEMODIR}/NEMOGCM/CONFIG/SHARED ]; then 
           ln -sf ${NEMODIR}/NEMOGCM/CONFIG/SHARED/*_ref ${exedir}/;
           ln -sf ${NEMODIR}/NEMOGCM/CONFIG/SHARED/*_def.xml ${exedir}/;
        fi
    else
        echo "WARNING: directory ${exedir} exists (not overwriting namelists)"
    fi

    #Copy executable and insert exe and valgrind commands
    if [[ ${MODE} == "STANDALONE" || "$MODE" == "POM1D" ]]; then
        ln -sf ${BFMDIR}/bin/${BFMEXE} ${exedir}/${BFMEXE}
        if [ "${EXECMD}" ]; then 
            execmd="${EXECMD} ${VALGRIND} ./${BFMEXE}"
        else
            execmd="${VALGRIND} ./${BFMEXE}"
        fi
    elif [[ "$MODE" == "NEMO" || "$MODE" == "NEMO_3DVAR" ]]; then 
        ln -sf ${NEMODIR}/NEMOGCM/CONFIG/${PRESET}/BLD/bin/${NEMOEXE} ${exedir}/${NEMOEXE}
        if [ "${EXECMD}" ]; then 
            execmd="${EXECMD} ${VALGRIND} ./${NEMOEXE}"
        else 
            execmd="${MPICMD} ${VALGRIND} ./${NEMOEXE}"
        fi
    fi

    #change values in runscript
    sed -e "s,_EXP_,${EXP},g"         \
        -e "s,_EXE_,${execmd},g"      \
        -e "s,_VERBOSE_,${VERBOSE},g" \
        -e "s,_PRESET_,${PRESET},g"   \
        -e "s,_QUEUE_,${QUEUE},g"     \
        -e "s,_PROC_,${PROC},g"     ${BFMDIR}/${SCRIPTS_PROTO}/${RUNPROTO} > ${exedir}/${RUNPROTO}_${EXP}
    printf "Go to ${exedir} and execute command:\n\tbsub < ./${RUNPROTO}_${EXP}\n"

fi
