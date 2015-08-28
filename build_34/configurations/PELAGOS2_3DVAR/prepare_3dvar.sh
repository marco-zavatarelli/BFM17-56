###############################################################################
##
##      Script to run 3DVAR
##
##      Authors:
##      Luca Visinelli (CMCC)
##      Pierluigi Di Pietro (INGV)
##
##      CMCC Centro Euro-Mediterraneo sui Cambiamenti Climatici
##      INGV Istituto Nazionale di Geofisica e Vulconologia
##
###############################################################################
#!/bin/ksh

set -x

##############################################################################
#
# NEMO INTERFACE 
#
##############################################################################

# Configuration of the 3DVAR:
# 102: Data Interpolation
# 103: Minimization and computation of misfits
nconf=_nconfig_

nemo_dom=_domains_         # Total number of processors
nameEOF=_nameEOF_         # Name of the EOF files

initialdate=_inidate_      # Initial date
currentdate=_date_         # Date of the first day in the current run
assimdate=_obsdate_        # Date at the center of the assimilation window
assim_hours=_assimhour_    # Number of hours for each assimilation

EXECDIR=_dir3dvar_         # Folder where the 3DVAR is run
en3dir=_INSITUdir_         # Folder containing the physics (T&S) data
BGCdir=_INSITU_BGCdir_     # Folder containing the BGC data

OBSDIC=_OBSDIC_
OBSALK=_OBSALK_
OBSCHL=_OBSCHL_
OBSOXY=_OBSOXY_
OBSpCO=_OBSpCO_

assim_hhmm=0000

###############################################################################
#
# END OF NEMO INTERFACE 
#
###############################################################################

export OMP_NUM_THREADS=_domains_
export OMP_STACKSIZE="1G"
export MPIMULTITASKMIX="ON"

yymm=`echo ${assimdate} | cut -c 1-6`
# current year and month
yy=`echo ${assimdate} | cut -c 1-4`
mm=`echo ${assimdate} | cut -c 5-6`
# previous year and month
mmm1=$(echo "$mm-1"|bc)
yym1=$(echo "$yy-1"|bc)
# future year and month
mmp1=$(echo "$mm+1"|bc)
yyp1=$(echo "$yy+1"|bc)

if [ $mmm1 -le 9 ]; then mmm1=0$mmm1; fi
if [ $mmp1 -le 9 ]; then mmp1=0$mmp1; fi

case $mm in 
	01) prev=${yym1}12
            nex=${yy}${mmp1};;
	12) prev=${yy}${mmm1}
	    nex=${yyp1}01;;
	02|03|04|05|06|07|08|09|10|11) prev=${yy}${mmm1}
	   nex=${yy}${mmp1};;
	*) echo failed at prev/next date evaluation ;;
esac

echo getting data from following months: $prev $yymm $nex

if [ $nconf -eq 102 ]
then
  echo 3DVAR operates as data interpolator: nconf = ${nconf}

#------------------------------------------------------------------------------
# Link physics data
#------------------------------------------------------------------------------

  rm INS_*

  ln -sf $en3dir/EN3_v2a_NoCWT_WijffelsTable1XBTCorr_Profiles_${prev}.nc ./INS_01
  ln -sf $en3dir/EN3_v2a_NoCWT_WijffelsTable1XBTCorr_Profiles_${yymm}.nc ./INS_02
  ln -sf $en3dir/EN3_v2a_NoCWT_WijffelsTable1XBTCorr_Profiles_${nex}.nc ./INS_03 
#------------------------------------------------------------------------------
# Link BGC data
#------------------------------------------------------------------------------

    if [ "${OBSDIC}" = "true" ]; then

      echo DIC assimilated

      fil1=${BGCdir}/DIC/GLODAP_CARINA_DIC_${prev}.nc
      if [ -e $fil1 ]; then
        ln -sf $fil1 ./INS_DIC_01
      else
        echo missing INS_DIC_01 file
      fi

      fil2=${BGCdir}/DIC/GLODAP_CARINA_DIC_${yymm}.nc
      if [ -e $fil2 ]; then
        ln -sf $fil2 ./INS_DIC_02
      else
        echo missing INS_DIC_02 file
      fi

      fil3=${BGCdir}/DIC/GLODAP_CARINA_DIC_${nex}.nc
      if [ -e $fil3 ]; then
        ln -sf $fil3 ./INS_DIC_03
      else
        echo missing INS_DIC_03 file
      fi

    fi


    if [ "${OBSALK}" = "true" ]; then

      echo ALK assimilated 

      fil1=${BGCdir}/ALK/GLODAP_CARINA_ALK_${prev}.nc
      if [ -e $fil1 ]; then
        ln -sf $fil1 ./INS_ALK_01
      else
        echo missing INS_ALK_01 file
      fi

      fil2=${BGCdir}/ALK/GLODAP_CARINA_ALK_${yymm}.nc
      if [ -e $fil2 ]; then
        ln -sf $fil2 ./INS_ALK_02
      else
        echo missing INS_ALK_02 file
      fi

      fil3=${BGCdir}/ALK/GLODAP_CARINA_ALK_${nex}.nc
      if [ -e $fil3 ]; then
        ln -sf $fil3 ./INS_ALK_03
      else
        echo missing INS_ALK_03 file
      fi

    fi

#------------------------------------------------------------------------------
# Link SLA data
#------------------------------------------------------------------------------

  rm ERS* ENVISAT* GFO* JASON* TP*

  for sat in e1 e2 en enn g2 j1 j1n j2 tp tpn
  do
    for tim in prev ${yymm} nex
    do
      dat=0;
      files=$(ls /users/home/ans049/DATA/SLA/$sat/$tim/*vfec*)
      for fil in $files
      do
        if [ -f $fil ]
        then
	  let  dat+=1
	  case $sat in
	    e1) 	ln -s $fil ERS10$dat;;
	    e2) 	ln -s $fil ERS20$dat;;
	    en|enn) 	ln -s $fil ENVISAT0$dat;;
	    g2) 	ln -s $fil GFO0$dat;;
	    j1|j1n) 	ln -s $fil JASON10$dat;;
	    j2) 	ln -s $fil JASON20$dat;;
	    tp|tpn) 	ln -s $fil TP0$dat;;
	    *) 	exit 'error';;
	  esac
        else
          echo file not present for $sat in $tim
        fi
      done
    done
  done

#------------------------------------------------------------------------------
# Run 3DVAR executable and produce obs_tab.dat
#------------------------------------------------------------------------------
  time ./MASTER -d ${assimdate} -t $assim_hhmm -w $assim_hours 

  mkdir ${assimdate}
  regs=0
  while [[ $regs -lt $nemo_dom ]]; do 
    mv INSOBS_$(printf "%04d" $regs).NC ./${assimdate}/INSOBS_$(printf "%04d" $regs)_${assimdate}.NC
    mv SLAOBS_$(printf "%04d" $regs).NC ./${assimdate}/SLAOBS_$(printf "%04d" $regs)_${assimdate}.NC
    regs=$(echo "$regs + 1"|bc)
  done 
  mv OBS_TAB.DAT OBS_TAB_${assimdate}.DAT
  mv LOG_3DVAR LOG_3DVAR_102_${assimdate}
  echo if present, nsobs, slaobs and obs_tab moved to obs_tab.dat${assimdate}

elif [ $nconf -eq 103 ]
then

  echo 3dvar operating as minimizer: nconf = $nconf

  if [ ${currentdate} -ne ${initialdate} ]
  then
    case $mm in
      01|02|03) ln -sf ${nameEOF}win.nc EOF.nc ;;
      04|05|06) ln -sf ${nameEOF}spr.nc EOF.nc ;;
      07|08|09) ln -sf ${nameEOF}sum.nc EOF.nc ;;
      10|11|12) ln -sf ${nameEOF}aut.nc EOF.nc ;;
           *)  exit error in linking eof ;;
    esac
    if [ ! -e OBS_TAB_${currentdate}.DAT ]
    then
      echo OBS_TAB file related to previous assim not found ${currentdate}
    else
      ln -sf OBS_TAB_${currentdate}.DAT OBS_TAB.DAT
    fi

    rm -rf INSMIS_0*.NC
    regs=0

    # NB The processor numbering in obsmisfit files starts from 1 instead of 0
    while [[ $regs -le $nemo_dom ]]; do
      if [ ! -e $EXECDIR/${currentdate}/INSMIS_$(printf "%04d" $regs)_${currentdate}.NC ]
      then
        echo insitu misfit file ${regs} not found ${currentdate}
      else
        ln -sf $EXECDIR/${currentdate}/INSMIS_$(printf "%04d" $regs)_${currentdate}.NC INSMIS_$(printf "%04d" $regs).NC
      fi
      if [ ! -e $EXECDIR/${currentdate}/SLAMIS_$(printf "%04d" $regs)_${currentdate}.NC ]
      then
        echo SLA misfit file ${regs} not found ${currentdate}
      else
        ln -sf $EXECDIR/${currentdate}/SLAMIS_$(printf "%04d" $regs)_${currentdate}.NC SLAMIS_$(printf "%04d" $regs).NC
      fi
      regs=$(echo "$regs + 1"|bc)
    done 

    ./shuffle_obs.x $nemo_dom  || exit -12


#------------------------------------------------------------------------------
# Run 3DVAR executable and produce ANINCR.nc
#------------------------------------------------------------------------------ 
    time ./MASTER -d ${currentdate} -t $assim_hhmm -w $assim_hours

    mv OBS_TAB.DAT OBS_TAB_corr_${currentdate}.DAT
    mv ANINCR.NC ANINCR_${currentdate}.NC
    mv LOG_3DVAR LOG_3DVAR_103_${currentdate}
    mv OBSSTAT.NC OBSSTAT_${currentdate}.NC
    echo correction produced and moved to ANINCR_${currentdate}.NC
  else  # Begin of the run
    echo no correction at run origin
    exit
  fi
else
  echo bad option $nconf
fi


#------------------------------------------------------------------------------
# EPILOGUE
#------------------------------------------------------------------------------
exit
#------------------------------------------------------------------------------
