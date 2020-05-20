#!/bin/sh
# Script to define local paths, compile and execute BFM_POM 1D,
# change restart files names

export BFMDIR="$(pwd)"
export NETCDF="/usr/local/netcdf"
export BFMDIR_RUN="$BFMDIR/bfm_run"

echo "BFMDIR is "${BFMDIR}

mkdir -p $BFMDIR_RUN/bfm17_pom1d
cd $BFMDIR_RUN/bfm17_pom1d
rm -rf  pom.exe pom_restart* bfm_restart*

cd $BFMDIR/build

./bfm_configure.sh -gcd -p BFM17_POM1D BFM17

cd $BFMDIR_RUN/bfm17_pom1d
cp $BFMDIR/build/configurations/BFM17_POM1D/*.nml .

echo "Namelists copied"

if [ "x$ID_run" != "x" ] && [ "y$next_ID_run" != "y" ] ; then
   ln -s pom_restart_$ID_run fort.70 #in
   ln -s pom_restart_$next_ID_run fort.71 #out
   
   ln -s bfm_restart_$ID_run in_bfm_restart.nc #in
   ln -s bfm_restart_$next_ID_run out_bfm_restart.nc #out
   echo "restart files linked"
fi

echo "Start model simulation"

time ./pom.exe
