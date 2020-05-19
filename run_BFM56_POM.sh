#!/bin/sh
# Script to define local paths, compile and execute BFM_POM 1D,
# change restart files names

export BFMDIR="$(pwd)"
export NETCDF="/usr/local/netcdf"
export BFMDIR_RUN="$BFMDIR/bfm_run"

mkdir $BFMDIR_RUN/bfm56_pom1d
cd $BFMDIR_RUN/bfm56_pom1d
rm -rf  pom.exe pom_restart* bfm_restart*

cd $BFMDIR/build

./bfm_configure.sh -gdc -p BFM56_POM1D

cd $BFMDIR_RUN/bfm56_pom1d
cp $BFMDIR/build/configurations/BFM56_POM1D/* .

ln -s pom_restart_$ID_run fort.70 #in
ln -s pom_restart_$next_ID_run fort.71 #out

ln -s bfm_restart_$ID_run in_bfm_restart.nc #in
ln -s bfm_restart_$next_ID_run out_bfm_restart.nc #out

time ./pom.exe
