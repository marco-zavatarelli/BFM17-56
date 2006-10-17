#!/bin/bash
# This script sets the BFM environmental variables used in the Makefile.
# Modify the following according to your system settings
# The script must be executed before compilation in the current shell
# (e.g. > . ./bfm_env.sh) or added to the .bashrc script in $HOME

echo "Setting BFM environment variables for host" $HOSTNAME
# BFM environmental variables
export BFMDIR=$HOME/_BFMRELEASE_
echo "Setting the BFM ROOTDIR to "$BFMDIR

# GOTM environmental variables #
export GOTMDIR=$HOME/gotm
echo "Setting the GOTM ROOTDIR to "$GOTMDIR
export GOTM_CASES=$HOME/bfm-run/gotm
echo "Setting the GOTM CASES to "$GOTMDIR

# NEMO environmental variables #

# POM environmental variables #

# Compilation #
export FORTRAN_COMPILER=IFORT
echo "Current Fortran compiler is "$FORTRAN_COMPILER
export COMPILATION_MODE=production
echo "Current compilation mode is "$COMPILATION_MODE
export NETCDFINC=/opt/netcdf/include
export NETCDFLIBDIR=/opt/netcdf/lib
echo "Linking NetCDF library from:"
echo $NETCDFLIBDIR
echo $NETCDFINC
