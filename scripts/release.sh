#!/bin/sh
# Script to create an official release of BFM
# Inputs:
# release.sh bfm_version release_version release_type

[ "$USER" = "vichi" ] || { echo "Only vichi can make new releases" ; exit 1; }

bfm_version=$1
release_version=$2
release_type=$3

release_name=bfm-$bfm_version-system-$release_version
dev_dir=`pwd`
tmp_dir=$HOME/tmp
tarfile=$release_name.tgz

TRUNK="http://www.bo.ingv.it/svn/bfm/trunk" 
CURRENT="http://www.bo.ingv.it/svn/bfm/tags/current_release" 
RELEASE="http://www.bo.ingv.it/svn/bfm/tags/releases/" 

# Create the tagged release
svn cp --message "Created a new release" $TRUNK $RELEASE/$release_name
if [ "$release_type" = "current" ] ; then
  svn cp --message "Made it current" $RELEASE/$release_name $CURRENT
fi
# Create the tgz file
svn export $RELEASE/$release_name $tmp_dir/
cd $tmp_dir 
# Create the scripts
sed -e "s/_BFMRELEASE_/$release_name/" $dev_dir/scripts/bfm_env_proto.sh > $release_name/bfm_env.sh
sed -e "s/_BFMRELEASE_/$release_name/" $dev_dir/scripts/bfm_env_proto.csh > $release_name/bfm_env.csh
rm -rf $release_name/scripts
tar -zcvf $tarfile $release_name && rm -rf $release_name

exit 0
