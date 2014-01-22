#!/usr/bin/perl -w

# DESCRIPTION
#   Execute tests for BFM model
#
# AUTHORS
#   Esteban Gutierrez esteban.gutierrez@cmcc.it
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

use 5.008_002;

use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

use lib './scripts';
use get_configuration; #parse configuration file
use classes;

#check if global variables are defined
if( ! length ${ENV{'BFMDIR'}} ) { print "ERROR, \$BFMDIR must be defined\n"; exit; }

#Fix values
my $BASE_DIR = "${ENV{'BFMDIR'}}/tools/bntest";
my $CONF_DIR = "${BASE_DIR}/configurations";
#Default values
my $temp_dir = "$BASE_DIR/tmp";
my $out_dir  = "$BASE_DIR/out";
my $list     = 0;
my $verbose  = 0;
my $help     = 0;
#Must parameters
my ($input_preset);
#configuration parameters
my (%user);

sub usage(){
    print "usage: $0 {-p [preset_file] -h -P} [-o [output_dir]] [-t [temporal_dir]] [-v]\n\n";
    print "This script execute tests for BFM model\n\n";
    print "MUST specify one of these options:\n";
    print "\t-f [preset_file]  test configuration preset\n";
    print "\t-P                list presets available\n";
    print "\t-h                print help message\n";
    print "ALTERNATIVE options are:\n";
    print "\t-o [output_dir]   output dir for generated files\n";
    print "\t-t [temporal_dir] output dir for temporal files\n";
    print "\t-v                verbose mode\n";

}

#check input arguments
use Getopt::Long qw(:config bundling noignorecase); # for getopts compat
GetOptions(
    'p=s'  => \$input_preset,
    't=s'  => \$temp_dir,
    'o=s'  => \$out_dir,  
    'v'    => \$verbose,
    'P'    => \$list,
    'h'    => \$help,
    ) or &usage() && exit;
if ( $list ){ 
    opendir my($dirlist), ${CONF_DIR} or die "Couldn't open dir '${CONF_DIR}': $!";
    my @presets = grep !/^\.\.?$/ && -d ${CONF_DIR}, readdir $dirlist;
    close $dirlist;
    print join("\n", @presets) . "\n";
    exit; 
}
if ( $help ){ &usage(); exit; }
if ( !$input_preset  ){ &usage(); exit; }

#create dirs
if( ! -d $temp_dir ){ 
    if($verbose){ print "Creating Temporal dir: $temp_dir\n"; }
    mkdir $temp_dir or die "Unable to create $temp_dir\n";
}

#read configuration file
if( $verbose ){ print "Reading configuration...\n"; }
my $lst_test = get_configuration("${CONF_DIR}/${input_preset}/configuration", $verbose);

foreach my $test (@$lst_test){
    if($verbose){ print "Executing: " . $test->getName() . "\n"; }
}

if( $verbose ){ print "Test finished\n"; }
