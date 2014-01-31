#!/usr/bin/perl -w

# DESCRIPTION
#   Testing environment for BFM model
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
use bfm_test_modules; #subroutines and functions for bfm_test
use classes;

#check if global variables are defined
if( ! length ${ENV{'BFMDIR'}} ) { print "ERROR, \$BFMDIR must be defined\n"; exit; }

#Fix values
my $BFMDIR     = "${ENV{'BFMDIR'}}";
my $BUILD_DIR  = "${BFMDIR}/build";
my $BASE_DIR   = "${BFMDIR}/tools/bfm_test";
my $CONF_DIR   = "${BASE_DIR}/configurations";
my $BFM_EXE    = "bfm_configure.sh";
#Default values
my ($generation, $execution, $analysis) = 0;
my $temp_dir  = "$BASE_DIR/tmp";
my $list      = 0;
my $overwrite = 0;
my $verbose   = 0;
my $help      = 0;
#Must parameters
my ($input_preset);
#configuration parameters
my (%user);

sub usage(){
    print "usage: $0 {-p [preset_file] -h -P} [-t [temporal_dir] -o -v]\n\n";
    print "This script generate, execute and analyze configurations for BFM model\n\n";
    print "INFORMATIVE options:\n";
    print "\t-P                list presets available\n";
    print "\t-h                print help message\n";
    print "for running you MUST specify one of these options:\n";
    print "\t-p [preset_file]  test configuration preset\n";
    print "\t-g                Generate preset\n";
    print "\t-x                Execute preset\n";
    print "\t-a                Analyze preset\n";
    print "ALTERNATIVE options are:\n";
    print "\t-t [temporal_dir] output dir for temporal files\n";
    print "\t-o                overwrite generated directories\n";
    print "\t-v                verbose mode\n";
    print "\n";
}

#check input arguments
use Getopt::Long qw(:config bundling noignorecase); # for getopts compat
GetOptions(
    'p=s'  => \$input_preset,
    'g'    => \$generation,
    'x'    => \$execution,
    'a'    => \$analysis,
    't=s'  => \$temp_dir,
    'o'    => \$overwrite,
    'v'    => \$verbose,
    'P'    => \$list,
    'h'    => \$help,
    ) or &usage() && exit;
if ( $list ){ # list dirs inside bfm configuration
    opendir my($dirlist), ${CONF_DIR} or die "Couldn't open dir '${CONF_DIR}': $!";
    my @presets = grep !/^\.\.?$/ && -d ${CONF_DIR}, readdir $dirlist;
    close $dirlist;
    print join("\n", @presets) . "\n";
    exit;
}
if ( $help ){ &usage(); exit; }
if ( !$input_preset || !($generation || $execution || $analysis) ){ &usage(); exit; }

#create dirs
if( ! -d $temp_dir ){ 
    if($verbose){ print "Creating Temporal dir: $temp_dir\n"; }
    mkdir $temp_dir or die "Unable to create $temp_dir\n";
}

#read configuration file
if( $verbose ){ print "Reading configuration...\n"; }
my $lst_test = get_configuration("${CONF_DIR}/${input_preset}/configuration", $BUILD_DIR, $verbose);

#for each test
foreach my $test (@$lst_test){
    if($verbose){ print "----------\n"; }
    #generate
    if( $generation ){
        if($verbose){ print "Generating " . $test->getName() . "\n"; }
        if( !generate_test($BUILD_DIR, $BFM_EXE, $overwrite, $temp_dir, $test) ){ next; }
    }
    #execute
    if( $execution ){
        if($verbose){ print "Executing " . $test->getName() . "\n"; }
        if( !execute_test($temp_dir, $test) ){ next; }
    }
    #analyze
    if( $analysis ){
        if($verbose){ print "Analyzing " . $test->getName() . "\n"; }
        if( !analyze_test($temp_dir, $test) ){ next; }
    }
}
if($verbose){ print "----------\n";    }

#print the summary
Test->printSummary($lst_test);

