#!/usr/bin/perl -w

# DESCRIPTION
#   Generate namelist files
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

use process_namelist; # get the namelist values in the hash
use F90Namelist; #get the namelists
use print_f90; # write the variables in the output file
use classes;


my ($in_dir, $input_nmls, $out_dir);


#fix values
my $VERBOSE = 0;
my $DEBUG   = 0;
my $HELP    = 0;

sub usage(){
    print "usage: $0 {-n [namelist]-t [output_dir]} [-v]\n\n";
    print "This script generate .nml files using lists of namelist in one file\n\n";
    print "MUST specify these OPTIONS:\n";
    print "\t-i [in_dir]          input directory where the namelist list files are\n";
    print "\t-n [namelist]        namelist files (separated by spaces)\n";
    print "\t-o [output_dir]      output directory for generated files\n";
    print "alternative OPTIONS are:\n";
    print "\t-v                   verbose mode\n";
    print "\t-d                   debug mode\n";
}

use Getopt::Long qw(:config bundling noignorecase); # for getopts compat
GetOptions(
    'i=s'  => \$in_dir,  
    'n=s'  => \$input_nmls,
    'o=s'  => \$out_dir,  
    'v'    => \$VERBOSE,
    'd'    => \$DEBUG,
    'h'    => \$HELP,
    ) or &usage() && exit;
if ( $HELP ){ &usage(); exit; }
if ( !$input_nmls || !$out_dir || !$in_dir ){ &usage(); exit; }

if( $VERBOSE ){ print "Processing: $input_nmls...\n"; }
my @lists = split(' ', $input_nmls);
foreach my $nml (@lists){
    my @lst_nml   = ();
    my @lst_com   = ();
    my %lst_index = ();

    #process namelists removing them from the input file
    if( $VERBOSE ){ print "Reading namelist $nml...\n"; }
    process_namelist("$in_dir/$nml", \@lst_nml, \@lst_com, $DEBUG);
    
    #write namelists output files
    if( $VERBOSE ){ print "Writing namelist $nml...\n"; }
    print_namelists( \@lst_nml, \@lst_com, \%lst_index, $out_dir, $DEBUG );
}

if( $VERBOSE ){ print "Print namelists finished\n"; }
