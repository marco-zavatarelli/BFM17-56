#!/usr/bin/perl -w

# DESCRIPTION
#   Generate .h .f90 and namelist files
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

use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

use process_memLayout; #reomve trailing white spaces
use process_namelist; # get the namelist values in the hash
use F90Namelist; #get the namelists
use print_f90; # write the variables in the output file
use classes;


my ($input_mem, $input_nml,, $proto_dir, $out_dir, @cpp_defs);

#fix values
my $VERBOSE = 0;
my $HELP    = 0;
my @PROTOS_NAME = qw(ModuleMem AllocateMem set_var_info_bfm init_var_bfm INCLUDE);
my @PROTOS_EXT  = qw(F90       F90         F90              F90          h      );

#structures to allocate the parameters read in memmory layour
my @lst_nml   = ();
my @lst_com   = ();
my %lst_group = ();
my %lst_param = ();
my %lst_sta   = ();
my %lst_const = ();

sub usage(){
    print "usage: $0 {-D[cpp_def] -r [mem_layout] -n [namelist] -f [prototype_dir] -t [output_dir]} [-v]\n\n";
    print "This script generate .F90, .h and .nml files using templates based on configuration files\n\n";
    print "MUST specify at least one these OPTIONS:\n";
    print "\t-D[cpp_def]          defines\n";
    print "\t-r [mem_layout]      memory layout\n";
    print "\t-n [namelist]        file containing namelists\n";
    print "\t-f [prototype_dir]   input dir for prototype files\n";
    print "\t-t [output_dir]      output dir for generated files\n";
    print "alternative OPTIONS are:\n";
    print "\t-v                   verbose mode\n";
}

use Getopt::Long qw(:config bundling noignorecase); # for getopts compat
GetOptions(
    'r=s'  => \$input_mem,
    'n=s'  => \$input_nml,
    'f=s'  => \$proto_dir,
    't=s'  => \$out_dir,  
    'D=s@' => \@cpp_defs,
    'v'    => \$VERBOSE,
    'h'    => \$HELP,
    ) or &usage() && exit;
if ( $HELP ){ &usage(); exit; }
if ( !$input_mem || !$input_nml || !$proto_dir || !$out_dir || !@cpp_defs ){ &usage(); exit; }

#process namelists removing them from the input file
if( $VERBOSE ){ print "Reading namelists...\n"; }
process_namelist($input_nml, \@lst_nml, \@lst_com);
if( $VERBOSE ){ foreach my $nml (@lst_nml){ print $nml->output; } }

#read memory layout file
if( $VERBOSE ){ print "Reading memory layout...\n"; }
process_memLayout($input_mem, \%lst_group, \%lst_param, \%lst_sta, \%lst_const, join(' ',@cpp_defs), $VERBOSE );

#check consistency between namelists and memory_layout
if( $VERBOSE ){ print "Checking namelist and memory layout...\n"; }
check_namelists( \@lst_nml, \%lst_group, \%lst_param, \%lst_const, $VERBOSE );

#write memory output files
if( $VERBOSE ){ print "Printing memory layout...\n"; }
foreach my $idx ( 0 .. $#PROTOS_NAME ){
    my $name = $PROTOS_NAME[$idx];
    my $ext  = $PROTOS_EXT[$idx];
    if( $VERBOSE ){ print "---------- ${name} ini \n"; }
    print_f90("$proto_dir/${name}.proto", "$out_dir/${name}.${ext}", \%lst_group, \%lst_param, \%lst_sta, \%lst_const, $VERBOSE);
    if( $VERBOSE ){ print "---------- ${name} end \n\n"; }
}

#write namelists output files
if( $VERBOSE ){ print "Writing namelists...\n"; }
print_namelists( \@lst_nml, \@lst_com, $out_dir, $VERBOSE );


#foreach my $key (keys %lst_group){ $lst_group{$key}->print(); }
#foreach my $key (keys %lst_param){ $lst_param{$key}->print(); }


if( $VERBOSE ){ print "Configuration files generation finished\n"; }
