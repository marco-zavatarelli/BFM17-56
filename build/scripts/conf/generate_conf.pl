#!/usr/bin/perl -w

#Author: Esteban Gutierrez esteban.gutierrez@cmcc.it

use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;

use process_memLayout; #reomve trailing white spaces
use process_namelist; # get the namelist values in the hash
use F90Namelist; #get the namelists
use print_f90; # write the variables in the output file
use classes;


my ($input, $proto_dir, $out_dir, @cpp_defs);

#fix values
my $VERBOSE = 0;
my $HELP    = 0;
my @PROTOS_NAME = qw(ModuleMem AllocateMem set_var_info_bfm init_var_bfm INCLUDE);
my @PROTOS_EXT  = qw(F90       F90         F90              F90          h      );

#structures to allocate the parameters read in memmory layour
my @lst_nml   = ();
my %lst_group = ();
my %lst_param = ();
my %lst_sta   = ();
my %lst_const = ();

sub usage(){
    print "usage: $0 {-D[cpp_def] -r [mem_layout] -n [namelist] -f [prototype_dir] -t [output_dir]} [-v]\n\n";
    print "This script generate .F90, .h and .nml files using templates based on configuration files\n\n";
    print "MUST specify at least one these OPTIONS:\n";
    print "\t-D[cpp_def]          defines\n";
    print "\t-r [mem_layout]      memory layout (with or without namelists inside)\n";
    print "\t-f [prototype_dir]   input dir for prototype files\n";
    print "\t-t [output_dir]      output dir for generated files\n";
    print "alternative OPTIONS are:\n";
    print "\t-v                   verbose mode\n";
}

use Getopt::Long qw(:config bundling noignorecase); # for getopts compat
GetOptions(
    'r=s'  => \$input,
    'f=s'  => \$proto_dir,
    't=s'  => \$out_dir,  
    'D=s@' => \@cpp_defs,
    'v'    => \$VERBOSE,
    'h'    => \$HELP,
    ) or &usage() && exit;
if ( $HELP ){ &usage(); exit; }
if ( !$input || !$proto_dir || !$out_dir || !@cpp_defs ){ &usage(); exit; }

my $nml_def_tmp = "$out_dir/nml_def_tmp";

#process namelists removing them from the input file
if( $VERBOSE ){ print "Reading namelists...\n"; }
open NML_DEF_TMP, ">", "$nml_def_tmp" or die "$nml_def_tmp cannot be opened: $!";
my $lines_def = process_namelist($input, \@lst_nml);
if( $VERBOSE ){ foreach my $nml (@lst_nml){ print $nml->output; } }
print NML_DEF_TMP $lines_def;
close(NML_DEF_TMP);

#read memory layout file
if( $VERBOSE ){ print "Reading memory layout...\n"; }
process_memLayout("$nml_def_tmp", \%lst_group, \%lst_param, \%lst_sta, \%lst_const, join(' ',@cpp_defs), $VERBOSE );

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
foreach my $nml (@lst_nml){
    if ( $nml->hash()->{'filename_nml_conf'} ){
        my $nml_name = "$out_dir/" . $nml->hash()->{'filename_nml_conf'}->{'value'}[0];
        $nml->remove('filename_nml_conf');
        open  NML_OUT, ">>", "$nml_name" or die "$nml_name cannot be opened: $!";
        print NML_OUT $nml->output;
        close NML_OUT;
        if( $VERBOSE ){ print "---------- $nml_name\n"; }
    }
}



if( $VERBOSE ){ print "Configuration files generation finished\n"; }
