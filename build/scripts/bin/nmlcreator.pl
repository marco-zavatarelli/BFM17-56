#!/usr/bin/perl -w

# DESCRIPTION
#   Apped to new configuration file the content of an old namelist file
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
use File::Basename;
use Getopt::Std;

my $HELP = 0;
my ($in_nml, $out_nml);


sub usage(){
    print "usage $0 -o [output_namelist] -i [input_namelist]\n\n";
    print "This script appends to new configuration file the content of an old namelist file\n";
    print "It adds the \"filename_nml_conf\" key with the name of the file at the end of each namelist read\n\n";
    print "\t-o [namelist]      output file where to append at the end the input namelist (complete PATH)\n";
    print "\t-i [namelist]      input file containing namelist/s to read (complete PATH)\n";
}


use Getopt::Long qw(:config bundling noignorecase); # for getopts compat
GetOptions(
    'i=s'  => \$in_nml,
    'o=s'  => \$out_nml,
    'h'    => \$HELP,
    ) or &usage() && exit;
if ( $HELP ){ &usage(); exit; }
if ( !$in_nml || !$out_nml ){ &usage(); exit; }


my $out = $ARGV[0];
my $in  = $ARGV[1];

open OUT , ">>" , $out_nml or die "cannot open output namelist file $out_nml: $!";
open IN , "<" , $in_nml or die "cannot open input file $in_nml: $!";
my @lines = <IN>;

foreach my $line (@lines){
    my($filename, $directories, $suffix) = fileparse($in_nml);
    $line =~ s/^\s*\/$/\tfilename_nml_conf = \'$filename\',\n\//;
    print OUT $line;
}

print "Done!\n";

close IN;
close OUT;
