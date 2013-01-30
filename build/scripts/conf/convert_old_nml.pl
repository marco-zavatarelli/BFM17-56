#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;

my $out = $ARGV[0];
my $in  = $ARGV[1];

open OUT , ">>" , $out or die "cannot open input file $out: $!";
open IN , "<" , $in or die "cannot open input file $in: $!";
my @lines = <IN>;
#print OUT $msg_ini;
foreach my $line (@lines){
    my($filename, $directories, $suffix) = fileparse($in);
    $line =~ s/^\s*\/$/\tfilename_nml_conf = \'$filename\',\n\//;
    print OUT $line;
}
close IN;
close OUT;

