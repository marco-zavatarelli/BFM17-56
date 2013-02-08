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

