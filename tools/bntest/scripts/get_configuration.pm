#!/usr/bin/perl -w

# DESCRIPTION
#   Extract anc check configuration file
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

package get_configuration;

use 5.008_002;

use strict;
use warnings;
use Exporter;
use Data::Dumper;

use classes;

########### VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(get_configuration);
########### VARIABLES ##########################


my $VERBOSE = 0;

sub get_configuration{
    my ($conf_file, $verbose) = @_;
    my %user_conf;
    $VERBOSE = $verbose;

    open CONFIG, "<$conf_file" or die "Could not open configuration file for ${conf_file}: $!";
    while (<CONFIG>) {
        chomp;                  # no newline
        s/#.*//;                # no comments
        s/^\s+//;               # no leading white
        s/\s+$//;               # no trailing white
        next unless length;     # anything left?
        my ($var, $value) = split(/\s*=\s*/, $_, 2);
        my @comasep = split(/\s*,\s*/, "$value");
        $user_conf{$var} = \@comasep;
    }
    close(CONFIG);
    check_conf(\%user_conf);
    return fill_conf(\%user_conf);
}

sub check_conf{
    my ($user_conf) = @_;

    #if($VERBOSE){ print Dumper($user_conf) . "\n"; }
    if( ! exists $$user_conf{TESTS} ){ 
        print "ERROR: Configuration file must contain \"TESTS\" variable\n"; 
        exit; 
    }
    my $num_tests = scalar(@{$$user_conf{TESTS}});
    if( exists $$user_conf{FORCINGS} && scalar(@{$$user_conf{FORCINGS}}) != $num_tests ){ 
        print "ERROR: \"FORCINGS\" var must have same number of elements as \"TESTS\" var ($num_tests elements)\n"; 
        exit; 
    }
}

sub fill_conf{
    my($user_conf) = @_;
    my @lst_test;
    my $name; my $forcing='';

    my $num_tests = $#{$$user_conf{TESTS}};
    foreach my $var (0..$num_tests){
        $name = ${$$user_conf{TESTS}}[$var];
        if( exists $$user_conf{FORCINGS} ){ $forcing = ${$$user_conf{FORCINGS}}[$var] }
        push(@lst_test, new Test($name, $forcing) );
    }
    if($VERBOSE){ Test->printAll(\@lst_test); }
 
    return \@lst_test;
}

1;
