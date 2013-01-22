#!/usr/bin/perl -w

#Author: Esteban Gutierrez esteban.gutierrez@cmcc.it

package read_namelist;

use strict;
use warnings;
use Exporter;
use F90Namelist;

use classes;


########### VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(read_namelist);
########### VARIABLES ##########################


########### FUNCTIONS ##########################

sub read_namelist{
    my $nml_val = shift; #input file
    my $lst_nml = shift; #output for names
    
    my $nl = Fortran::F90Namelist->new() or die "Couldn't get object\n";

    # Read one namelist from file
    #open(NAMELIST , "< t/$nml_val") or die "Couldn't open file: $?\n";
    print "FILE: $nml_val\n";
    # $nl->parse(file     => "$nml_val",
    #            merge    => 0,
    #            all      => 1,);
    open NML_VAL, "<", "$nml_val" or die "$nml_val cannot be opened: $!";
    
    $nl->parse(file => \*NML_VAL);
    print "NAME: " . $nl->name() . " SLOTS: " . $nl->nslots() . "\n";
    print "Content: ", join(",  ", @{$nl->slots}), "\n";

    # foreach my $key ( sort keys %{$nl->hash()} ){
    #     print "$key => ";
    #     print "@{${${$nl->hash()}{$key}}{'value'}}" . "; ";
    #     # foreach my $key1 ( sort keys %{${$nl->hash()}{$key}} ){
    #     #     print $key1 . "=" . ${${$nl->hash()}{$key}}{$key1} . "; ";
    #     # }
    #     print "\n";
    # }

    print "--------------\n";
    my $nl1 = Fortran::F90Namelist->new() or die "Couldn't get object\n";
    $nl1->parse(file => \*NML_VAL);
    print "NAME: " . $nl1->name() . " SLOTS: " . $nl1->nslots() . "\n";
    print "Content ", join(",  ", @{$nl1->slots}), "\n";
}

1;
