#!/usr/bin/perl -w

#Author: Esteban Gutierrez esteban.gutierrez@cmcc.it

package read_namelist;

use strict;
use warnings;
use Exporter;
use F90Namelist;
use Data::Dumper;

use classes;


########### VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(read_namelist);
########### VARIABLES ##########################


########### FUNCTIONS ##########################

sub read_namelist{
    my $nml_val = shift; #input file
    my $lst_nml = shift; #output for names
    
    my $nl = F90Namelist->new(debug => 0) or die "Couldn't get object\n";

    # Read one namelist from file
    open(NAMELIST , "< $nml_val") or die "Couldn't open file: $nml_val. $?\n";
    my @lines = <NAMELIST>;
    close(NAMELIST);

    #process each namelist in the file
    my $block = '';
    foreach my $line (@lines){
        $line =~ s/!.*//;
        if ( $line ){
            if( $line =~ m/^\s*(\&.*)/ ){
                $block = $1;
            }else{
                $block .= $line;
                if( $line =~ m/^\// ){
                    $nl->parse(text => $block);                
                    print "NAME: " . $nl->name() . " SLOTS: " . $nl->nslots() . "\n";
                    print Dumper($nl->hash()) , "\n";
                    $block = '';
                }
            }
        }
    }

    print "F90 format:\n", $nl->output();
    exit;

}

&read_namelist( $ARGV[0] );

1;
