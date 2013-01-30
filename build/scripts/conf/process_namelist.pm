#!/usr/bin/perl -w

#Author: Esteban Gutierrez esteban.gutierrez@cmcc.it

package process_namelist;

use strict;
use warnings;
use Exporter;
use F90Namelist;
use Data::Dumper;

use classes;


########### VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(process_namelist);
########### VARIABLES ##########################


########### FUNCTIONS ##########################

sub process_namelist{
    my $nml_val = shift; #input file
    my $lst_nml = shift; #output for names
    
    my $index = 0;

    # Read one namelist from file
    open(NAMELIST , "< $nml_val") or die "Couldn't open file: $nml_val. $?\n";
    my @lines = <NAMELIST>;
    close(NAMELIST);

    #process each namelist in the file
    my $block     = '';
    my $lines_not = '';
    foreach my $line (@lines){
        $line =~ s/!.*//;
        if ( $line ){
            if( $line =~ m/^\s*(\&.*)/ ){
                $block = $1;
            }else{
                if( $block ){ 
                    $block .= $line;
                }else{
                    #lines which are not processed
                    $lines_not .= $line;
                }

                #is the end of the namelist?
                if( $line =~ m/^\s*\// ){
                    $$lst_nml[$index] = F90Namelist->new(debug => 0) or die "Couldn't get object\n";
                    $$lst_nml[$index]->parse(text => $block);
                    #print "OUTPUT:\n", $$lst_nml[$index]->output();
                    $index++;
                    $block = '';
                }
            }
        }
    }
    return $lines_not;
}

1;
