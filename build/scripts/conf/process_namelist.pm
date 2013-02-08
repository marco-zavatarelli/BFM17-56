# DESCRIPTION
#   Process namelist files (read, check and generation)
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

package process_namelist;

use strict;
use warnings;
use Exporter;
use F90Namelist;
use Data::Dumper;

use classes;


########### VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(process_namelist check_namelists print_namelists);
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


sub check_namelists{
    my ($lists_ref, $groups_ref, $params_ref, $const_ref, $VERBOSE ) = @_;
    my %lookup = map {(lc $_, $$groups_ref{$_})} keys %$groups_ref; #lowecase all the group names
    my @const  = keys %$const_ref;

    foreach my $list (@$lists_ref){
        #check _parameters lists
        if( $list->{NAME} =~ /(.*)_parameters$/ ){
            my $nml_name = $1;
            my $grp_name = "${nml_name}plankton";
            #check if the group exists in the memory layout for the namelist
            if( exists $lookup{$grp_name} ){
                if ( $VERBOSE ){ print "\t$nml_name\n"; }
                #check all the parameters which are part of this group
                my $clm_num = 0;
                my @params_grp = ();
                foreach my $param (sort keys %$params_ref){
                    my $prm_grp_name = $$params_ref{$param}->getGroup();
                    if( $prm_grp_name && lc($prm_grp_name) eq $grp_name ){
                        if( $VERBOSE ){ print "\t$param -> $grp_name\n"; }
                        push ( @params_grp, $param );
                        #check the number parameters of this group
                        #which will be the number of columns should exist in namelist params
                        $clm_num++;
                    }
                }

                if( $clm_num > 0 ){
                    #add new parameter for output comment with parameters in group
                    $list->add_elements(\@params_grp);
                    #remove external elements in the list or add 0's
                    foreach my $element ( @{$list->slots} ){
                        if( $VERBOSE ){ print "\t\t$element\n"; }
                        if( $element eq "filename_nml_conf" ){
                            #avoid this element
                        }elsif( $element =~ /\w+\((\d+)\,\:\)/ ){
                            #element type array "name(number,:)"
                            my $values_num = $1;
                            if ( $values_num > $clm_num ){ 
                                print "WARNING: ($values_num > $clm_num) removing element $element in namelist $nml_name\n";
                                $list->remove($element);  
                            }
                        }else{
                            #element type normal "name"
                            my $values_num = $#{${$list->hash}{$element}{value}}+1;
                            if( $values_num > $clm_num ){
                                print "WARNING: ($values_num > $clm_num) removing element values from $element in namelist $nml_name\n";
                                splice(@{${$list->hash}{$element}{value}} , $clm_num, ($values_num - $clm_num) );
                                splice(@{${$list->hash}{$element}{typesv}}, $clm_num, ($values_num - $clm_num) );
                            }elsif( $values_num < $clm_num ){
                                print "WARNING: ($values_num < $clm_num) adding zero values to element $element in namelist $nml_name\n";
                                push( @{${$list->hash}{$element}{value}} , '0.0' x ($clm_num - $values_num) );
                                push( @{${$list->hash}{$element}{typesv}}, (${${$list->hash}{$element}{typesv}}[0]) x ($values_num - $clm_num) );
                            }
                        }
                    }
                }
            }
        }

        #check bfm_save list
        if( $list->{NAME} eq "bfm_save_nml" ){
            foreach my $element ( @{$list->slots} ){
                if( $element eq "ave_save" ){
                    foreach my $value ( @{${$list->hash}{$element}{value}} ){
                        if ( $value =~ /(.*)\(ii(.*)\)/ ){
                            if( ! exists $$params_ref{$2} ){ print "WARNING: output $value does not exists\n"; }
                            $value = $1;
                        }

                        if( ! exists $$params_ref{$value} ){
                            #check if it is part of constituent
                            my @parts = ( $value =~ /(.*)(\D)/ );
                            if( ! exists $$params_ref{$parts[0]} || ! exists ${$$params_ref{$parts[0]}->getComponents()}{$parts[1]} ){
                                print "WARNING: output $value does not exists\n"; 
                            }
                        }
                    }
                }
            }
        }
    }
}


sub print_namelists{
    my ( $lst_nml, $out_dir, $VERBOSE ) = @_ ;

    foreach my $nml (@$lst_nml){
        if ( $nml->hash()->{'filename_nml_conf'} ){
            my $nml_name = "$out_dir/" . $nml->hash()->{'filename_nml_conf'}->{'value'}[0];
            $nml->remove('filename_nml_conf');
            open  NML_OUT, ">>", "$nml_name" or die "$nml_name cannot be opened: $!";

            if( $nml->elements() ){
                my @line_tmp = ( "   ", @{$nml->elements} );
                print NML_OUT "!" . pack('(A20)*', @line_tmp ) . "\n";
            }
             
            foreach my $line ( split(/\n/,$nml->output) ){
                if( $line =~ "^[&\/].*" ){
                    #print header or footer 
                    print NML_OUT $line . "\n"
                }else{
                    my @parts = ( $line =~ /^\s*(.*\=)(.*)/ );
                    my @line_tmp = ( "    $parts[0]" , split( ',', $parts[1]) );
                    print NML_OUT pack( '(A20)*', @line_tmp ) . "\n";
               }
            }
            
            close NML_OUT;
            if( $VERBOSE ){ print "---------- $nml_name\n"; }
        }
    }
}

1;
