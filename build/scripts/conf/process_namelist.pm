#!/usr/bin/perl -w

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

sub calculateMaxLen{
    my ( $ref_line, $ref_max_array ) = @_;

    foreach my $idx (0..$#{$ref_line}){
        if ( ! $$ref_max_array[$idx] || length($$ref_line[$idx]) > $$ref_max_array[$idx] ){
            $$ref_max_array[$idx] = length($$ref_line[$idx])
        }
    }
}


sub process_namelist{
    my $nml_val = shift; #input file
    my $lst_nml = shift; #output for names
    my $lst_com = shift; #output for commentaries
    
    my $index = 0;

    # Read one namelist from file
    open(NAMELIST , "< $nml_val") or die "Couldn't open file: $nml_val. $?\n";
    my @lines = <NAMELIST>;
    close(NAMELIST);

    #process each namelist in the file
    my $block     = '';
    my $comm      = '';

    foreach my $line (@lines){
        if( $line =~ m/^\!\s*NAMELIST (.*)/ ){
            $comm = "! $1 ";
        }
        elsif( $line =~ m/^\s*(\&.*)/ ){
            $block = $1;
            $$lst_com[$index] = $comm;
            $comm  = '';
        }else{
            if( $block ){ 
                $block .= $line;
            }

            if( $comm ){
                $comm .= $line;
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
                                my @temp_new_values = map '0.0', 1..($clm_num - $values_num);
                                my @temp_new_typesv = map ${${$list->hash}{$element}{typesv}}[0], 1..($clm_num - $values_num);
                                push( @{${$list->hash}{$element}{value}} , @temp_new_values );
                                push( @{${$list->hash}{$element}{typesv}}, @temp_new_typesv );

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
                        my $tmp = $value;
                        if ( $tmp =~ /(.*)\(ii(.*)\)/ ){
                            if( ! exists $$params_ref{$2} ){ print "WARNING: output $tmp does not exists\n"; }
                            $tmp = $1;
                        }

                        if( ! exists $$params_ref{$tmp} ){
                            #check if it is part of constituent
                            my @parts = ( $tmp =~ /(.*)(\D)/ );
                            if( ! exists $$params_ref{$parts[0]} || ! exists ${$$params_ref{$parts[0]}->getComponents()}{$parts[1]} ){
                                print "WARNING: output $tmp does not exists\n"; 
                            }
                        }
                    }
                }
            }
        }
    }
}


sub print_namelists{
    my ( $lst_nml, $lst_com, $out_dir, $VERBOSE ) = @_ ;

    my $index = 0;
    foreach my $nml (@$lst_nml){
        if ( $nml->hash()->{'filename_nml_conf'} ){
            #print Dumper ($lst_nml) , "\n";
            #first get column sizes to print with a beauty format
            #insert all elements in a table
            my $nml_name = "$out_dir/" . $nml->hash()->{'filename_nml_conf'}->{'value'}[0];
            $nml->remove('filename_nml_conf');
            my @max_len_array = ();
            my @tbl = ();

            if( $nml->elements() ){
                my @line_tmp = ( "!", " " ,@{$nml->elements} );
                calculateMaxLen(\@line_tmp, \@max_len_array);
                push( @tbl, [@line_tmp] );
            }
             
            foreach my $line ( split(/\n/,$nml->output) ){
                if( $line =~ "^[&\/].*" ){
                    #print header or footer 
                    my @line_tmp = ( $line );
                    calculateMaxLen(\@line_tmp, \@max_len_array);
                    push( @tbl, [@line_tmp] );
                }else{
                    my @parts = ( $line =~ /^\s*(.*)\=(.*)/ );
                    my @line_tmp = ();
                    push( @line_tmp, "    $parts[0]", "=", split( ',', $parts[1]) );
                    calculateMaxLen(\@line_tmp, \@max_len_array);
                    push( @tbl, [@line_tmp] );
               }
            }

            #print the formated output to the file
            open  NML_OUT, ">>", "$nml_name" or die "$nml_name cannot be opened: $!";
            print NML_OUT $$lst_com[$index++];
            my @pad_len = map { "%-${_}s  " } @max_len_array;
            foreach my $line (@tbl){
                foreach my $idx ( 0..$#{$line} ){
                    printf NML_OUT "$pad_len[$idx]", ${$line}[$idx];
                }
                printf NML_OUT "\n" ;
            }
            printf NML_OUT "\n\n\n" ;
            close NML_OUT;
            if( $VERBOSE ){ print "---------- $nml_name\n"; }
        }
    }
}

1;
