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
my $reg_array='(\w+)(\w{3})\((\d+)\,\:\)';
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
    my @const  = keys %$const_ref;

    foreach my $list (@$lists_ref){
        #check _parameters lists
        if( $list->{NAME} =~ /(.*)_parameters$/ ){
            if ( $VERBOSE ){ print "\tLIST: $list->{NAME}\n"; }
            my $nml_name = $1;
            my $grp_name = "${nml_name}Plankton";
            #check if the group exists in the memory layout for the namelist
            if( exists $$groups_ref{$grp_name} ){
                if ( $VERBOSE ){ print "\t\tFound correspondance with $grp_name\n"; }
                my $list_name = $list->name;
                #check all the parameters which are part of this group
                my @params_grp = ();
                foreach my $param ( sort { $$params_ref{$a}->getIndex() cmp $$params_ref{$b}->getIndex() } keys %$params_ref ){
                    my $prm_grp_name = $$params_ref{$param}->getGroup();
                    if( $prm_grp_name && $prm_grp_name eq $grp_name ){
                        push ( @params_grp, $param );
                    }
                }
                #add new parameter for output comment with parameters in group
                $list->add_elements(\@params_grp);
                my $grp_size = scalar(@params_grp);
                my %pred_terms  = ();
                my %pred_line   = ();
                my %pred_cols   = ();

                #remove external elements in the list or add 0's
                foreach my $element ( @{$list->slots} ){
                    #get the number of values inside the line
                    my $found_group = 0;
                    my $columns = 0;

                    if( $element eq "filename_nml_conf" ){
                        #avoid this element
                        $columns = 1;
                    }elsif( $element =~ /$reg_array/ ){
                        #element type array "nameACRONYM(number,:) = word , word , word"
                        my $prename = "$1$2";
                        my $prename_index = "$1$2$3";
                        my $acro = $2;
                        my $index_num = $3;
                        #print "$prename - $prename_index - $acro - $index_num\n";
                        if ( $index_num > $grp_size ){ 
                            print "WARNING: ($index_num > $grp_size) removing element $element in namelist $nml_name\n";
                            $list->remove($element);  
                        }

                        #search for the group which has the acronym
                        foreach my $group_name (keys %$groups_ref){
                            if( $$groups_ref{$group_name}->getAcro() eq $acro ){ 
                                #calculate the number of elements belong to this group
                                my @params_grp_inside = ();
                                foreach my $param ( sort { $$params_ref{$a}->getIndex() cmp $$params_ref{$b}->getIndex() } keys %$params_ref ){
                                    my $prm_grp_name = $$params_ref{$param}->getGroup();
                                    if( $prm_grp_name && $prm_grp_name eq $group_name ){
                                        push ( @params_grp_inside, $param );
                                        $columns++;
                                    }
                                }
                                #add new parameter for output comment with parameters in subgroups
                                $list->add_subElements($acro, \@params_grp_inside);
                                $found_group = 1;
                                last;
                            }
                        }
                        if( ! $found_group ){ print "WARNING: in param $element not found group for acronym $acro\n"; last; }

                        $pred_terms{$prename_index}  = $columns;
                        $pred_line{$prename}         = exists $pred_line{$prename} ? ($pred_line{$prename}+1) : 1;
                        $pred_cols{$prename}         = $columns;
                        #print "$prename - $index_num - $prename_index\n";

                    }else{
                        #element type normal "name"
                        $columns = $grp_size;
                    }


                    #check the number parameters of this group
                    #should be the number of columns exist in namelist params
                    my $values_num = $#{${$list->hash}{$element}{value}}+1;
                    if( $values_num > $columns ){
                        print "WARNING: ($values_num > $columns) removing element values from $element in namelist $nml_name\n";
                        splice(@{${$list->hash}{$element}{value}} , $columns, ($values_num - $columns) );
                        splice(@{${$list->hash}{$element}{typesv}}, $columns, ($values_num - $columns) );
                    }elsif( $values_num < $columns ){
                        print "WARNING: ($values_num < $columns) adding zero values to element $element in namelist $nml_name\n";
                        my @temp_new_values = map '0.0', 1..($columns - $values_num);
                        my @temp_new_typesv = map ${${$list->hash}{$element}{typesv}}[0], 1..($columns - $values_num);
                        push( @{${$list->hash}{$element}{value}} , @temp_new_values );
                        push( @{${$list->hash}{$element}{typesv}}, @temp_new_typesv );
                    }
                }

                #check the number of arrays inside the namelist
                #should be the number of elements of the group
                foreach my $pred (keys %pred_line){
                    my $values_line = $pred_line{$pred};
                    my $values_cols = $pred_cols{$pred};
                    #print "$values_line < $grp_size @params_grp\n";
                    if( $values_line < $grp_size  ){
                        foreach my $index (1 .. $grp_size ){
                            my $name = "$pred$index";
                            my @temp_new_columns = map '0.0', 1..($values_cols);
                            if ( ! exists $pred_terms{$name} ){
                                print "WARNING: ($values_line < $grp_size) adding new predator term $pred($index,:) in namelist $nml_name\n";
                                my $new_line = "&$list_name $pred($index,:) = @temp_new_columns \/";
                                $list->parse({text => $new_line, merge => 1});
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
            #print Dumper ($nml) , "\n";
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
                    #check if is an array to add the comments
                    if( $line =~ /$reg_array/ ){
                        my $acro  = $2;
                        my $group = ${$nml->elements}[($3-1)];
                        my @line_tmp = ( "!   $group", " " ,@{$nml->subElements($acro)} );
                        calculateMaxLen(\@line_tmp, \@max_len_array);
                        push( @tbl, [@line_tmp] );
                    }


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
