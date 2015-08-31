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
my $reg_comm='\s*\!\s*(.+)'; 
my $MAX_PARAMS_PER_LINE = 9;
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

sub getGroupComponents{
    my ( $groups_ref, $params_ref, $group_acro,  ) = @_;
    my @components = ();

    foreach my $group_name (keys %$groups_ref){
        if( $$groups_ref{$group_name}->getAcro() eq $group_acro ){ 
            #calculate the number of elements belong to this group
            foreach my $param ( sort { $$params_ref{$a}->getIndex() cmp $$params_ref{$b}->getIndex() } keys %$params_ref ){
                my $prm_grp_name = $$params_ref{$param}->getGroup();
                if( $prm_grp_name && $prm_grp_name eq $group_name ){
                  push( @components, $param );  
                }
            }
        }
    }

    return \@components;
}


sub check_cols{
    #check the number parameters of this group
    #should be the number of columns exist in namelist params
    my ($element, $list, $group_compo, $comment_compo, $values ) = @_;

    my $nml_name   = $list->name;
    my $columns    = scalar(@$group_compo);
    my $values_num = scalar(@{$list->get_values($element)});

    if( $values_num > $columns ){
        print "WARNING: ($values_num > $columns) removing element values from $element in namelist $nml_name\n";
        $list->remove_values($element, $columns, ($values_num - $columns));
    }elsif( $values_num < $columns ){
        print "WARNING: ($values_num < $columns) adding zero values to element $element in namelist $nml_name\n";

        my $index = 0;
        my @final_array = ();
        my %lst_values  = ();
        @{\%lst_values}{@$group_compo} = map '0.0', 1..$columns; # create the list of components with 0
        foreach my $value (@$comment_compo){
            $lst_values{$value} = @{$values}[$index] if ( exists $lst_values{$value} );
            $index++;
        }
        foreach my $key (@$group_compo){
            push( @final_array, $lst_values{$key} );            
        }

        # print "$element - $nml_name:\n";
        # print "\tFINAL: @final_array\n";
        # print "\tINITIAL: @$values\n";
        # print "\tGROUP: @$group_compo\n";
        # print "\tINITAL: @$comment_compo\n";
        $list->change_values($element, \@final_array);
    }
}


sub process_namelist{
    my ($nml_val, $lst_nml, $lst_com, $debug ) = @_;
    # nml_val - input file
    # lst_nml - output for names
    # lst_com - output for commentaries
    # debug   - print debug messages
    
    my $index = 0;

    # Read one namelist from file
    open(NAMELIST , "< $nml_val") or die "Couldn't open file: $nml_val. $?\n";
    my @lines = <NAMELIST>;
    close(NAMELIST);

    #process each namelist in the file
    my $block          = '';
    my $comm           = '';
    my $sub_comm       = '';
    my %sub_comm_param = ();

    foreach my $line (@lines){
        if( $line =~ m/^\!\s*NAMELIST (.*)/ ){
            #start namelist comments
            $comm = "! $1 ";
        }
        elsif( $line =~ m/^\s*(\&.*)/ ){
            #start the namelist values
            $block = $1;
            $$lst_com[$index] = $comm;
            $comm  = '';
        }else{
            if( $block ){ 
                if ( $line =~ m/$reg_comm/ ){ 
                    #found comment inside namelist (should be added to the next param
                    $sub_comm = $1 if ( ! $sub_comm );
                }else{
                    if ( $sub_comm ){
                        #found param for comment inside namelist
                        my $param_temp = ($line =~ m/\s*([^\s\=]+)/)[0];
                        my @comm_temp = split('\s+',$sub_comm);
                        $sub_comm_param{$param_temp} = \@comm_temp;
                        $sub_comm='';
                        if( $debug ){ print "$param_temp -> @comm_temp\n"; }
                    }
                }
                $block .= $line;
            }

            if( $comm ){
                #it is a namelist comment (before starts namelists values)
                $comm .= $line;
            }

            #is the end of the namelist?
            if( $line =~ m/^\s*\// ){
                if( $debug ){ print "BLOCK:\n$block\n"; }
                $$lst_nml[$index] = F90Namelist->new(debug => 0) or die "Couldn't get object\n";
                $$lst_nml[$index]->parse(text => $block);
                #add commentaries of parameters
                foreach my $key (keys %sub_comm_param){ 
                    if( $$lst_nml[$index]->get_parameter($key) ){
                        $$lst_nml[$index]->add_subComments($key, $sub_comm_param{$key});
                    }else{
                        my $comm_ref = $sub_comm_param{$key};
                        print "WARNING: problems adding comment \"@{$comm_ref}\" for param \"$key\" in namelist \"" . $$lst_nml[$index]->name . "\"\n";
                    }
                }
                if( $debug ){ print "OUTPUT:\n", $$lst_nml[$index]->output(); }
                $index++;
                $block = '';
                $comm = '';
                $sub_comm = '';
                %sub_comm_param = ();
            }
        }
    }
}


sub check_namelists{
    my ($lists_ref, $groups_ref, $params_ref, $const_ref, $index_ref, $debug ) = @_;
    my @const  = keys %$const_ref;

    foreach my $list (@$lists_ref){
        #check _parameters lists
        if( $list->name =~ /(.*)_parameters$/ ){
            if ( $debug ){ print "\tLIST: " . $list->name . "\n"; }
            my $nml_name = $1;
            my $grp_name = "${nml_name}Plankton";
            my $check = 0;
            #check if the group exists in the memory layout for the namelist
            if( exists $$groups_ref{"${nml_name}Plankton"}){
                $grp_name = "${nml_name}Plankton";
                $check = 1;
            }elsif( exists $$groups_ref{"${nml_name}"}){
                $grp_name = "${nml_name}";
                $check = 1;
            }

            if( $check ){
                if ( $debug ){ print "\t\tFound correspondance with $grp_name\n"; }
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
                my %pred_terms   = ();
                my %pred_line    = ();
                my %pred_cols    = ();
                my %pred_comm    = ();
                my $last_comment = ''; 

                #check every element of the namelist
                foreach my $element ( @{$list->get_all_parameters()} ){
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
                        
                        #search for the components which belongs to the group
                        my $params_grp_inside = getGroupComponents($groups_ref, $params_ref, $acro);
                        if ( scalar(@$params_grp_inside) > 0 ){
                            $columns = scalar(@$params_grp_inside);
                            #add new parameter for output comment with parameters in subgroups
                            $list->add_subElements($acro, $params_grp_inside);
                            $pred_terms{$prename_index}  = $columns;
                            $pred_line{$prename}         = exists $pred_line{$prename} ? ($pred_line{$prename}+1) : 1;
                            $pred_cols{$prename}         = $columns;
                            $pred_comm{$prename}         = $list->subComments($element) if ( ! exists $pred_comm{$prename} );

                            #check if the comment of the group exists
                            if( ! $pred_comm{$prename} ){ print "ERROR: Commentary is missing in $element in namelist $nml_name\n"; exit 1; }


                            check_cols($element, $list, $params_grp_inside, $pred_comm{$prename}, $list->get_values($element) );
                        }else{
                            print "WARNING: in param $element not found group for acronym $acro\n"; last; 
                        }
                    }else{
                        #element type normal "name"
                        $last_comment = $list->subComments($element) if $list->subComments($element);
                        #check if the comment exists
                        if( ! $last_comment ){ print "ERROR: Commentary is missing in $element in namelist $nml_name\n"; exit 1; }

                        $list->subComments($element);
                        $columns = $grp_size;
                        check_cols($element, $list, \@params_grp, $last_comment, $list->get_values($element) );
                    }
                }

                #check the number of arrays inside the namelist
                #should be the number of elements of the subgroup
                foreach my $pred (keys %pred_line){
                    my $values_line = $pred_line{$pred};
                    my $values_cols = $pred_cols{$pred};
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
        elsif( $list->{NAME} eq "bfm_save_nml" ){
            foreach my $element ( @{$list->slots} ){
                if( $element eq "ave_save" ){
                    foreach my $value ( @{$list->get_values($element)} ){
                        my $tmp = $value;
                        if ( $tmp =~ /(.*)\(ii(.*)\)/ ){
                            if( ! exists $$params_ref{$2} ){ print "WARNING: output $tmp does not exists\n"; }
                            $tmp = $1;
                        }

                        if ( $tmp =~ /(jbot|jsur|jriv)(.*)/ ){ $tmp = $2; }

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

        #check bfm_init list
        elsif( $list->{NAME} eq "bfm_init_nml" ){
            foreach my $element ( @{$list->slots} ){
                if( $element =~ m/InitVar\((\w+)\)\%(\w+)/ ){
                    my $ele_name = $1;
                    my $ele_end  = $2;
                    if( exists $$index_ref{$ele_name} ){ 
                        my $element_new="InitVar($$index_ref{$ele_name})%${ele_end}";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
            }
        }
        #check namtrc_dta list
        elsif( $list->{NAME} eq "namtrc_dta" ){
            foreach my $element ( @{$list->slots} ){
                if( $element =~ m/sn_trcdta\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="sn_trcdta($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
                elsif( $element =~ m/rn_trfac\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="rn_trfac($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
            }
        }
        #check namtrc_bc list
        elsif( $list->{NAME} eq "namtrc_bc" ){
            foreach my $element ( @{$list->slots} ){
                if( $element =~ m/sn_trcsbc\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="sn_trcsbc($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
                elsif( $element =~ m/rn_trsfac\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="rn_trsfac($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
                elsif( $element =~ m/sn_trccbc\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="sn_trccbc($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
                elsif( $element =~ m/rn_trcfac\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="rn_trcfac($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
                elsif( $element =~ m/sn_trcobc\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="sn_trcobc($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
                elsif( $element =~ m/rn_trofac\((\w+)\)/ ){
                    my $ele_name = $1;
                    if( exists $$index_ref{$ele_name} ){
                        my $element_new="rn_trofac($$index_ref{$ele_name})";
                        $list->change_name($element, $element_new);
                    }else{ print "WARNING: parameter \"$ele_name\" in \"$element\" does not exists\n";  }
                }
            }
        }
    }
}


sub print_namelists{
    my ( $lst_nml, $lst_com, $lst_index, $out_dir, $debug ) = @_ ;

    my $index = 0;
    foreach my $nml (@$lst_nml){
        if ( $nml->get_parameter('filename_nml_conf') ){
            #print Dumper ($nml) , "\n";
            #insert all elements in a table
            my $nml_name = "$out_dir/" . ${$nml->get_values('filename_nml_conf')}[0];
            $nml->remove('filename_nml_conf');
            my @max_len_array = ();
            my @tbl = ();         
            my $pred_terms = '';

            #first get column sizes to print with a beauty format
            foreach my $line ( split(/\n/,$nml->output) ){
                if( $line =~ "^&.*" ){
                    #print header
                    my @line_tmp = ( $line );
                    calculateMaxLen(\@line_tmp, \@max_len_array);
                    push( @tbl, [@line_tmp] );
                    if( $nml->elements() ){
                        my @line_tmp = ( "!", " " ,@{$nml->elements} );
                         calculateMaxLen(\@line_tmp, \@max_len_array);
                         push( @tbl, [@line_tmp] );
                    }
                }elsif( $line =~ "^\/.*" ){
                    #print footer
                    my @line_tmp = ( $line );
                    calculateMaxLen(\@line_tmp, \@max_len_array);
                    push( @tbl, [@line_tmp] );

                }else{
                    #check if is an array to add the comments
                    if( $line =~ /$reg_array/ ){
                        my $acro  = $2;
                        if( defined $nml->elements ){ 
                            my $group = ${$nml->elements}[($3-1)];
                            if ( $pred_terms ne $acro ){
                                $pred_terms = $acro;
                                my @line_tmp = ( "! ", " " ,@{$nml->subElements($acro)} );
                                calculateMaxLen(\@line_tmp, \@max_len_array);
                                push( @tbl, [@line_tmp] )
                            };
                            push( @tbl, [( "!   $group" )] );
                        }else{
                            print "WARNING: Group $2 not defined\n";
                        }
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
            if( $nml->name eq 'bfm_init_nml' ){
                #print the list of param indix
                print NML_OUT "! ";
                my $params_per_line = 0;
                print NML_OUT "Index of parameters for using inside InitVar structure\n!   ";
                foreach my $par ( sort { $$lst_index{$a} <=> $$lst_index{$b} } keys %$lst_index ){
                    if ( $params_per_line == $MAX_PARAMS_PER_LINE ){
                        print NML_OUT "\n!   ";
                        $params_per_line = 0;
                    }
                    print NML_OUT "$par=$$lst_index{$par}, ";
                    $params_per_line++;
                }
                print NML_OUT "\n";
            }
            my @pad_len = map { "%-${_}s  " } @max_len_array;
            foreach my $line (@tbl){
                foreach my $idx ( 0..$#{$line} ){
                    printf NML_OUT "$pad_len[$idx]", ${$line}[$idx];
                }
                printf NML_OUT "\n" ;
            }
            printf NML_OUT "\n\n\n" ;
            close NML_OUT;
            if( $debug ){ print "---------- $nml_name\n"; }
        }
    }
}

1;
