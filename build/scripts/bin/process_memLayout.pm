#!/usr/bin/perl -w

# DESCRIPTION
#   Generate the memory layout from configuration file using templates
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

package process_memLayout;

use strict;
use warnings;
use Exporter;
use Data::Dumper;

use classes;


########### REGULAR EXPRESSIONS ##########################
my $XPR_GLOBAL_COMMENT = '([^#]*)#{0,1}(.*)'; # dont process commentaries

#my $XPR_START_BLOCK = '^([123])d-([^\s\n]+)(?:\s+-if-exist\s+){0,1}([^-\n]*)(?:-Z\s+){0,1}(.*)'; #block to indicate dimension type and other characteristics of the parameter
my $XPR_START_BLOCK = '^([123])d-([^\s\n-]+)(?:-(pel|ben|ice)){0,1}(?:\s+-if-exist\s+){0,1}([^-\n]*)(?:-Z\s+){0,1}(.*)'; #block to indicate dimension type and other characteristics of the parameter
my $XPR_END_BLOCK   = '^end';

#my $XPR_START_GROUP = '^group\s+([^:]*):(.*)'; # group name : units
my $XPR_START_GROUP = '^group\s+([^\(\:\s]+)\s*\(?([^\)]*)\)?\s*:\s*(.*)'; # group name (acronym) : units
my $XPR_END_GROUP   = '^end';

my $XPR_ACRO        = '^\s*[a-zA-Z_\-0-9]+\s*$';

#my $XPR_PARAM      = '^([^:]+):{0,1}([^:]*):{0,1}(.*)'; #name : comment : units
my $XPR_PARAM       = '^(?!group\s+)([^:]+):{0,1}([^:]*):{0,1}(.*)'; #name : comment : units

my $XPR_NAME        = '^([^\[\(\=\s]+)(.*)';
my $XPR_NAME_CONS   = '^\[(.+)\](.*)';
my $XPR_NAME_QUOTA  = '^\((.+)\)(.*)';
my $XPR_NAME_FUNC   = '^\=(.+)';

my $XPR_FUNC_SYM    = '(\+|-){0,1}([^\+-]+)';
my $XPR_FUNC_TERM   = '(.*)(->|<-)(.*)'; # ( a ->/<- b )
my $XPR_FUNC_STAT   = '(\+|-){0,1}\({0,1}([^\)]+)\){0,1}';

my $XPR_UNITS_NAME  = '(\D)';
my $XPR_UNITS_VALUE = '([^:]+)';
########### REGULAR EXPRESSIONS ##########################


########### VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(process_memLayout);
########### VARIABLES ##########################


########### FUNCTIONS ##########################

#remove starting and trailing whitespaces
sub trim{
    my $temp = $_[0];
    $temp =~ s/^\s+//;
    $temp =~ s/\s+$//;
    chomp($temp);
    return $temp;
}
sub trimPar{
    my $temp = $_[0];
    $temp = trim($temp);
    $temp =~ s/^\s*\(\s*//;
    $temp =~ s/s*\)s*$//;
    chomp($temp);
    return $temp;
}

sub check_directive{
    my ( $include, $directives, $blk_include ) = @_;

    trim($include);
    my @includes = split(/\s+/, $include);
    #print join("<>",@includes) . "\n";
    foreach my $inc (@includes){
        #print "\t$directives ----> $inc ----> return ";
        if ( !($directives =~ /${inc}(\s|$)/) ){
            #print "0\n"; 
            return 0; 
        }
    }
    #print "1\n"; 
    $$blk_include = $include;
    return 1;
}

sub insert_term{
    my ($func_matrix, $sign, $term, $lst_param) = @_;

    my($compo1, $direction, $compo2, $params);
    #print "TERM: $term\n";
    if( $term =~ /$XPR_FUNC_TERM/ ){
        $compo1 = trim($1);
        $direction = $2;
        $compo2 = trim($3);
    }else{ print "ERROR malformed func:\"$term\"\n"; exit 1; }

    #get direction
    my $ss=0;
    if(    $direction eq "->" ){ $ss = 1;                                                            }
    elsif( $direction eq "<-" ){ $ss = 0;                                                            }
    else{                        print "ERROR in symbol:\"$direction\" in func:\"$term\"\n"; exit 1; }

    #get component 1 and sign1
    my @compo1_array = ();
    my @sign1_array  = ();
    #print "COMPO1: $compo1\n";
    while ( $compo1 =~ /$XPR_FUNC_SYM/g ){
        my $sign1_ele = "+";
        if( $1 ){ $sign1_ele = $1 };
        my $element  = $2;
        #print "Element: $element Sign1: $sign1_ele\n";
        
        if( $sign eq '-' ){
            if( $sign1_ele eq '+'){ $sign1_ele='-'; }
            else{                   $sign1_ele='+'; }
        }

        if( $element =~ /(.)\.(.)/ ){ 
            my $first = $1;
            my $compo = $2;
            foreach my $param ( values(%$$lst_param) ) {
                if( $param->getSigla() =~ /^${first}.$/ ){
                    foreach my $compoPar ( keys(%{$param->getComponents()}) ){
                        if( $compoPar eq $compo ){ 
                            push( @compo1_array, $param->getSigla() . $compo ); 
                            push(@sign1_array, $sign1_ele);
                        }
                    }
                }
            }
        }else{ 
            push( @compo1_array, $element  ); 
            push( @sign1_array,   $sign1_ele );
        }
        #print "ARRAY1: @compo1_array\n";   
    }

    #get component 2 and sign2
    my @compo2_array = ();
    my @sign2_array   = ();
    #print "COMPO2: $compo2\n";
    while ( $compo2 =~ /$XPR_FUNC_SYM/g ){
        my $sign2_ele = "+";
        if( $1 ){ $sign2_ele = $1 };
        my $element  = $2;
        #print "Element: $element Sign2: $sign2_ele\n";
        
        if( $sign eq '-' ){
            if( $sign2_ele eq '+'){ $sign2_ele='-'; }
            else{                   $sign2_ele='+'; }
        }

        if( $element =~ /(.)\.(.)/ ){ 
            my $first = $1;
            my $compo = $2;
            foreach my $param ( values(%$$lst_param) ) {
                if( $param->getSigla() =~ /^${first}.$/ ){
                    foreach my $compoPar ( keys(%{$param->getComponents()} ) ){
                        if( $compoPar eq $compo ){ 
                            push( @compo2_array, $param->getSigla() . $compo ); 
                            push(@sign2_array, $sign2_ele);
                        }
                    }
                }
            }
        }else{ 
            push( @compo2_array, $element  ); 
            push( @sign2_array,   $sign2_ele );
        }
        #print "ARRAY2: @compo2_array\n";   
    }

    #combine all in one array to avoid doing vector iterations and simplify maths
    my (@compo1, @compo2, @sign1, @sign2, @dir); 

    for my $index1 ( 0 .. $#compo1_array ) {
        my $elem1 = $compo1_array[$index1]; 
        my $sign1 = $sign1_array[$index1];
        for my $index2 ( 0 .. $#compo2_array ) {
            my $elem2 = $compo2_array[$index2]; 
            my $sign2 = $sign2_array[$index2];

            push( @compo1, $elem1 );
            push( @compo2, $elem2 );
            push( @sign1 , $sign1 );
            push( @sign2 , $sign2 );
            push( @dir   , $ss    );
        }
    }

    $$func_matrix{dir}    = \@dir;
    $$func_matrix{compo1} = \@compo1;
    $$func_matrix{sign1}  = \@sign1;
    $$func_matrix{compo2} = \@compo2;
    $$func_matrix{sign2}  = \@sign2;
}


sub process_memLayout{
    my $nml_def   = shift; #input file
    my $lst_group = shift; #output for groups
    my $lst_param = shift; #output for parameters
    my $lst_sta   = shift; #output for stadistics about parameters
    my $lst_const = shift; #output for constituents
    my $cpp_defs  = shift; #include directives
    my $DEBUG   = shift; #debug mode

    my $blk_dim;
    my $blk_type;
    my $blk_subtype;
    my $blk_group;
    my $blk_include;
    my $blk_is_include;
    my $blk_z;


    #open the default file
    open NML_DEF, "<", "$nml_def" or die "$nml_def cannot be opened: $!";
    my @lines_def = <NML_DEF>;
    close(NML_DEF);

    #iterate through default namelist substituting the parameters for new ones
    my $blk_group_idx = 0;
    my $param_idx     = 0;
    my $const_idx     = 1;
    undef $blk_dim;
    undef $blk_type;
    $blk_subtype = 'pel';
    undef $blk_group;
    undef $blk_include;
    $blk_is_include = 1;;
    undef $blk_z;
    foreach my $line_raw (@lines_def){
        #first remove any commentary from line
        my @line_decomposed = ( $line_raw =~ /$XPR_GLOBAL_COMMENT/ );

        my $line = trim($line_decomposed[0]);

        if( $line =~ /$XPR_START_BLOCK/ ){ 
            if( $1 ) { $blk_dim        = $1;                                            }
            if( $2 ) { $blk_type       = $2;                                            }
            if( $3 ) { $blk_subtype    = $3;                                            }
            if( $4 ) { $blk_is_include = check_directive($4, $cpp_defs, \$blk_include); }
            if( $5 ) { $blk_z          = $5;                                            } 
        }
        elsif( $line =~ /$XPR_END_BLOCK/ && !$blk_group ){
            undef $blk_dim; undef $blk_type; $blk_subtype='pel'; undef $blk_include; $blk_is_include=1; undef $blk_z; 
        }
        elsif( $line =~ /$XPR_START_GROUP/ ){
            #extract name and units
            if( $DEBUG ){ print "Processing GROUP: $line_raw"; }
            my ( $obj, $name, $units, $par_compo, @par_const );
            my ( $in_name, $in_acro, $in_units ) = ($1, $2, $3);
            if ( ! process_name($in_name, $blk_dim, $blk_subtype, \$name, \$units, undef, undef, undef, undef) ){ print "WARNING: Group not found: $in_name\n"; next; }
            if( ! ($in_acro =~ /$XPR_ACRO/) ){ print "ERROR: Group \"$name\" has not a valid Acronym: $in_acro \n"; exit 1; }
            process_units($in_units, $units, undef, \$par_compo, undef, undef, \@par_const);
            
            $blk_group=$name;
            #add to groups list
            if( $blk_is_include ){
                $obj = new Group( $blk_group, $in_acro, $blk_dim, $blk_type, $blk_subtype, $par_compo, $blk_include, $blk_z, $blk_group_idx);
                $blk_group_idx++;
                $$lst_group{$blk_group} = $obj;
                #add constituents to list
                foreach my $const (@par_const){
                    if ( ! exists $$lst_const{$const} ){ 
                        $$lst_const{$const} = $const_idx;
                        #print "\$\$lst_const{$const} = $const_idx\n";
                        $const_idx++;
                    }
                }
                if ( $DEBUG ){ $obj->print(); }
            }#else{ print "INFO: not included GROUP $blk_group\n" }
        }
        elsif( $line =~ /$XPR_END_GROUP/ )  { 
            undef $blk_group; 
        }
        elsif( $line =~ /$XPR_PARAM/ ){
            #extract name, comment and units
            if( $DEBUG ){ print "Processing PARAM: $line_raw"; }
            if( $blk_is_include ){
                my ( $units, $size );
                my ( $obj, $par_name, $par_unit, $par_type, $par_compo, $par_compoEx, $par_comm, $par_func, $par_group, $par_quota, @par_const );
                
                if ( ! process_name($1, $blk_dim, $blk_subtype, \$par_name, \$units, \$par_func , \$par_quota, \$lst_param, \$lst_sta) ){ print "WARNING: Parameter not found: $1\n"; next; }

                if( $par_quota && $units ){
                    my $par_name_src = $par_name;
                    $par_comm  = trim($2);
                    my $temp_unit =$3;
                    my @compo_array = ( $units =~ /$XPR_UNITS_NAME/g ); 
                    $units = undef;
                    $par_group = $blk_group;
                    if ( $par_group ){ $size = process_units($temp_unit, $units, \$par_unit, \$par_compo, \$par_compoEx, $$lst_group{$par_group}, \@par_const); }
                    else{              $size = process_units($temp_unit, $units, \$par_unit, \$par_compo, \$par_compoEx, undef                  , \@par_const); }
                    $par_type    = 'diaggrp';

                    foreach my $compo_tmp ( @compo_array ){
                        $par_name = $par_name_src . $compo_tmp;

                        $obj = new Parameter( $par_name, $blk_dim, $par_type, $blk_subtype, $blk_include, $blk_z, $par_unit, $par_compo, $par_compoEx, $par_comm, $par_func, $par_group, $par_quota, $param_idx );
                        #print " " . $obj->getSigla() . " - " . $obj->getIndex() . "\n";
                        $param_idx++;
                        $$lst_param{$par_name} = $obj;
                        $$lst_sta{"${par_type} ${blk_dim}d ${blk_subtype}"} += $size;
                        #@{$lst_const}{@par_const} = 0; # create the list of constituents
                        if ( $DEBUG ){ $obj->print(); }
                    }
                }else{
                    $par_comm  = trim($2);
                    $par_group = $blk_group;
                    if ( $par_group ){ $size = process_units($3, $units, \$par_unit, \$par_compo, \$par_compoEx, $$lst_group{$par_group}, \@par_const); }
                    else{              $size = process_units($3, $units, \$par_unit, \$par_compo, \$par_compoEx, undef                  , \@par_const); }
                    
                    if( $par_quota )                                                       { $par_type = 'diaggrp'                      } #name(quota)
                    elsif( $3 && ($blk_type eq 'variable') && ( $blk_dim =~ /(2|3)/ ) )    { $par_type = 'diagnos'; } #if 3d/2d has units and is a variable => insert in diagnos group
                    else                                                                   { $par_type = $blk_type;                     };#normal parameter
                    
                    $obj = new Parameter( $par_name, $blk_dim, $par_type, $blk_subtype, $blk_include, $blk_z, $par_unit, $par_compo, $par_compoEx, $par_comm, $par_func, $par_group, $par_quota, $param_idx );
                    $param_idx++;
                    #print " " . $obj->getSigla() . " - " . $obj->getIndex() . "\n";
                    $$lst_param{$par_name} = $obj;
                    $$lst_sta{"${par_type} ${blk_dim}d ${blk_subtype}"} += $size;
                    #@{$lst_const}{@par_const} = 0; # create the list of constituents
                    if ( $DEBUG ){ $obj->print(); }
                }
            }#else{ print "INFO: not included PARAM $par_name\n" }
        }#else { print "INFO: line without processing: $line_raw"; }
    }


    #assing values to constituents
    #$$lst_const{'c'} = 1; $$lst_const{'n'} = 2; $$lst_const{'p'} = 3; $$lst_const{'f'} = 5; $$lst_const{'l'} = 4; $$lst_const{'s'} = 6; $$lst_const{'h'} = 7;
    #$$lst_const{'c'} = 1; $$lst_const{'n'} = 2; $$lst_const{'p'} = 3; $$lst_const{'l'} = 4; $$lst_const{'s'} = 5; $$lst_const{'h'} = 6;
    #$$lst_const{'c'} = 1; $$lst_const{'n'} = 2; $$lst_const{'p'} = 3; $$lst_const{'l'} = 4; $$lst_const{'s'} = 5;

}


sub process_name{
    my ( $line, $blk_dim, $blk_subtype, $name_ref, $units_ref, $func_ref, $quota_ref, $lst_param, $lst_sta ) = @_;

    #initialize
    $$name_ref=undef; $$units_ref=undef; $$func_ref=undef; $$quota_ref=undef;

    #get only the name part
    if( $line && $line =~ /$XPR_NAME/ ){ $$name_ref = $1; $line = $2; }else{ return 0; }
    #get constituents part
    if( $line =~ /$XPR_NAME_CONS/ ){ $$units_ref = $1; $line = $2; }
    #get quota part
    if( $line =~ /$XPR_NAME_QUOTA/ ){ $$quota_ref = $1; $line = $2; }
    #get function part
    if( $line && $line =~ /$XPR_NAME_FUNC/ ){ 
        my $xpr = trim($1);
        my ( @function_dir, @function_compo1, @function_sign1, @function_compo2, @function_sign2 );
        while( $xpr =~ /$XPR_FUNC_STAT/g ){
            my ($sign, $term , %func_matrix);
            if( $1 ){ $sign=trimPar($1); }else{ $sign = '+'; }
            if( $2 ){ $term=trimPar($2); }else{ $term = ""; }
            
            #print "SIGN: $sign TERM: $term\n";
            #insert the term
            insert_term(\%func_matrix, $sign, $term, $lst_param);
            push( @function_dir   , @{$func_matrix{dir}} ); 
            push( @function_compo1, @{$func_matrix{compo1}} ); 
            push( @function_sign1  , @{$func_matrix{sign1}} ); 
            push( @function_compo2, @{$func_matrix{compo2}} ); 
            push( @function_sign2  , @{$func_matrix{sign2}} ); 
        }
        ${$$func_ref}{xpr}    = $xpr;
        ${$$func_ref}{dir}    = \@function_dir;
        ${$$func_ref}{compo1} = \@function_compo1;
        ${$$func_ref}{sign1}  = \@function_sign1;
        ${$$func_ref}{compo2} = \@function_compo2;
        ${$$func_ref}{sign2}  = \@function_sign2;
        ${$$lst_sta}{"select ${blk_dim}d ${blk_subtype}"} += scalar(@function_compo1); 
        #print Dumper(\%{$$func_ref}) , "\n";
    }

    #if($$name_ref){  print " NAME: $$name_ref"; }; if($$units_ref){ print " UNITS: $$units_ref"; }; if($$quota_ref){ print " QUOTA: $$quota_ref"; }; if($$func_ref){  print " FUNC: $$func_ref"; }; print "\n";
    return 1;
}


sub process_units{
    my ( $values, $names, $unit_ref, $compo_ref, $compoEx_ref, $group_obj, $par_const ) = @_;
    my ( @name_list, @value_list );

    if ( $names ){
        #print "UNIT NAMES $names , UNIT VALUES $values\n";
        #get the name and values
        @name_list  = ( $names =~ /$XPR_UNITS_NAME/g );
        #has several components (ex: carbon, nitrogen,...) each one with their unit
        my %child_compo = ();
        my %child_compo_ex = ();
        
        if( $name_list[0] eq '-' ){
            #the units have to be extracted from parent
            if( ! $group_obj ){ print "ERROR: group parent does not exists\n"; exit 1; }
            %child_compo = %{$group_obj->getComponents()}; #copy all elements from parent
            splice @name_list, 0, 1; #remove first element (-)
            foreach my $name (@name_list){ 
                if( exists $child_compo{$name} ){
                    ${child_compo_ex}{$name} = $child_compo{$name}; #add to excluded componentes
                    delete($child_compo{$name}); #remove from included components
                }
            }
        }else{
            #units and values are listed
            @value_list = ( $values =~ /$XPR_UNITS_VALUE/g );
            if( $#value_list != $#name_list ){ 
                print "ERROR: Components and units have not same size ($#name_list != $#value_list) !!\n"; 
                print "\tlists: (@value_list != @name_list)\n";
                exit 1;
            }
            my $index = 0;
            foreach my $name (@name_list){ 
                $child_compo{$name} = trim($value_list[$index++]);
                push( @$par_const, $name );
            }
        }
        $$unit_ref    = undef;
        $$compo_ref   = \%child_compo;
        if( keys(%child_compo_ex) != 0 ){ $$compoEx_ref=\%child_compo_ex; }else{ $$compoEx_ref=undef; }
        return keys(%child_compo);
    }else{
        if ( $group_obj ){
            #copy parent units
            $$unit_ref    = undef;
            $$compo_ref   = \%{$group_obj->getComponents()};
            $$compoEx_ref = undef;
            return keys(%{$group_obj->getComponents()});
        }else{
            #only one unit which is the parameter unit
            $$unit_ref    = trim($values);
            $$compo_ref   = undef;
            $$compoEx_ref = undef;
            return 1;
        }
    }
}


1;
