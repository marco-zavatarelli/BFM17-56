#!/usr/bin/perl -w

# DESCRIPTION
#   Generate f.90 and .h files
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

package print_f90;

use strict;
use warnings;
use Exporter;
use F90Namelist;
use List::Util qw[min max];
use Data::Dumper;

use classes;



########### VARIABLES TO EXPORT ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(print_f90);
########### VARIABLES TO EXPORT ##########################


########### VARIABLES GLOBAL ##########################
my $VERBOSE = 0;
my $SPACE = "        ";
my $MAX_CHAR_PER_LINE = 76;
my $STRING_INDEX=1;
my %STRING_INDEX_ARRAY = ();
my ( $LST_PARAM, $LST_GROUP, $LST_STA, $LST_CONST, $LST_INDEX, $LST_NML );
########### VARIABLES GLOBAL ##########################





########### REGULAR EXPRESSIONS ##########################
my $dispatch = {
    #COMMENT LINE FOR PROGRAMMERS
    '\%\!.*' => \&func_COMMENT,
    #DESCRIPTION
    '\%(3|2)d-(diagnos|state)-(pel|ben)-desc(?:\s|\n)'=> \&func_DESC ,
    '\%(3|2)d-(diaggrp)-(pel|ben)-desc(?:\s|\n)' => \&func_DESC_DIAGG ,
    #NR
    '([^\%]*)\%(3|2)d-(flux|diagnos|state)-(pel|ben)-nr(?:\s|\n)'=> \&func_NR ,
    #ARRAY
    '([^\%]*)\%(pel|ben)-field-array\s*(\S*)(?:\s|\n)'=> \&func_ARRAY_FIELD ,
    #PP
    '([^\%]*)\%(3|2)d-(diagnos|state)-(pel|ben)-pp(?:\s|\n)'=> \&func_PP ,
    '([^\%]*)\%(3|2)d-(diaggrp)-(pel|ben)-pp(?:\s|\n)'=>       \&func_PP_DIAGG ,
    '\%(3|2)d-(diaggrp)-(pel|ben)-assign-pp(?:\s|\n)'=>        \&func_PP_ASSIGN ,
    #POINTER
    '\%(3|2)d-(diagnos|state)-(pel|ben)-pointer(?:\s|\n)'=>    \&func_POINT ,
    '([^\%]*)\%(3|2)d-(diaggrp)-(pel|ben)-pointer(?:\s|\n)'=>    \&func_POINT_DIAGG ,
    '([^\%]*)\%(pel|ben)-field-pointer\s*(\S*)(?:\s|\n)'=>    \&func_POINT_FIELD ,
    '\%(3|2)d-(diagnos|state)-(pel|ben)-alloc-pointer(?:\s|\n)' =>    \&func_POINT_ALLOC ,
    '\%(3|2)d-(diaggrp)-(pel|ben)-alloc-pointer(?:\s|\n)' =>    \&func_POINT_ALLOC_DIAGG ,
    '\%(pel|ben)-alloc-pointer-field\s*(\S*)(?:\s|\n)' =>    \&func_POINT_ALLOC_FIELD ,
    #HEADER
    '\%(3|2)d-(diagnos|state)-(pel|ben)-header(?:\s|\n)' =>    \&func_HEADER ,
    '\%(3|2)d-(diaggrp)-(pel|ben)-header(?:\s|\n)' =>    \&func_HEADER_DIAGG ,
    '\%(pel|ben)-field-header\s*(\S*)(?:\s|\n)' =>    \&func_HEADER_FIELD ,
    #STRING
    '\%(3|2)d-(diagnos|state|flux)-(pel|ben)-string(?:\s|\n)'=>    \&func_STRING ,
    '\%(3|2)d-(diaggrp)-(pel|ben)-string(?:\s|\n)'=>    \&func_STRING_DIAGG ,
    '\%(pel|ben)-field-string\s*(\S*)(?:\s|\n)'=>    \&func_STRING_FIELD ,
    '\%(3|2)d-(diagnos|state|flux)-(pel|ben)-string-index(?:\s|\n)'=>    \&func_STRING_INDEX ,
    '\%(pel|ben)-string-index-field\s*(\S*)(?:\s|\n)'=>    \&func_STRING_INDEX_FIELD ,
    #ALLOC
    '\%(3|2)d-(diagnos|state)-(pel|ben)-alloc(?:\s|\n)' => \&func_ALLOC ,
    '\%(1|3|2)d-(intvar|variable|variables)-(pel|ben)-alloc(?:\s|\n)'=>    \&func_ALLOC_INTVAR ,
    #FLUX
    '\%(3|2)d-(flux)-(pel|ben)-alloc(?:\s|\n)' => \&func_FLUX_ALLOC ,
    '\%(3|2)d-(flux)-(pel|ben)-fill(?:\s|\n)' => \&func_FLUX_FILL ,
    #CONSTITUENT
    '([^\%]*)\%constituent(?:\s|\n)' => \&func_CONSTITUENT ,
    #GROUP    
    '\%(3|2)d-group-(pel|ben)-header(?:\s|\n)'                => \&func_GROUP_HEADER ,
    '([^\%]*)\%(3|2)d-group-(pel|ben)-parameter(?:\s|\n)'     => \&func_GROUP_PARAMETER ,
    '([^\%]*)\%group-(pel|ben)-calc(?:\s|\n)'                 => \&func_GROUP_CALC,
    '([^\%]*)\%(3|2)d-group-(pel|ben)-function-name(?:\s|\n)' => \&func_GROUP_FUNCTION_NAME ,
    '\%(3|2)d-groupfunctions-(pel|ben)(?:\s|\n)'              => \&func_GROUP_FUNCTIONS ,
    #INTVAR/VARIABLE
    '([^\%]*)\%(3|2|1)d-(intvar|variable)-(pel|ben)(?:\s|\n)' => \&func_VARIABLE ,
    #INIT
    '\%value-init-calc-(pel|ben)\s*(\S*)(?:\s|\n)'       => \&func_INIT_VALUE_CALC ,
    '([^\%]*)\%(3|2)d-(state)-(pel|ben)-Initpp(?:\s|\n)' => \&func_INIT_PP ,
    '\%init-func-constituents(?:\s|\n)'                  => \&func_INIT_FUNC_CONST ,
    '\%(3|2)d-(state)-(pel|ben)-InitDefault(?:\s|\n)'      => \&func_INIT_DEFAULT ,
    '\%(3|2)d-(state)-(pel|ben)-InitSets(?:\s|\n)'         => \&func_INIT_SETS ,
    '\%(3|2)d-(state)-(pel|ben)-InitInternal(?:\s|\n)'     => \&func_INIT_INTERNAL ,
    '\%(3|2)d-state-(pel|ben)-func-zeroing(?:\s|\n)'       => \&func_INIT_FUNC_ZERO ,
};
########### REGULAR EXPRESSIONS ##########################


########### FUNCTIONS ##########################

sub print_f90 {

    my ($templ_mem, $output, $lst_group, $lst_param, $lst_sta, $lst_const, $lst_index, $lst_nml, $verbose) = @_;

    $LST_PARAM = $lst_param;
    $LST_GROUP = $lst_group;
    $LST_STA   = $lst_sta;
    $LST_CONST = $lst_const;
    $LST_INDEX = $lst_index;
    $LST_NML = $lst_nml;
    $VERBOSE   = $verbose;

    #open the template file
    open TEMPL_MEM, "<", "$templ_mem" or die "$templ_mem cannot be opened: $!";
    my @lines_template = <TEMPL_MEM>;
    close(TEMPL_MEM);

    #open the output file
    open my $OUT, ">", "$output" or die "$output cannot be opened: $!";

    foreach my $line_raw (@lines_template){
        my $replaced = 0;
        foreach my $key (sort keys %$dispatch){
            if( $line_raw =~ /$key/ ){
                my ($arg1, $arg2, $arg3, $arg4);
                $arg1 = $1; $arg2=$2; $arg3=$3; $arg4=$4;
                if   ($arg4) { &{$$dispatch{$key}}($OUT, $arg1, $arg2, $arg3, $arg4); }
                elsif($arg3) { &{$$dispatch{$key}}($OUT, $arg1, $arg2, $arg3);        }
                elsif($arg2) { &{$$dispatch{$key}}($OUT, $arg1, $arg2);               }
                elsif($arg1) { &{$$dispatch{$key}}($OUT, $arg1);                      }
                else         { &{$$dispatch{$key}}($OUT);                             }
                print $OUT "\n";
                $replaced=1;
                if( $VERBOSE ){ print " $line_raw"; }
            }
            if ($replaced){ last; }
        }
        if(!$replaced){ print $OUT $line_raw; }
    }

    #close the output file
    close($OUT);

}


sub printList{
    my ($file, $line, $initialLine, $separator) = @_;
    #split the line if necessary and print it
    
    if( !$separator ){ $separator = ", "; }
    my $len_sep = length($separator);

    
    #print first line and first element
    print $file "$initialLine";
    print $file $${line}[0];
    my $count = length($initialLine) + length($${line}[0]);
    
    #print the rest of elements with separator
    foreach my $index ( 1 .. $#{$line} ){
        my $ele = $${line}[$index];
        $count += length($ele) + $len_sep;
        if( $count > $MAX_CHAR_PER_LINE ){
            #inset a new line
            print $file "${separator}&\n${SPACE} ${ele}"; 
            $count = length($ele) + $len_sep;
        }else{ print $file "${separator}${ele}"; }
    }
}



sub decreaseUnit{
    my ( $unit ) = @_;

    my @params = ( $unit =~ m/(.*)([0-9]+$)/);
    my $len = $params[1] - 1;
    if( $len == 0 ){
        return $params[0] . "/day";
    }else{
        return $params[0] . $len . "/day";
    }
}


sub sizeGroup{
    my ( $dim ) = @_;

    my $size = 0;
    foreach my $root (keys %$LST_PARAM){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() && $param->getQuota() ){
            #print $root . " -> " . $param->getQuota() . "\n";
            foreach my $member (keys %$LST_PARAM){
                my $param2 = $$LST_PARAM{$member};
                if( defined($param2->getGroup()) && $param2->getGroup() eq $param->getQuota() ){ $size++; }
            }
        }
    }

    return $size;
}


sub sizeQuota{
    my ( $param ) = @_;

    my $size = 0;
    if( $param->getQuota() ){
        my $group_name = $param->getQuota();
        foreach my $name2 (keys %$LST_PARAM){
            my $param2 = $$LST_PARAM{$name2};
            if ( $param2->getGroup() && $group_name eq $param2->getGroup() ){ $size++; }
        }
    }else{ $size++; }

    return $size;
}

sub sizeDimType {
    my ($dim, $type, $subt) = @_;

    my $size = 0;

    foreach my $name (keys %$LST_PARAM){
        my $param       = $$LST_PARAM{$name};
        my $dim_var     = $param->getDim();
        my $type_var    = $param->getType();
        my $subtype_var = $param->getSubtype();
        if( $type_var eq 'diaggrp' ){ $type_var = 'diagnos'; }

        if( $dim_var && $dim_var == $dim 
            && $type_var eq $type 
            && $subt eq $subtype_var ){ 
            my $sizeQuota = sizeQuota($param);
            if( $param->getComponents() ){
                $size += (keys(%{$param->getComponents()}) * ($sizeQuota));
                #print " ${name}x" . (keys(%{$param->getComponents()}) * ($sizeQuota));
            }else{
                $size += $sizeQuota;
                #print " $name";
            }
        }
    }

    #print " -- $size $dim $type -- ";
    return $size;
}

# sub sizeWithCompo {
#     my ($param) = @_;

#     if( $param->getComponents() ){ return keys( %{$param->getComponents()} );  }
#     else{ return 1; }
# }

sub checkDimType {
    my ($dim, $type, $subt) = @_;

    foreach my $key (keys %$LST_PARAM){
        my $param       = $$LST_PARAM{$key};
        my $dim_var     = $param->getDim();
        my $type_var    = $param->getType();
        my $subtype_var = $param->getSubtype();

        if( $dim_var == $dim 
            && $type_var eq $type 
            && $subtype_var eq $subt ){ 
            #print " $dim_var == $dim && $type_var eq $type\n";
            return 1; 
        }   
    }

    #print " != $dim && != $type\n";
    return 0;
}


sub func_COMMENT {}


###########################################  ALLOCATE_MEM FUNCTIONS ######################


sub func_ALLOC {
    my ( $file, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_ALLOC: "; }

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type, $subt)  ){ return; }
    
    my $j = "";
    my $TYPE = uc($type);
    if ( $dim == 2 ){ $j = "_XY"; }

    print $file "${SPACE} \n";
    print $file "${SPACE}allocate(D${dim}${TYPE}(1:NO_D${dim}_BOX_${TYPE}S,1:NO_BOXES$j),stat=status)\n";
    print $file "${SPACE}if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}$TYPE\")\n";
    print $file "${SPACE}D${dim}${TYPE} = ZERO\n";
    if ( $type eq "state" && $subt eq "pel" ) {
        print $file "#ifndef EXPLICIT_SINK\n";
        print $file "${SPACE}  allocate(D${dim}SOURCE(1:NO_D${dim}_BOX_STATES,1:NO_BOXES$j),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SOURCE\")\n";
        print $file "${SPACE}  D${dim}SOURCE = ZERO\n";
        print $file "#else\n";
        print $file "${SPACE}  allocate(D${dim}SOURCE(1:NO_D${dim}_BOX_STATES,1:NO_D${dim}_BOX_STATES,1:NO_BOXES$j),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SOURCE\")\n";
        print $file "${SPACE}  D${dim}SOURCE = ZERO\n";
        print $file "${SPACE}  allocate(D${dim}SINK(1:NO_D${dim}_BOX_STATES,1:NO_D${dim}_BOX_STATES,1:NO_BOXES$j) ,stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SINK\")\n";
        print $file "${SPACE}  D${dim}SINK = ZERO\n";
        print $file "#endif\n";
        print $file "${SPACE}allocate(D${dim}STATETYPE(1:NO_D${dim}_BOX_${TYPE}S ),stat=status)\n";
        print $file "${SPACE}if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\",\"D${dim}STATETYPE\")\n";
        print $file "${SPACE}D${dim}STATETYPE = ZERO\n";
        print $file "#ifdef BFM_NEMO\n";
        print $file "${SPACE}  allocate(D${dim}STATEOBC(1:NO_D${dim}_BOX_STATES),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\",\"D${dim}STATEOBC\")\n";
        print $file "${SPACE}  D${dim}STATEOBC = ZERO\n";
        print $file "#endif\n";
    }
}

sub func_POINT_ALLOC  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_POINT_ALLOC: "; }

    my $line = "";
    my $TYPE = uc($type);

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $param->getSubtype eq $subt ){

            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        $line .= "$nameC => D${dim}${TYPE}(pp$nameC,:); $nameC=ZERO\n";
                    }
                }
            }else{
                $line .= "$name => D${dim}${TYPE}(pp$name,:); $name=ZERO\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}

sub func_PP_ASSIGN  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_PP_ASSIGN: "; }
    
    my $line = "";
    my $TYPE = uc($type);

    #calculate lenght
    my $len=0;

    #foreach my $root (sort keys %$LST_PARAM){
    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() 
            && $param->getQuota()
            && $param->getSubtype eq $subt ){
            #print $root . " -> " . $param->getQuota() . "\n";
            foreach my $member ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                my $param2 = $$LST_PARAM{$member};
                if( defined($param2->getGroup()) && $param2->getGroup() eq $param->getQuota() ){
                    $line .= "pp${root}(ii${member})=" . ++($$LST_STA{"diagnos ${dim}d $subt"}) . "\n";
                }
            }
        }
    }

    if( $line ){ print $file $line; }
}

sub func_POINT_ALLOC_DIAGG  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_PP_ALLOC_DIAG: "; }
    
    my $line = "";

    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() 
            && $param->getQuota()
            && $param->getSubtype eq $subt ){
            #print $root . " -> " . $param->getQuota() . "\n";
            my @namesMem = ();
            foreach my $member ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                my $param2 = $$LST_PARAM{$member};
                if( defined($param2->getGroup()) && ($param2->getGroup() eq $param->getQuota()) ){
                    push (@namesMem, $param2->getSigla() );
                }
            }
            if( $namesMem[0] && $namesMem[$#namesMem]){
                my $mem1 = $namesMem[0];
                my $mem2 = $namesMem[$#namesMem];
                $line .= "$root => D${dim}DIAGNOS(pp${root}(ii${mem1}): pp${root}(ii${mem2}),:)\n";
                $line .= "$root=ZERO\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_POINT_ALLOC_FIELD  {
    my ( $file, $subt, $spec) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_POINT_ALLOC_FIELD: "; }

    my $line_ini = "";
    my $line_par = "";
    
    my $SPEC       = uc($spec);
    my $spec_short = substr($spec,0,3);

    my $n = 0;
    my $m = sizeDimType(3, 'state', $subt);
    if( exists $$LST_STA{"state 3d ${subt} index"} ){
        $n = $$LST_STA{"state 3d ${subt} index"};
        $$LST_STA{"state 3d ${subt} index"} += $m;
    }else{
        $n = $$LST_STA{"diagnos 2d ${subt}"};
        $$LST_STA{"state 3d ${subt} index"} = $n + $m;
    }
    $$LST_STA{"state 2d ${subt} ${spec} index"} = $n;

    $line_ini .= "${SPACE}PEL${SPEC} => D2DIAGNOS(${n}+1:${n}+${m},:); PEL${SPEC}=ZERO\n";
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( 3 == $param->getDim() 
            && 'state' eq $param->getType()
            && $subt eq $param->getSubtype() ){

            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        $line_par .= "${SPACE}j${spec_short}${nameC} => D2DIAGNOS(${n}+pp${nameC},:); j${spec_short}${nameC}=ZERO\n";
                    }
                }
            }else{
                $line_par .= "${SPACE}j${spec_short}${name} => D2DIAGNOS(${n}+pp${name},:); j${spec_short}${name}=ZERO\n";
            }
        }
    }

    if( $line_par ){ print $file $line_ini . $line_par; }
}


sub func_ALLOC_INTVAR  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_ALLOC_INTVAR: "; }
    
    my $line = "";
    my $j = "";
    my $TYPE = uc($type);
    if ( $dim == 2 ){ $j = "_XY"; }

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){

            my $mode = 0;
            my $l = "";
            my $groupindex = "";
            if( $param->getDim() > 1 ){ $mode += 2; }
            if( $param->getQuota()   ){ 
                $mode += 1; 
                my @group = ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP );
                $groupindex = firstidx { $_ eq $param->getQuota() } @group;
                $l = "ii";
            }

            if    ( $mode == 3 ){ $line .= "${SPACE}allocate(${name}(1:${l}${groupindex},1:NO_BOXES${j}),stat=status); ${name} = ZERO\n"; }
            elsif ( $mode == 2 ){ $line .= "${SPACE}allocate(${name}(1:NO_BOXES${j}),stat=status); ${name} = ZERO\n";                     }
            elsif ( $mode == 1 ){ $line .= "${SPACE}allocate(${name}(1:${l}${groupindex}),stat=status); ${name} = ZERO\n";                }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_FLUX_ALLOC  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_FLUX_ALLOC: "; }

    my $line = "";
   
    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type, $subt)  ){ return; }

    if( $dim == 3 ){
        $line .= "  allocate( D${dim}FLUX_MATRIX(1:NO_D${dim}_BOX_STATES, 1:NO_D${dim}_BOX_STATES),stat=status )\n";
        $line .= "  allocate( D${dim}FLUX_FUNC(1:NO_D${dim}_BOX_FLUX, 1:NO_BOXES),stat=status )\n";
        $line .= "  D${dim}FLUX_FUNC = 0\n";
    }elsif( $dim == 2 ){
        $line .= "  allocate( D${dim}FLUX_MATRIX(1:NO_D${dim}_BOX_STATES, 1:NO_D${dim}_BOX_STATES),stat=status )\n";
        $line .= "  allocate( D${dim}FLUX_FUNC(1:NO_D${dim}_BOX_FLUX),stat=status )\n";
        $line .= "  D${dim}FLUX_FUNC = 0\n";
    }

    if( $line ){ print $file $line; }
}


sub func_FLUX_FILL  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_FLUX_ALLOC: "; }
    my $line = "";


    if( exists $$LST_STA{"flux ${dim}d $subt"} && exists $$LST_STA{"select ${dim}d $subt"} ){

        my %d3flux_func             = ();
        my %d3flux_func_dir         = ();
        my %d3flux_matrix_index     = ();
        my %d3flux_matrix_index_dir = ();
        my $flxindex                = 0;

        #get the fluxes which contain origin-dest
        $flxindex = 1;
        foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
            my $param = $$LST_PARAM{$name};
            if( ${dim} == $param->getDim() 
                && ${type} eq $param->getType()
                && ${subt} eq $param->getSubtype() ){

                my $function = $param->getFunction();
                my $sigla = $param->getSigla();
                for my $indexC ( 0 .. $#{$$function{compo1}} ) {
                    my $sign1  = $$function{sign2}[$indexC];
                    my $dir    = $$function{dir}[$indexC];
                    my $compo1 = $$function{compo1}[$indexC];
                    my $compo2 = $$function{compo2}[$indexC];
                    if( $compo2 eq '*' ){ $compo2 = $compo1; }

                    if( $dir ){ # A-> B
                        if ( exists $d3flux_matrix_index{"pp$compo1, pp$compo2"} ){ 
                            $d3flux_matrix_index{"pp$compo1, pp$compo2"}++ } 
                        else{ 
                            $d3flux_matrix_index{"pp$compo1, pp$compo2"} = 1  
                        }
                        push( @{$d3flux_func{$name}}, "D${dim}FLUX_MATRIX(pp$compo1, pp$compo2)%p(" . $d3flux_matrix_index{"pp$compo1, pp$compo2"} . ")=${sign1}${flxindex}\n " );
                    }else{ # B-> A
                        if ( exists $d3flux_matrix_index{"pp$compo2, pp$compo1"} ){ 
                            $d3flux_matrix_index{"pp$compo2, pp$compo1"}++ } 
                        else{ 
                            $d3flux_matrix_index{"pp$compo2, pp$compo1"} = 1
                        };
                        push( @{$d3flux_func{$name}}, "D${dim}FLUX_MATRIX(pp$compo2, pp$compo1)%p(" . $d3flux_matrix_index{"pp$compo2, pp$compo1"} . ")=${sign1}${flxindex}\n " );
                    }

                    if( $compo1 eq $compo2 ){ # A == B
                        if ( exists $d3flux_matrix_index_dir{"pp$compo1, pp$compo2"} ){ 
                            $d3flux_matrix_index_dir{"pp$compo1, pp$compo2"}++ } 
                        else{ 
                            $d3flux_matrix_index_dir{"pp$compo1, pp$compo2"} = 1  
                        }
                        push( @{$d3flux_func_dir{$name}}, "D${dim}FLUX_MATRIX(pp$compo1, pp$compo2)%dir(" . $d3flux_matrix_index{"pp$compo1, pp$compo2"} . ")=${dir}\n " );
                    }
                }
                $flxindex++;
            }
        }
        #print Dumper(\%d3flux_matrix_index_dir) , "\n"; 
        #print Dumper(\%d3flux_matrix_index) , "\n"; 
        #print Dumper(\%d3flux_func_dir) , "\n"; 
        #print Dumper(\%d3flux_func) , "\n"; 

        # allocate
        foreach my $key ( keys %d3flux_matrix_index ){
            $line .= "  allocate( D${dim}FLUX_MATRIX($key)%p( 1:$d3flux_matrix_index{$key} ) )\n";
            if( exists $d3flux_matrix_index_dir{$key} ){
                $line .= "  allocate( D${dim}FLUX_MATRIX($key)%dir( 1:$d3flux_matrix_index_dir{$key} ) )\n";
            }
        }

        $line .= "\n\n";

        #fill by the order of the functions
        foreach my $key ( keys %d3flux_func ){
            $line .= "  ! $key = " . ${$$LST_PARAM{$key}->getFunction()}{xpr} . "\n";
            $line .= "  @{$d3flux_func{$key}} ";
            if( exists $d3flux_func_dir{$key} ){
                $line .= "@{$d3flux_func_dir{$key}} ";
            }
            $line .= "\n";
        }

    }


    if( $line ){ print $file $line; }
}




###########################################  MODULE_MEM FUNCTIONS ######################


sub func_DESC  {
    my ( $file, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_DESC: "; }

    my $line     = "";
    my $line_par = "";


    $line = sprintf "! %10s %60s %15s", "${dim}d name", "description", "unit\n";
    $line .= "! ";
    foreach (1..10){ $line .=  "-"; } $line .=  " ";
    foreach (1..60){ $line .=  "-"; } $line .=  " ";
    foreach (1..15){ $line .=  "-"; } $line .=  "\n";
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        my $comment = $param->getComment();
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){

            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        my $unit = ${$param->getComponents()}{$const};
                        $line_par .= sprintf "! %10s %60s %15s\n", $nameC, $comment, $unit;
                    }
                }
            }else{
                my $unit = 
                    $line_par .= sprintf  "! %10s %60s %15s\n", $name, $comment, $param->getUnit();
            }
        }
    }

    if( $line_par ){ print $file $line . $line_par; }

}


sub func_ARRAY_FIELD  {
    my ( $file, $before, $subt, $spec) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_ARRAY_FIELD: "; }

    #save the spec in stadistics for further calculations of number of variables
    if( ! exists $$LST_STA{"spec num $subt"} ){ $$LST_STA{"spec num $subt"} = 1; }
    else{ $$LST_STA{"spec num $subt"} += 1; }

    my $spec_upper = uc($spec);
    my $subt_short = uc(substr($subt,0,3));
    
    print $file "$before ${subt_short}${spec_upper}\n";
}

sub func_NR  {
    my ( $file, $before, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_NR: "; }

    my $title  = "${type} ${dim}d ${subt}";
    my $number = 0;
    if( exists $$LST_STA{$title} ){
        if( $title eq "diagnos 2d ${subt}" ){
            $number = $$LST_STA{"${title}"} + sizeGroup($dim) + ( $$LST_STA{"state 3d ${subt}"} * $$LST_STA{"spec num ${subt}"} );
        }
        elsif( $title eq "diagnos 3d ${subt}" ){ 
            $number=sizeDimType($dim, $type, $subt); 
        }
        else{ 
            $number=$$LST_STA{$title}; 
        }
    }else{ $number = 0; }

    print $file "$before${number}\n";
}


sub func_PP  {
    my ( $file, $before, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_PP: "; }

    my @line_par = ();

    if( ! checkDimType($dim, $type, $subt)  ){ return; }

    my $index = 1;
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType() 
            && $subt eq $param->getSubtype() ){

            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        push( @line_par, "pp${nameC}=".$index++ );
                    }
                }
                if( $param->getComponentsEx() && keys(%{$param->getComponentsEx()}) != 0 ){
                    foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                        if( exists ${$param->getComponentsEx()}{$const} ){
                            my $nameC = $name . $const;
                            push( @line_par, "pp${nameC}=0" );
                        }
                    }
                }
            }else{
                push( @line_par, "pp${name}=".$index++);
            }
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $before); print $file "\n"; }

}


sub func_POINT  {
    my ( $file, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_POINT: "; }
    
    my $line     = "";
    my @line_par = ();

    if( ! checkDimType($dim, $type, $subt)  ){ return; }

    $line .= "${SPACE}real(RLEN),public,dimension(:),pointer  :: ";
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){

            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        push( @line_par, ${nameC} );
                    }
                }
            }else{
                push( @line_par, ${name} );
            }
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $line); print $file "\n"; }
}


sub func_CONSTITUENT  {
    my ( $file, $before, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_CONSTITUENT: "; }
    
    my @line_par = ();

    foreach my $key ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
        push( @line_par, "ii" . uc($key) . "=" . $$LST_CONST{$key} );
    }

    if( $#line_par > 0 ){ printList($file, \@line_par, $before); print $file "\n"; }
}


sub func_GROUP_PARAMETER  {
    my ( $file, $before, $dim, $subt ) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_GROUP_PARAMETER: "; }

    foreach my $group_name ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ) {
        my $group = $$LST_GROUP{$group_name};
        if( $dim == $group->getDim() 
            && $subt eq $group->getSubtype() ){
            my @elements = ();
            my $index = 1;
            foreach my $param_name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                my $param = $$LST_PARAM{$param_name};
                if( $param->getGroup() && $group_name eq $param->getGroup() ){ push(@elements, ("ii$param_name=". $index++ ) ); }
            }
            if( $#elements == -1 ){
                print $file "${before}ii$group_name=". scalar(@elements) . "\n" ;
            }else{
                print $file "${before}ii$group_name=". scalar(@elements) . ", " . join(", ",@elements) . "\n" ;
            }
        }
    }
}


sub func_GROUP_CALC {
    my ( $file, $before, $subt ) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_GROUP_CALC: "; }

    my $line = '';
    my $name = '';

    foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$groupname};
        if( $subt eq $group->getSubtype()){

            $line .= "${before}Calc" . $group->getSigla() . "(ii" . $group->getSigla() . ") = .TRUE.\n";
        }
    }

    if( $line ){ print $file $line; }
}




sub func_VARIABLE  {
    my ( $file, $before, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_VARIABLE: "; }

    my $line_par = "";

    foreach my $param_name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$param_name};
        if( $param->getDim() == $dim 
            && ($param->getType() eq $type)
            && ($param->getSubtype() eq $subt) ){

            my $mode  = 0;
            if( $param->getQuota() ){    $mode++; }
            if( $param->getDim () > 1 ){ $mode++; }
            
            if    ($mode == 0) { $line_par .= "${before}                             :: $param_name ! "; }
            elsif ($mode == 1) { $line_par .= "${before},dimension(:),allocatable    :: $param_name ! "; }
            elsif ($mode == 2) { $line_par .= "${before},dimension(:,:),allocatable  :: $param_name ! "; }

            if( $param->getComment() ){ $line_par .= $param->getComment();           }
            if( $param->getUnit() ){    $line_par .= " (" . $param->getUnit()  . ")";}
            $line_par .= "\n";
        }
    }

    if( $line_par ){
        print $file $line_par;
    }
    
}

sub func_DESC_DIAGG  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_DESC_DIAGG: "; }

    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() 
            && $param->getQuota()
            && $subt eq $param->getSubtype() ){
            #print $root . " -> " . $param->getQuota() . "\n";
            foreach my $member ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                my $param2 = $$LST_PARAM{$member};
                if( defined($param2->getGroup()) && $param2->getGroup() eq $param->getQuota() ){
                    printf $file "! %10s %60s %15s\n", "${root}(ii${member})", $param2->getComment(), $param->getUnit();
                }
            }
        }
    }
}


sub func_PP_DIAGG  {
    my ( $file, $before, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_PP_DIAGG: "; }

    my @line_par = ();

    if( ! checkDimType($dim, $type, $subt)  ){ return; }

    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() 
            && $subt eq $param->getSubtype() 
            && $param->getQuota() ){
            push(@line_par,"pp${root}(ii". $param->getQuota() .")" );
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $before); print $file "\n"; }
}


sub func_POINT_DIAGG  {
    my ( $file, $before, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_POINT_DIAGG: "; }
    
    my $line     = "";
    my @line_par = ();

    if( ! checkDimType($dim, $type, $subt)  ){ return; }

    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() 
            && $subt eq $param->getSubtype() 
            && $param->getQuota() ){

            push( @line_par, "${root}" );
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $before); print $file "\n"; }

}

sub func_POINT_FIELD  {
    my ( $file, $before, $subt, $spec) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_POINT_FIELD: "; }

    my @line_par = ();

    my $SPEC = substr($spec,0,3);
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $param->getDim() == 3
            && ($param->getType() eq 'state') 
            && ($param->getSubtype() eq $subt) ){

            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        push( @line_par, "j${SPEC}${nameC}");
                    }
                }
            }else{
                push( @line_par, "j${SPEC}${name}" );
            }
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $before); print $file "\n"; }
}


sub func_GROUP_FUNCTION_NAME  {
    my ( $file, $before, $dim, $subt ) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_GROUP_FUNCTION_NAME: "; }

    my @line_par = ();

    foreach my $group_name ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group = $$LST_GROUP{$group_name};
        if( $group->getDim() == $dim
            && $group->getSubtype() eq $subt ){
            push( @line_par, "pp${group_name}" );
            push( @line_par, "${group_name}" );
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $before); }
}


sub func_GROUP_FUNCTIONS  {
    my ( $file, $dim, $subt ) = @_;
    if ( $VERBOSE ){ print "ModuleMem -> FUNCTION CALLED func_GROUPFUNCTIONS: "; }

    foreach my $pre ("pp", ""){
        foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
            my $group   = $$LST_GROUP{$groupname};

            if( $group->getDim() == $dim
                && $subt eq $group->getSubtype()){
                my @members = ();
                foreach my $paramname ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                    my $param = $$LST_PARAM{$paramname};
                    if ( $param->getGroup() && ($param->getGroup() eq $groupname) ){ push(@members, $param); }
                }
                #order to show in alpahabetic order
                #@members = sort { $a cmp $b } @members;

                my $line = "\n";
                $line .= "${SPACE}" ."function $pre${groupname}(n,constituent,cmax)\n";
                $line .= "${SPACE}" ."\n";
                $line .= "${SPACE}" ."  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
                $line .= "${SPACE}" ."  ! Implicit typing is never allowed\n";
                $line .= "${SPACE}" ."  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
                $line .= "${SPACE}" ."  IMPLICIT NONE\n";
                $line .= "${SPACE}" ."\n";
                if ( $pre eq "pp" ) {
                    $line .= "${SPACE}" ."  integer ::$pre$groupname\n";
                } else {
                    $line .= "${SPACE}" ."  real(RLEN),dimension(:),pointer ::$pre$groupname\n";
                }
                $line .= "${SPACE}" ."  integer, intent(IN)            :: n\n";
                $line .= "${SPACE}" ."  integer, intent(IN)            :: constituent\n";
                $line .= "${SPACE}" ."  integer, intent(IN), optional  :: cmax\n";
                $line .= "${SPACE}" ."  \n";
                if ( $pre eq "pp" ) {
                    my @refers    = ();
                    my @refersMax = ();
                    foreach my $member (@members){
                        my @components_mem = ( sort { ($$LST_CONST{$a} cmp $$LST_CONST{$b}) } keys %{$member->getComponents()} );
                        my $maxkey = $$LST_CONST{ $components_mem[ $#components_mem ] }; #insert the minimun component
                        my $minkey = $components_mem[0];                                 #insert the index of the maximun component

                        push( @refers, "pp" . $member->getSigla() . ${minkey} );
                        push( @refersMax, ${maxkey} );
                    }
                    if( scalar(@refers)    == 0 ){ push(@refers,    0); }
                    if( scalar(@refersMax) == 0 ){ push(@refersMax, 0); }
                    $line .= "${SPACE}" ."  integer,dimension(" . max(1,scalar(@members))      . ") :: referto=(/"         . join(',',@refers)    . "/)\n";
                    $line .= "${SPACE}" ."  integer,dimension(" . max(1,scalar(@members))      . ") :: const_max=(/"       . join(',',@refersMax) . "/)\n";

                    my @consAdd = (0) x scalar(keys %$LST_CONST);
                    my $num_ele = 0;
                    
                    foreach my $elem ( sort { ($$LST_CONST{$a} cmp $$LST_CONST{$b}) } keys %{$group->getComponents()} ){
                        #print "$elem [". ($$LST_CONST{$elem}-1) ."]: $num_ele ";
                        $consAdd[($$LST_CONST{$elem}-1)] = $num_ele++;
                    }
                    $line .= "${SPACE}" ."  integer,dimension(" . scalar(keys %$LST_CONST) . ") :: constituent_add=(/" . join(',',@consAdd)   . "/)\n";
                    $line .= "${SPACE}" ."\n";
                    $line .= "${SPACE}" ."  if ( constituent <=const_max(n) ) then\n";
                    $line .= "${SPACE}" ."   $pre$groupname=referto(n)+ constituent_add(constituent)\n";
                    $line .= "${SPACE}" ."  else\n";
                    $line .= "${SPACE}" ."   $pre$groupname=0\n";
                    $line .= "${SPACE}" ."  endif\n";
                    $line .= "${SPACE}" ."  if ( present(cmax) ) $pre$groupname = const_max(n)\n";
                }else{
                    $line .= "${SPACE}" ."  $pre$groupname => D${dim}STATE(pp${groupname}(n,constituent),:)\n";
                }
                $line .= "${SPACE}" ."\n";
                $line .= "${SPACE}" ."END function\n";
                print $file $line;
            }

        }
    }
}

###########################################  SET_VAR_INFO_BFM ######################


sub func_STRING  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING: $dim, $type, $subt"; }

    my $line  = "";

    if( ! defined $STRING_INDEX_ARRAY{"${type}_${subt}_${dim}_S"} ){
        $STRING_INDEX_ARRAY{"${type}_${subt}_${dim}_S"} = $STRING_INDEX;
    }

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType() 
            && $subt eq $param->getSubtype() ){

            my $comm  = $param->getComment();
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        my $unitC = ${$param->getComponents()}{$const};
                        $line .= "${SPACE}var_names($STRING_INDEX)=\"$nameC\"\n";
                        $line .= "${SPACE}var_long($STRING_INDEX)=\"$comm\"\n";
                        $line .= "${SPACE}var_units($STRING_INDEX)=\"$unitC\"\n";
                        if( $type eq 'diagnos' || $type eq 'state'){ $$LST_INDEX{"$nameC"} = "$STRING_INDEX"; }
                        $STRING_INDEX++;
                    }
                }
            }else{
                my $unit  = $param->getUnit();
                $line .= "${SPACE}var_names($STRING_INDEX)=\"$name\"\n";
                $line .= "${SPACE}var_long($STRING_INDEX)=\"$comm\"\n";
                $line .= "${SPACE}var_units($STRING_INDEX)=\"$unit\"\n";
                if( $type eq 'diagnos' || $type eq 'state'){ $$LST_INDEX{"$name"} = "$STRING_INDEX"; }
                $STRING_INDEX++;
            }
        }
    }

    $STRING_INDEX_ARRAY{"${type}_${subt}_${dim}_E"} = $STRING_INDEX - 1;

    if( $line ){ print $file $line; }
}


sub func_STRING_DIAGG  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_DIAGG: $dim, $type, $subt"; }

    my $line  = "";
    my $index = 1;

    if( ! defined $STRING_INDEX_ARRAY{"${type}_${subt}_${dim}_S"} ){
        $STRING_INDEX_ARRAY{"${type}_${subt}_${dim}_S"} = $STRING_INDEX;
    }

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $param->getDim() == $dim 
            && $param->getType eq $type 
            && $param->getSubtype eq $subt ){

            my $group_name = $param->getQuota();
            my $unit       = $param->getUnit();
            my $comm       = $param->getComment();
            foreach my $name2 ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                my $param2 = $$LST_PARAM{$name2};
                if ( $param2->getGroup() && $group_name eq $param2->getGroup() ){ 
                    my $param_name = $param2->getSigla();
                    my $nameP = $name . "(ii$param_name)";
                    my $commP = $comm;$commP =~ s/$group_name/$param_name\($group_name\)/;

                    $line .= "${SPACE}var_names($STRING_INDEX)=\"$nameP\"\n";
                    $line .= "${SPACE}var_long($STRING_INDEX)=\"$commP\"\n";
                    $line .= "${SPACE}var_units($STRING_INDEX)=\"$unit\"\n";
                    $STRING_INDEX++;
                }
            }
        }
    }

    $STRING_INDEX_ARRAY{"${type}_${subt}_${dim}_E"} = $STRING_INDEX - 1;

    if( $line ){ print $file $line; }
}


sub func_STRING_FIELD  {
    my ( $file, $subt, $spec ) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_FIELD: $subt, $spec"; }

    my $line;
    
    my $SPEC       = uc($spec);
    my $spec_short = substr($spec,0,3);

    if( ! defined $STRING_INDEX_ARRAY{"diagnos_${subt}_${spec}_2_S"} ){
        $STRING_INDEX_ARRAY{"diagnos_${subt}_${spec}_2_S"} = $STRING_INDEX;
    }

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( 3 == $param->getDim() 
            && 'state' eq $param->getType() 
            && $subt eq $param->getSubtype() ){

            my $comm  = $param->getComment();
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        my $unitC = decreaseUnit( ${$param->getComponents()}{$const} );
                        $line .= "${SPACE}var_names($STRING_INDEX)=\"j${spec_short}${nameC}\"\n";
                        $line .= "${SPACE}var_long($STRING_INDEX)=\"flux of $comm at $SPEC\"\n";
                        $line .= "${SPACE}var_units($STRING_INDEX)=\"$unitC\"\n";
                        $STRING_INDEX++;
                    }
                }
            }else{
                my $unit  = decreaseUnit($param->getUnit());
                $line .= "${SPACE}var_names($STRING_INDEX)=\"j${spec_short}${name}\"\n";
                $line .= "${SPACE}var_long($STRING_INDEX)=\"flux of $comm at $SPEC\"\n";
                $line .= "${SPACE}var_units($STRING_INDEX)=\"$unit\"\n";
                $STRING_INDEX++;
            }
        }
    }

    $STRING_INDEX_ARRAY{"diagnos_${subt}_${spec}_2_E"} = $STRING_INDEX - 1;

    if( $line ){ print $file $line; }
}


sub func_STRING_INDEX  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_INDEX:"; }

    my $line;

    my $search_S = "";
    my $search_E = "";
    my $out      = "";

    my $subt_short = ucfirst(substr($subt,0,3));

    my $type_short = ucfirst(substr($type,0,4));
    if( $type_short eq "Stat" ){ $type_short = 'State'; } #fix because "state" has 5 chars
    if( $dim == 2 ){ $type_short .= '2d'; } #fix because name 2D not collapse with 3D
    $search_S = "${type}_${subt}_${dim}_S";
    $search_E = "${type}_${subt}_${dim}_E";
    $out = "${type_short}";

    if( $type eq "diagnos" ){
        #fix because "diagnos" ends with "diaggrp" values
        $search_E = "diaggrp_${subt}_${dim}_E";
    }
    if( $subt eq 'pel' && $type eq 'state' && $dim == 2 ){
        #fix because "Pel state 2D" ends with "Pel diagnos river 2D" values
        $search_E = "diagnos_${subt}_river_2_E";
    }

    if( defined $STRING_INDEX_ARRAY{"${search_S}"} && defined $STRING_INDEX_ARRAY{"${search_E}"} ){
        $line .= "${SPACE}st${subt_short}${out}S=" . $STRING_INDEX_ARRAY{"${search_S}"} . "\n"; 
        $line .= "${SPACE}st${subt_short}${out}E=" . $STRING_INDEX_ARRAY{"${search_E}"} . "\n";
    }

    if( $line ){ print $file $line; }
}

sub func_STRING_INDEX_FIELD  {
    my ( $file, $subt, $spec ) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_INDEX_FIELD:"; }

    my $line;

    my $search_S = "";
    my $search_E = "";
    my $out      = "";

    my $subt_short = ucfirst(substr($subt,0,3));

    my $spec_short = ucfirst(substr($spec,0,3));
    my $dim = 2; #fix beacuse "specs" are 2d

    $search_S = "diagnos_${subt}_${spec}_${dim}_S";
    $search_E = "diagnos_${subt}_${spec}_${dim}_E";
    $out = "${spec_short}";

    if( $spec eq "surface" ){
        #fix because "surface" start with "diagnos" values 
        $search_S = "diagnos_${subt}_${dim}_S";
    }    

    if( defined $STRING_INDEX_ARRAY{"${search_S}"} && defined $STRING_INDEX_ARRAY{"${search_E}"} ){
        $line .= "${SPACE}st${subt_short}${out}S=" . $STRING_INDEX_ARRAY{"${search_S}"} . "\n"; 
        $line .= "${SPACE}st${subt_short}${out}E=" . $STRING_INDEX_ARRAY{"${search_E}"} . "\n";
    }

    if( $line ){ print $file $line; }
}


###########################################  INIT_VAR_BFM ######################


sub func_INIT_VALUE_CALC {
    my ( $file, $subt, $value ) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_VALUE_CALC: "; }

    my $line = '';

    if ( !$value ){ $value = 0; }

    if($subt eq 'pel' ){ $line = "${SPACE}CalcPelagicFlag = " . ${value} . "\n"; }
    else               { $line = "${SPACE}CalcBenthicFlag = " . ${value} . "\n"; }


    foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$groupname};
        if( $subt eq $group->getSubtype()){

            $line .= "${SPACE}" . join(',', 'Calc' . $group->getSigla()) . " = " . ${value} . "\n";
        }
    }

    if( $line ){ print $file $line; }
}


sub func_INIT_PP  {
    my ( $file, $before, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_PP: "; }

    my @line_par = ();

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const . '0';
                        push( @line_par, ${nameC} );
                    }
                }
            }else{
                push( @line_par, (${name} . '0') );
            }
        }
    }

    if( $#line_par >= 0 ){ 
        printList($file, \@line_par, $before ); print $file "\n";
    }
}

sub func_INIT_FUNC_CONST {
    my ( $file ) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_FUNC_CONST: "; }

    my $line = '';
    my @constList = ();
    my @constNoCList = ();
    my @constOptionalList = ();
    my @constOptionalRatioList = ();

    foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
        push(@constList, $const);
        if($const ne $constList[0]){
            push( @constNoCList, $const );
            push( @constOptionalList, $const . $constList[0] );
            push( @constOptionalRatioList, $const . $constList[0] . '_ratio' );
        }
    }

    $line .="   subroutine init_constituents( ". join(',',@constList)  . "," . join( ',', @constOptionalList )    . ")\n";
    $line .="     use global_mem, only: RLEN,ZERO\n";
    $line .="     IMPLICIT NONE\n";
    $line .="     real(RLEN),dimension(:),intent(in)             :: " . $constList[0]                              . "\n";
    $line .="     real(RLEN),intent(in),optional                 :: " . join( ',', @constOptionalList )            . "\n";
    $line .="     real(RLEN),dimension(:),intent(inout),optional :: " . join( ',', @constList[1 .. $#constList]  ) . "\n";
    $line .="     real(RLEN)                                     :: " . join( ',', @constOptionalRatioList       ) . "\n";
    $line .="     \n";

    foreach my $constOpt (@constOptionalList){
        $line .= "     " . "    " . $constOpt . "_ratio = " . $constOpt . "_ratio_default\n";
        $line .= "     " . "    if (present(" . $constOpt . ")) then\n";
        $line .= "     " . "      if (" . $constOpt . ">ZERO) " . $constOpt . "_ratio = " . $constOpt . "\n";
        $line .= "     " . "    end if\n\n";
    }

    foreach my $constNoC (@constNoCList){
        $line .= "     " . "    if (present(" . $constNoC . ")) then\n";
        $line .= "     " . "      where (" . $constNoC . "==ZERO)\n";
        $line .= "     " . "        " . $constNoC . " = " . $constNoC . "c_ratio*c\n";
        $line .= "     " . "      end where\n";
        $line .= "     " . "    end if\n";
    }

    $line .="   end subroutine init_constituents\n";


    if( $line ){ print $file $line; }
}



sub func_INIT_FUNC_ZERO {
    my ( $file, $dim, $subt ) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_FUNC_ZERO: "; }
    my $line = '';


    foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$groupname};

        if( $group->getDim() == $dim
            && $subt eq $group->getSubtype()){

            $line .= "${SPACE}do j = 1, ii" . $group->getSigla()    . "\n";
            $line .= "${SPACE}  if (.NOT.Calc" . $group->getSigla() ."(j)) then\n";
            $line .= "${SPACE}    iiLastElement=pp" . $group->getSigla() . "(j,1,cmax=1)\n";
            $line .= "${SPACE}    do i = 1,iiLastElement\n";
            $line .= "${SPACE}      D3STATE(pp" . $group->getSigla() . "(j,i),:) = p_small\n";
            $line .= "${SPACE}      D3STATETYPE(pp" . $group->getSigla() . "(j,i)) = OFF\n";
            $line .= "#if defined key_obcbfm\n";
            $line .= "${SPACE}      D3STATEOBC(pp" . $group->getSigla() . "(j,i)) = NOOBCSTATES\n";
            $line .= "#endif\n";
            $line .= "${SPACE}    end do\n";
            $line .= "${SPACE}  end if\n";
            $line .= "${SPACE}end do\n\n";
        }
    }

    if( $line ){ print $file $line; }
}

sub func_INIT_DEFAULT  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_DEFAULT: "; }

    my $line = '';

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        $line .= "${SPACE}" . $name . $const . "0 = _ZERO_\n";
                    }
                }
            }else{
                $line .= "${SPACE}" . $name . "0 = _ZERO_\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_INIT_SETS  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INITS_SETS: "; }

    my $line = '';

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        $line .= "${SPACE}  " . $name . $const . " = " . $name . $const . "0\n";
                    }
                }
            }else{
                $line .= "${SPACE}  " . $name . " = " . $name . "0\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_INIT_INTERNAL {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_INTERNAL: "; }

    my $line = '';
    my @constList         = ();
    my @constNoC          = ();
    my @constOptionalList = ();

    foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
        push(@constList, $const);
        if($const ne $constList[0]){ 
            push( @constNoC         , $const                 );
            push( @constOptionalList, $const . $constList[0] ); 
        }
    }


    foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group     = $$LST_GROUP{$groupname};
        my $groupAcro = $$LST_GROUP{$groupname}->getAcro();

        if( $group->getDim() == $dim
            && $subt eq $group->getSubtype()){

            $line .= "${SPACE}do i = 1 , ( ii". $groupname ." )\n";
            $line .= "${SPACE}  call init_constituents( c=" . $groupname  . "(i,iiC), &\n${SPACE}";
            my @temp_line;
            foreach my $const (@constNoC){
                if( exists ${$$LST_GROUP{$groupname}->getComponents()}{$const} ){
                    push( @temp_line, "    " . $const . "=D3STATE(pp" . $groupname . "(i,ii" . uc($const) . "),:)" );
                }
            }
            my @temp_line2;
            foreach my $constOpt (@constOptionalList){
                #get the constituents active inside the group
                if( exists ${$$LST_GROUP{$groupname}->getComponents()}{ (split('',$constOpt))[0] } ){
                    my $temp_compo = "p_q" . $constOpt . $groupAcro;
                    #if the optional initialization element exists in the namelist => add to initialize constituents
                    foreach my $list (@$LST_NML){
                        foreach my $param ( @{$list->slots()} ){
                            if( $param eq $temp_compo){ 
                                push( @temp_line2, $constOpt . "=" . $temp_compo . "(i)" );
                            }
                        }
                    }
                }
            }
            if( scalar(@temp_line2) ){ push( @temp_line, "    " . join(",  ", @temp_line2 ) ); }

            $line .= join(",  &\n${SPACE}", @temp_line );
            $line .= " )\n";
            $line .= "${SPACE}end do\n";
        }
    }

    if( $line ){ print $file $line; }
}



###########################################  INCLUDE ######################



sub func_HEADER  {
    my ( $file, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_HEADER: "; }

    my $line = '';
    my $TYPE = uc($type);

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        $line .= "#define " . $nameC . "(A) D" . $dim . $TYPE . "(pp" . $nameC . ",A)\n";
                    }
                }
            }else{
                $line .= "#define " . $name . "(A) D" . $dim . $TYPE . "(pp" . $name . ",A)\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_GROUP_HEADER  {
    my ( $file, $dim, $subt ) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_GROUP_HEADER: "; }

    my $line = '';

    foreach my $name ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$name};
        if( $dim == $group->getDim() 
            && $subt eq $group->getSubtype()){
            $line .= "#define " . $name . "(A,B) D" . $dim . "STATE(pp" . $name . "(A,B),:)\n";
        }
    }

    if( $line ){ print $file $line; }
}


sub func_HEADER_FIELD  {
    my ( $file, $subt, $spec) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_HEADER_FIELD: "; }

    my $line_ini = "";
    my $line_par = "";
    
    my $SPEC       = uc($spec);
    my $spec_short = substr($spec,0,3);

    my $n = 0;
    if( exists $$LST_STA{"state 2d $subt ${spec} index"} ){
        $n = $$LST_STA{"state 2d $subt ${spec} index"};
    }


    $line_ini .= "#define PEL${SPEC}(A,B) D2DIAGNOS($n+A,B)\n";
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( 3 == $param->getDim() 
            && 'state' eq $param->getType()
            && $subt eq $param->getSubtype() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        $line_par .= "#define j${spec_short}${nameC}(B) D2DIAGNOS($n+pp$nameC,B)\n";
                    }
                }
            }else{
                $line_par .= "#define j${spec_short}${name}(B) D2DIAGNOS($n+pp${name},B)\n";
            }
        }
    }

    if( $line_par ){ print $file $line_ini . $line_par; }
}



sub func_HEADER_DIAGG  {
    my ( $file, $dim, $type, $subt) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_HEADER_DIAGG: "; }

    my $line = '';
    my $TYPE = uc($type);

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){
            $line .= "#define " . $name . "(A,B) D" . $dim . "DIAGNOS(pp" . $name . "(A),B)\n";
        }
    }

    if( $line ){ print $file $line; }
}
