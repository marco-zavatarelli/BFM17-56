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
my $DEBUG = 0;
my $SPACE = "   ";
my $MAX_CHAR_PER_LINE = 76;
my $STRING_INDEX=1;
my %STRING_INDEX_ARRAY = ();
my ( $LST_PARAM, $LST_GROUP, $LST_STA, $LST_CONST, $LST_INDEX, $LST_NML );
########### VARIABLES GLOBAL ##########################





########### REGULAR EXPRESSIONS ##########################
my $dispatch = {
    #COMMENT LINE FOR DEVELOPERS
    '\%\!.*' => \&func_COMMENT,
    #DESCRIPTION
    '\%(3|2)d-(diagnos|state)-(pel|ben|ice)-desc(?:\s|\n)'=> \&func_DESC ,
    '\%(3|2)d-(diaggrp)-(pel|ben|ice)-desc(?:\s|\n)' => \&func_DESC_DIAGG ,
    #NR
    '([^\%]*)\%(3|2)d-(flux|diagnos|state)-(pel|ben|ice)-nr(?:\s|\n)'=> \&func_NR ,
    #ARRAY
    '([^\%]*)\%(pel|ben|ice)-field-array\s*(\S*)(?:\s|\n)'=> \&func_ARRAY_FIELD ,
    #PP
    '([^\%]*)\%(3|2)d-(diagnos|state)-(pel|ben|ice)-pp(?:\s|\n)'=> \&func_PP ,
    '([^\%]*)\%(3|2)d-(diaggrp)-(pel|ben|ice)-pp(?:\s|\n)'=>       \&func_PP_DIAGG ,
    '\%(3|2)d-(diaggrp)-(pel|ben|ice)-assign-pp(?:\s|\n)'=>        \&func_PP_ASSIGN ,
    #POINTER
    '\%(3|2)d-(diagnos|state)-(pel|ben|ice)-pointer(?:\s|\n)'=>    \&func_POINT ,
    '([^\%]*)\%(3|2)d-(diaggrp)-(pel|ben|ice)-pointer(?:\s|\n)'=>    \&func_POINT_DIAGG ,
    '([^\%]*)\%(pel|ben|ice)-field-pointer\s*(\S*)(?:\s|\n)'=>    \&func_POINT_FIELD ,
    '\%(3|2)d-(diagnos|state)-(pel|ben|ice)-alloc-pointer(?:\s|\n)' =>    \&func_POINT_ALLOC ,
    '\%(3|2)d-(diaggrp)-(pel|ben|ice)-alloc-pointer(?:\s|\n)' =>    \&func_POINT_ALLOC_DIAGG ,
    '\%(pel|ben|ice)-alloc-pointer-field\s*(\S*)(?:\s|\n)' =>    \&func_POINT_ALLOC_FIELD ,
    #HEADER
    '\%(3|2)d-(diagnos|state)-(pel|ben|ice)-header(?:\s|\n)' =>    \&func_HEADER ,
    '\%(3|2)d-(diaggrp)-(pel|ben|ice)-header(?:\s|\n)' =>    \&func_HEADER_DIAGG ,
    '\%(pel|ben|ice)-field-header\s*(\S*)(?:\s|\n)' =>    \&func_HEADER_FIELD ,
    #STRING
    '\%(3|2)d-(diagnos|state|flux)-(pel|ben|ice)-string(?:\s|\n)'       => \&func_STRING ,
    '\%(3|2)d-(diaggrp)-(pel|ben|ice)-string(?:\s|\n)'                  => \&func_STRING_DIAGG ,
    '\%(pel|ben|ice)-field-string\s*(\S*)(?:\s|\n)'                     => \&func_STRING_FIELD ,
    '\%(3|2)d-(diagnos|state|flux)-(pel|ben|ice)-string-index(?:\s|\n)' => \&func_STRING_INDEX ,
    '\%(pel|ben|ice)-string-index-field\s*(\S*)(?:\s|\n)'               => \&func_STRING_INDEX_FIELD ,
    '\%startend-string-index\s*(\S*)(?:\s|\n)'                          => \&func_STRING_INDEX_STARTEND ,
    #ALLOC
    '\%(3|2)d-(diagnos|state)-(pel|ben|ice)-alloc(?:\s|\n)' => \&func_ALLOC ,
    '\%(1|3|2)d-(intvar|variable|variables)-(pel|ben|ice)-alloc(?:\s|\n)'=>    \&func_ALLOC_INTVAR ,
    #FLUX
    '\%(3|2)d-(flux)-(pel|ben|ice)-alloc(?:\s|\n)' => \&func_FLUX_ALLOC ,
    '\%(3|2)d-(flux)-(pel|ben|ice)-fill(?:\s|\n)' => \&func_FLUX_FILL ,
    #CONSTITUENT
    '([^\%]*)\%constituent(?:\s|\n)' => \&func_CONSTITUENT ,
    #GROUP    
    '\%(3|2)d-group-(pel|ben|ice)-header(?:\s|\n)'                => \&func_GROUP_HEADER ,
    '([^\%]*)\%(3|2)d-group-(pel|ben|ice)-parameter(?:\s|\n)'     => \&func_GROUP_PARAMETER ,
    '([^\%]*)\%group-(pel|ben|ice)-calc(?:\s|\n)'                 => \&func_GROUP_CALC,
    '([^\%]*)\%(3|2)d-group-(pel|ben|ice)-function-name(?:\s|\n)' => \&func_GROUP_FUNCTION_NAME ,
    '\%(3|2)d-groupfunctions-(pel|ben|ice)(?:\s|\n)'              => \&func_GROUP_FUNCTIONS ,
    #INTVAR/VARIABLE
    '([^\%]*)\%(3|2|1)d-(intvar|variable)-(pel|ben|ice)(?:\s|\n)' => \&func_VARIABLE ,
    #INIT
    '([^\%]*)\%value-init-calc-(pel|ben|ice)\s*(\S*)(?:\s|\n)' => \&func_INIT_VALUE_CALC ,
    '([^\%]*)\%(3|2)d-(state)-(pel|ben|ice)-Initpp(?:\s|\n)'   => \&func_INIT_PP ,
    '\%init-ratiodefault-constituents(?:\s|\n)'                => \&func_INIT_RATIO_CONST ,
    '\%init-func-constituents(?:\s|\n)'                        => \&func_INIT_FUNC_CONST ,
    '\%init-(pel|ben|ice)-namelist\s*(\d*)\s*(\d*)(?:\s|\n)'   => \&func_INIT_NAMELIST ,
    '\%(3|2)d-init-(pel|ben|ice)-output-variables(?:\s|\n)'    => \&func_INIT_OUTPUT_VARIABLES ,
    '\%(3|2)d-(state)-(pel|ben|ice)-InitDefault(?:\s|\n)'      => \&func_INIT_DEFAULT ,
    '\%(3|2)d-(state)-(pel|ben|ice)-InitSets(?:\s|\n)'         => \&func_INIT_SETS ,
    '\%(3|2)d-(state)-(pel|ben|ice)-InitInternal(?:\s|\n)'     => \&func_INIT_INTERNAL ,
    '\%(3|2)d-state-(pel|ben|ice)-func-zeroing(?:\s|\n)'       => \&func_INIT_FUNC_ZERO ,
};
########### REGULAR EXPRESSIONS ##########################


########### FUNCTIONS ##########################

sub print_f90 {

    my ($templ_mem, $output, $lst_group, $lst_param, $lst_sta, $lst_const, $lst_index, $lst_nml, $debug) = @_;

    $LST_PARAM = $lst_param;
    $LST_GROUP = $lst_group;
    $LST_STA   = $lst_sta;
    $LST_CONST = $lst_const;
    $LST_INDEX = $lst_index;
    $LST_NML = $lst_nml;
    $DEBUG   = $debug;

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
                if( $DEBUG ){ print " $line_raw"; }
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
    my ( $dim, $subt ) = @_;

    my $size = 0;
    foreach my $root (keys %$LST_PARAM){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() && $subt eq $param->getSubtype() && $param->getQuota() ){
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
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_ALLOC: "; }

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type, $subt)  ){ return; }
    
    my $j = '';
    my $TYPE = uc($type);
    if ( $dim == 2 ){ $j = '_XY'; }
    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    print $file "${SPACE} \n";
    print $file "${SPACE}allocate(D${dim}${TYPE}${SUBTYPE}(1:NO_D${dim}_BOX_${TYPE}S${SUBTYPE},1:NO_BOXES${j}),stat=status)\n";
    print $file "${SPACE}if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}${TYPE}${SUBTYPE}\")\n";
    print $file "${SPACE}D${dim}${TYPE}${SUBTYPE} = ZERO\n";
    if ( $type eq "state" ) {
        print $file "#ifndef EXPLICIT_SINK\n";
        print $file "${SPACE}  allocate(D${dim}SOURCE${SUBTYPE}(1:NO_D${dim}_BOX_STATES${SUBTYPE},1:NO_BOXES${j}),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SOURCE${SUBTYPE}\")\n";
        print $file "${SPACE}  D${dim}SOURCE${SUBTYPE} = ZERO\n";
        print $file "#else\n";
        print $file "${SPACE}  allocate(D${dim}SOURCE${SUBTYPE}(1:NO_D${dim}_BOX_STATES${SUBTYPE},1:NO_D${dim}_BOX_STATES${SUBTYPE},1:NO_BOXES${j}),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SOURCE${SUBTYPE}\")\n";
        print $file "${SPACE}  D${dim}SOURCE${SUBTYPE} = ZERO\n";
        print $file "${SPACE}  allocate(D${dim}SINK${SUBTYPE}(1:NO_D${dim}_BOX_STATES${SUBTYPE},1:NO_D${dim}_BOX_STATES${SUBTYPE},1:NO_BOXES${j}) ,stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SINK${SUBTYPE}\")\n";
        print $file "${SPACE}  D${dim}SINK${SUBTYPE} = ZERO\n";
        print $file "#endif\n";
        print $file "${SPACE}allocate(D${dim}STATETYPE${SUBTYPE}(1:NO_D${dim}_BOX_${TYPE}S${SUBTYPE} ),stat=status)\n";
        print $file "${SPACE}if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\",\"D${dim}STATETYPE${SUBTYPE}\")\n";
        print $file "${SPACE}D${dim}STATETYPE${SUBTYPE} = ZERO\n";
        print $file "#ifdef BFM_NEMO\n";
        print $file "${SPACE}  allocate(D${dim}STATEOBC${SUBTYPE}(1:NO_D${dim}_BOX_STATES${SUBTYPE}),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\",\"D${dim}STATEOBC${SUBTYPE}\")\n";
        print $file "${SPACE}  D${dim}STATEOBC${SUBTYPE} = ZERO\n";
        print $file "#endif\n";
    }
}

sub func_POINT_ALLOC  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_POINT_ALLOC: "; }

    my $line = "";
    my $TYPE = uc($type);
    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $param->getSubtype eq $subt ){

            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        $line .= "$nameC => D${dim}${TYPE}${SUBTYPE}(pp$nameC,:); $nameC=ZERO\n";
                    }
                }
            }else{
                $line .= "$name => D${dim}${TYPE}${SUBTYPE}(pp$name,:); $name=ZERO\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}

sub func_PP_ASSIGN  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_PP_ASSIGN: "; }
    
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
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_PP_ALLOC_DIAG: "; }
    
    my $line = "";
    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

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
                $line .= "$root => D${dim}DIAGNOS${SUBTYPE}(pp${root}(ii${mem1}): pp${root}(ii${mem2}),:)\n";
                $line .= "$root=ZERO\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_POINT_ALLOC_FIELD  {
    my ( $file, $subt, $spec) = @_;
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_POINT_ALLOC_FIELD: "; }

    if( $subt ne 'pel' ){ print "ERROR: Subtype can only be \'pel\' for field $spec\n"; exit; }

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
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_ALLOC_INTVAR: "; }
    
    my $line = "";
    my $j = "";
    my $TYPE = uc($type);
    if ( $dim == 2 ){ $j = "_XY"; }

    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

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
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_FLUX_ALLOC: "; }

    my $line = "";
   
    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type, $subt)  ){ return; }

    if( $dim == 3 ){
        $line .= "  allocate( D${dim}FLUX_MATRIX${SUBTYPE}(1:NO_D${dim}_BOX_STATES${SUBTYPE}, 1:NO_D${dim}_BOX_STATES${SUBTYPE}),stat=status )\n";
        $line .= "  allocate( D${dim}FLUX_FUNC${SUBTYPE}(1:NO_D${dim}_BOX_FLUX${SUBTYPE}, 1:NO_BOXES),stat=status )\n";
        $line .= "  D${dim}FLUX_FUNC${SUBTYPE} = 0\n";
    }elsif( $dim == 2 ){
        $line .= "  allocate( D${dim}FLUX_MATRIX${SUBTYPE}(1:NO_D${dim}_BOX_STATES${SUBTYPE}, 1:NO_D${dim}_BOX_STATES${SUBTYPE}),stat=status )\n";
        $line .= "  allocate( D${dim}FLUX_FUNC${SUBTYPE}(1:NO_D${dim}_BOX_FLUX${SUBTYPE}),stat=status )\n";
        $line .= "  D${dim}FLUX_FUNC${SUBTYPE} = 0\n";
    }

    if( $line ){ print $file $line; }
}


sub func_FLUX_FILL  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $DEBUG ){ print "AllocateMem -> FUNCTION CALLED func_FLUX_ALLOC: "; }
    my $line = "";

    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

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
                        push( @{$d3flux_func{$name}}, "D${dim}FLUX_MATRIX${SUBTYPE}(pp$compo1, pp$compo2)%p(" . $d3flux_matrix_index{"pp$compo1, pp$compo2"} . ")=${sign1}${flxindex}\n " );
                    }else{ # B-> A
                        if ( exists $d3flux_matrix_index{"pp$compo2, pp$compo1"} ){ 
                            $d3flux_matrix_index{"pp$compo2, pp$compo1"}++ } 
                        else{ 
                            $d3flux_matrix_index{"pp$compo2, pp$compo1"} = 1
                        };
                        push( @{$d3flux_func{$name}}, "D${dim}FLUX_MATRIX${SUBTYPE}(pp$compo2, pp$compo1)%p(" . $d3flux_matrix_index{"pp$compo2, pp$compo1"} . ")=${sign1}${flxindex}\n " );
                    }

                    if( $compo1 eq $compo2 ){ # A == B
                        if ( exists $d3flux_matrix_index_dir{"pp$compo1, pp$compo2"} ){ 
                            $d3flux_matrix_index_dir{"pp$compo1, pp$compo2"}++ } 
                        else{ 
                            $d3flux_matrix_index_dir{"pp$compo1, pp$compo2"} = 1  
                        }
                        push( @{$d3flux_func_dir{$name}}, "D${dim}FLUX_MATRIX${SUBTYPE}(pp$compo1, pp$compo2)%dir(" . $d3flux_matrix_index{"pp$compo1, pp$compo2"} . ")=${dir}\n " );
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
            $line .= "  allocate( D${dim}FLUX_MATRIX${SUBTYPE}($key)%p( 1:$d3flux_matrix_index{$key} ) )\n";
            if( exists $d3flux_matrix_index_dir{$key} ){
                $line .= "  allocate( D${dim}FLUX_MATRIX${SUBTYPE}($key)%dir( 1:$d3flux_matrix_index_dir{$key} ) )\n";
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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_DESC: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_ARRAY_FIELD: "; }

    #save the spec in stadistics for further calculations of number of variables
    if( ! exists $$LST_STA{"spec num $subt"} ){ $$LST_STA{"spec num $subt"} = 1; }
    else{ $$LST_STA{"spec num $subt"} += 1; }

    my $spec_upper = uc($spec);
    my $subt_short = uc(substr($subt,0,3));
    
    print $file "$before ${subt_short}${spec_upper}\n";
}

sub func_NR  {
    my ( $file, $before, $dim, $type, $subt) = @_;
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_NR: "; }

    my $title  = "${type} ${dim}d ${subt}";
    my $number = 0;
    if( exists $$LST_STA{$title} ){
        if( $title eq "diagnos 2d ${subt}" ){
            #print "TITLE(${title}):" . $$LST_STA{"${title}"} . " + sizeGroup($dim, $subt): " . sizeGroup($dim, $subt) 
            #    . " + ( state-3d-${subt}: " . $$LST_STA{"state 3d ${subt}"} . " * spec-num-${subt}: " . $$LST_STA{"spec num ${subt}"} . ")\n";
            if( exists($$LST_STA{"spec num ${subt}"}) ){
                $number = $$LST_STA{"${title}"} + sizeGroup($dim, $subt) + ( $$LST_STA{"state 3d ${subt}"} * $$LST_STA{"spec num ${subt}"} );
            }else{
                $number = $$LST_STA{"${title}"} + sizeGroup($dim, $subt);
            }
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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_PP: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_POINT: "; }
    
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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_CONSTITUENT: "; }
    
    my @line_par = ();

    foreach my $key ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
        push( @line_par, "ii" . uc($key) . "=" . $$LST_CONST{$key} );
    }

    if( $#line_par > 0 ){ 
        #Create global constant for the last constituent index
        push( @line_par, "iiLastElement=" . scalar( keys %$LST_CONST ) );
        printList($file, \@line_par, $before); print $file "\n";
    }
}


sub func_GROUP_PARAMETER  {
    my ( $file, $before, $dim, $subt ) = @_;
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_GROUP_PARAMETER: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_GROUP_CALC: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_VARIABLE: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_DESC_DIAGG: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_PP_DIAGG: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_POINT_DIAGG: "; }
    
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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_POINT_FIELD: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_GROUP_FUNCTION_NAME: "; }

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
    if ( $DEBUG ){ print "ModuleMem -> FUNCTION CALLED func_GROUPFUNCTIONS: "; }

    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

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
                $line .= "${SPACE}" ."function $pre${groupname}(n,constituent)\n";
                $line .= "${SPACE}" ."\n";
                $line .= "${SPACE}" ."  IMPLICIT NONE\n";
                $line .= "${SPACE}" ."\n";
                if ( $pre eq "pp" ) {
                    $line .= "${SPACE}" ."  integer ::$pre$groupname\n";
                } else {
                    $line .= "${SPACE}" ."  real(RLEN),dimension(:),pointer ::$pre$groupname\n";
                }
                $line .= "${SPACE}" ."  integer, intent(IN)            :: n\n";
                $line .= "${SPACE}" ."  integer, intent(IN)            :: constituent\n";
                $line .= "${SPACE}" ."  \n";
                if ( $pre eq "pp" ) {

                    my @line_pointers = ();
                    foreach my $member (@members){
                        my @refers    = ();
                        my %const_mem = %{$member->getComponents()};

                        foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                            if( exists $const_mem{$const}){
                                push( @refers, "pp" . $member->getSigla() . ${const} . " ");
                            }else{
                                push( @refers, "0     " );
                            }
                        }
                        push( @line_pointers, join(",",@refers) );
                    }

                    $line .= "${SPACE}" . "  integer,dimension(" . max(1,scalar(@members)) . " * iiLastElement ) :: pointers = (/ & \n";
                    $line .= "${SPACE}" .    "${SPACE}       " 
                       . join(", &\n${SPACE}       ", @line_pointers) . "  &\n"
                       .      "${SPACE}       /)\n\n";

                    $line .= "${SPACE}" ."  IF( n > " . max(1,scalar(@members)) . " .OR. n == 0 ) THEN\n";                    
                    $line .= "${SPACE}" ."    pp" . ${groupname}  . " = 0\n";
                    $line .= "${SPACE}" ."  ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN\n";
                    $line .= "${SPACE}" ."    pp" . ${groupname}  . " = 0\n";
                    $line .= "${SPACE}" ."  ELSE\n";
                    $line .= "${SPACE}" ."    pp" . ${groupname}  . " = pointers( ( (n-1) * iiLastElement ) + constituent )\n";
                    $line .= "${SPACE}" ."  ENDIF\n";

                }else{
                    $line .= "${SPACE}" ."  if ( pp" . ${groupname} . "(n,constituent) > 0 ) then\n";
                    $line .= "${SPACE}" ."    $pre$groupname => D${dim}STATE${SUBTYPE}(pp" . ${groupname} . "(n,constituent),:)\n";
                    $line .= "${SPACE}" ."  else\n";
                    $line .= "${SPACE}" ."    $pre$groupname => NULL()\n";
                    $line .= "${SPACE}" ."  end if\n";
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
    if ( $DEBUG ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING: $dim, $type, $subt"; }

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
    if ( $DEBUG ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_DIAGG: $dim, $type, $subt"; }

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
    if ( $DEBUG ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_FIELD: $subt, $spec"; }

    my $line;
    
    my $SPEC       = uc($spec);
    my $spec_short = substr($spec,0,3);

    if( ! defined $STRING_INDEX_ARRAY{"diagnos_${subt}_${spec}_2_S"} ){
        $STRING_INDEX_ARRAY{"diagnos_${subt}_${spec}_2_S"} = $STRING_INDEX;
        #print Dumper(\%STRING_INDEX_ARRAY) . "\n";
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
    if ( $DEBUG ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_INDEX:"; }

    my $line;

    my $search_S = "";
    my $search_E = "";
    my $out      = "";

    my $subt_short = ucfirst(substr($subt,0,3));

    my $type_short = ucfirst(substr($type,0,4));
    if( $type_short eq "Stat" ){ 
        $type_short = 'State';  #fix because "state" has 5 chars
    }else{
        if( $dim == 2 ){ $type_short .= '2d'; } #fix because name 2D not collapse with 3D
    }
    $search_S = "${type}_${subt}_${dim}_S";
    $search_E = "${type}_${subt}_${dim}_E";
    $out = "${type_short}";

    if( $type eq "diagnos" ){
        #fix because "diagnos" ends with "diaggrp" values
        if ( defined $STRING_INDEX_ARRAY{"diaggrp_${subt}_${dim}_E"} ){
            $search_E = "diaggrp_${subt}_${dim}_E";
        }
    }
    # if( $subt eq 'pel' && $type eq 'state' && $dim == 2 ){
    #     #fix because "Pel state 2D" ends with "Pel diagnos river 2D" values
    #     $search_E = "diagnos_${subt}_river_2_E";
    # }

    if( defined $STRING_INDEX_ARRAY{"${search_S}"} && defined $STRING_INDEX_ARRAY{"${search_E}"} ){
        $line .= "${SPACE}st${subt_short}${out}S=" . $STRING_INDEX_ARRAY{"${search_S}"} . "\n"; 
        $line .= "${SPACE}st${subt_short}${out}E=" . $STRING_INDEX_ARRAY{"${search_E}"} . "\n";
    }else{
        print "WARNING: \"${search_S}\" and \"${search_E}\" do not exist. No print value in set_var_info\n";
        #print Dumper(\%STRING_INDEX_ARRAY) . "\n";
    }

    if( $line ){ print $file $line; }
}

sub func_STRING_INDEX_FIELD  {
    my ( $file, $subt, $spec ) = @_;
    if ( $DEBUG ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_INDEX_FIELD:"; }

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

    if( defined $STRING_INDEX_ARRAY{"${search_S}"} && defined $STRING_INDEX_ARRAY{"${search_E}"} ){
        $line .= "${SPACE}st${subt_short}${out}S=" . $STRING_INDEX_ARRAY{"${search_S}"} . "\n"; 
        $line .= "${SPACE}st${subt_short}${out}E=" . $STRING_INDEX_ARRAY{"${search_E}"} . "\n";
    }

    if( $line ){ print $file $line; }
}

sub func_STRING_INDEX_STARTEND  {
    my ( $file, $subt ) = @_;
    if ( $DEBUG ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_INDEX_STARTEND:"; }

    my $line = '';
    my $subt_short = '';
    if( ! $subt ){ $subt = ''; }
    else{ $subt_short = ucfirst(substr($subt,0,3)); }

    my $max=0;
    my $min=0;
    foreach my $key ( keys %STRING_INDEX_ARRAY ){
        if( $key =~ m/.*${subt}.*_E/ ){
            if( $max < $STRING_INDEX_ARRAY{$key} ){ 
                $max = $STRING_INDEX_ARRAY{$key};
            }
        }
        if( $key =~ m/.*${subt}.*_S/ ){
            if( !$min || $min > $STRING_INDEX_ARRAY{$key} ){
                $min = $STRING_INDEX_ARRAY{$key};
            }
        }
    }

    if( $min && $max ){
            $line .= "${SPACE}st${subt_short}Start=". $min . "\n";
            $line .= "${SPACE}st${subt_short}End=". $max . "\n";
    }

    if( $line ){ print $file $line; }
    #print Dumper(\%STRING_INDEX_ARRAY) . "\n";
}



###########################################  INIT_VAR_BFM ######################


sub func_INIT_VALUE_CALC {
    my ( $file, $before, $subt, $value ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_VALUE_CALC: "; }

    my $line = '';

    if ( !$value ){ $value = 0; }

    foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$groupname};
        if( $subt eq $group->getSubtype()){

            $line .= "${before}" . join(',', 'Calc' . $group->getSigla()) . " = " . ${value} . "\n";
        }
    }

    if( $line ){ print $file $line; }
}


sub func_INIT_PP  {
    my ( $file, $before, $dim, $type, $subt ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_PP: "; }

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type, $subt)  ){ return; }

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

sub func_INIT_RATIO_CONST {
    my ( $file ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_RATIO_CONST: "; }

    my $line = '';
    my $defratio = '';
    my @constList = ();
    my @constNoCList = ();
    my @constOptionalList = ();
    my @constOptionalRatioList = ();

    foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
        push(@constList, $const);
        if($const ne $constList[0]){
            push( @constNoCList, $const );
            push( @constOptionalList, $const . $constList[0] );
        }
    }
    foreach my $constOpt (@constOptionalList){
     $defratio = "ZERO";
     if($constOpt eq "nc") {$defratio = "0.0126_RLEN    ! Redfield" };
     if($constOpt eq "pc") {$defratio = "0.7862e-3_RLEN ! Redfield" };
     if($constOpt eq "sc") {$defratio = "0.0145_RLEN    ! Redfield" };
     if($constOpt eq "lc") {$defratio = "0.03_RLEN      ! standard diatom value" };
     if($constOpt eq "fc") {$defratio = "3.e-04_RLEN    ! standard diatom value" };
     $line .="   real(RLEN),parameter :: " . $constOpt . "_ratio_default = " . $defratio . "\n";
    }

    if( $line ){ print $file $line; }
}


sub func_INIT_FUNC_CONST {
    my ( $file ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_FUNC_CONST: "; }

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


sub func_INIT_NAMELIST {
    my ( $file, $subt, $num1, $num2 ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_NAMELIST: "; }

    my $line = '';

    my $subtype= '_' . $subt;
    if ( $subtype eq '_pel' ){ $subtype = '' } #fix because pel is default and vars has no suffix

    $line .= "${SPACE}icontrol=0\n";
    $line .= "${SPACE}open(namlst,file=bfm_init_fname${subtype},action='read',status='old',err=${num1})\n";
    $line .= "${SPACE}read(namlst,nml=bfm_init_nml${subtype},err=${num2})\n";
    $line .= "${SPACE}close(namlst)\n";
    $line .= "${SPACE}icontrol=1\n";
    $line .= "${num1} if ( icontrol == 0 ) then\n";
    $line .= "${SPACE}  LEVEL3 'I could not open ',trim(bfm_init_fname${subtype})\n";
    $line .= "${SPACE}  LEVEL3 'The initial values of the BFM variables are set to ZERO'\n";
    $line .= "${SPACE}  LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname${subtype})\n";
    $line .= "${SPACE}  icontrol=1\n";
    $line .= "${SPACE}end if\n";
    $line .= "${num2} if ( icontrol == 0 ) then\n";
    $line .= "${SPACE}  FATAL \'I could not read bfm_init_nml${subtype}\'\n";
    $line .= "${SPACE}  stop \'init_var_bfm\'\n";
    $line .= "${SPACE}end if\n";
   
    if( $line ){ print $file $line; }
}


sub func_INIT_OUTPUT_VARIABLES {
    my ( $file, $dim, $subt ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_OUTPUT_VARIABLES: "; }

    my $line = '';

    my $Subtype= ucfirst($subt);
    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    $line .= "${SPACE}   write(Flun,155) \'ID\',\'Var\',\'Unit\',\'Long Name\',\'Flag\'\n";
    $line .= "${SPACE}   do n=st${Subtype}StateS,st${Subtype}StateE\n";
    $line .= "${SPACE}     write(Flun,156) n,trim(var_names(n)),trim(var_units(n)) &\n";
    $line .= "${SPACE}       ,trim(var_long(n)),D${dim}STATETYPE${SUBTYPE}(n-st${Subtype}StateS+1)\n";
    $line .= "${SPACE}   end do\n";
   
    if( $line ){ print $file $line; }
}


sub func_INIT_FUNC_ZERO {
    my ( $file, $dim, $subt ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_FUNC_ZERO: "; }
    my $line = '';

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, 'state', $subt)  ){ return; }

    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$groupname};

        if( $group->getDim() == $dim
            && $subt eq $group->getSubtype()){

            $line .= "${SPACE}do j = 1, ii" . $group->getSigla()    . "\n";
            $line .= "${SPACE}  if (.NOT.Calc" . $group->getSigla() ."(j)) then\n";
            $line .= "${SPACE}    do i = 1,iiLastElement\n";
            $line .= "${SPACE}      if ( pp" . $group->getSigla() . "(j,i) /= 0 ) then \n";
            $line .= "${SPACE}        D${dim}STATE${SUBTYPE}(pp" . $group->getSigla() . "(j,i),:) = p_small\n";
            $line .= "${SPACE}        D${dim}STATETYPE${SUBTYPE}(pp" . $group->getSigla() . "(j,i)) = OFF\n";
            $line .= "#if defined key_obcbfm\n";
            $line .= "${SPACE}        D${dim}STATEOBC${SUBTYPE}(pp" . $group->getSigla() . "(j,i)) = NOOBCSTATES\n";
            $line .= "#endif\n";
            $line .= "${SPACE}      end if\n";
            $line .= "${SPACE}    end do\n";
            $line .= "${SPACE}  end if\n";
            $line .= "${SPACE}end do\n\n";
        }
    }

    if( $line ){ print $file $line; }
}

sub func_INIT_DEFAULT  {
    my ( $file, $dim, $type, $subt ) = @_;
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_DEFAULT: "; }

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
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INITS_SETS: "; }

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
    if ( $DEBUG ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_INTERNAL: "; }

    my $line = '';
    my @constList         = ();
    my @constNoC          = ();
    my @constOptionalList = ();

    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
        push(@constList, $const);
        if($const ne $constList[0]){ 
            push( @constNoC         , $const                 );
            push( @constOptionalList, $const . $constList[0] ); 
        }
    }


    foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group         = $$LST_GROUP{$groupname};
        my $groupAcro     = $$LST_GROUP{$groupname}->getAcro();
        my $groupname_nml = $groupname; $groupname_nml =~ s/Plankton//; $groupname_nml .= "_parameters";

        if( $group->getDim() == $dim
            && $subt eq $group->getSubtype()
            && exists ${$$LST_GROUP{$groupname}->getComponents()}{'c'}
            && ! ( keys %{$$LST_GROUP{$groupname}->getComponents()} == 2 && exists ${$$LST_GROUP{$groupname}->getComponents()}{'h'} ) ){

            $line .= "${SPACE}do i = 1 , ( ii". $groupname ." )\n";
            $line .= "${SPACE}  call init_constituents( c=" . $groupname  . "(i,iiC), &\n${SPACE}";
            my @temp_line;
            foreach my $const (@constNoC){
                if( exists ${$$LST_GROUP{$groupname}->getComponents()}{$const} ){
                    push( @temp_line, "    " . $const . "=D" . $dim . "STATE" . $SUBTYPE . "(pp" . $groupname . "(i,ii" . uc($const) . "),:)" );
                }
            }
            my @temp_line2;
            foreach my $constOpt (@constOptionalList){
                #get the constituents active inside the group
                if( exists ${$$LST_GROUP{$groupname}->getComponents()}{ (split('',$constOpt))[0] } ){
                    my $temp_compo = "p_q" . $constOpt . $groupAcro;
                    #print $temp_compo . " " . $groupname_nml . " ";
                    #if the optional initialization element exists in the namelist => add to initialize constituents
                    foreach my $list (@$LST_NML){
                        # search for namelists starting with groupname_paramters*
                        if( $list->name() =~ m/^$groupname_nml/ ){
                            foreach my $param ( @{$list->slots()} ){
                                if( $param eq $temp_compo){ 
                                    push( @temp_line2, $constOpt . "=" . $temp_compo . "(i)" );
                                }
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
    if ( $DEBUG ){ print "INCLUDE -> FUNCTION CALLED func_HEADER: "; }

    my $line = '';
    my $TYPE = uc($type);

    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        $line .= "#define ${nameC}(A) D${dim}${TYPE}${SUBTYPE}(pp${nameC},A)\n";
                    }
                }
            }else{
                $line .= "#define ${name}(A) D${dim}${TYPE}${SUBTYPE}(pp${name},A)\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_GROUP_HEADER  {
    my ( $file, $dim, $subt ) = @_;
    if ( $DEBUG ){ print "INCLUDE -> FUNCTION CALLED func_GROUP_HEADER: "; }

    my $line = '';

    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    foreach my $name ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$name};
        if( $dim == $group->getDim() 
            && $subt eq $group->getSubtype()){
            $line .= "#define ${name}(A,B) D${dim}STATE${SUBTYPE}(pp${name}(A,B),:)\n";
        }
    }

    if( $line ){ print $file $line; }
}


sub func_HEADER_FIELD  {
    my ( $file, $subt, $spec) = @_;
    if ( $DEBUG ){ print "INCLUDE -> FUNCTION CALLED func_HEADER_FIELD: "; }

    if( $subt ne 'pel' ){ print "ERROR: Subtype can only be \'pel\' for field $spec\n"; exit; }

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
    if ( $DEBUG ){ print "INCLUDE -> FUNCTION CALLED func_HEADER_DIAGG: "; }

    my $line = '';

    my $TYPE = uc($type);
    my $SUBTYPE= '_' . uc($subt);
    if ( $SUBTYPE eq '_PEL' ){ $SUBTYPE = '' } #fix because pel is default and vars has no suffix

    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() 
            && $type eq $param->getType()
            && $subt eq $param->getSubtype() ){
            $line .= "#define ${name}(A,B) D${dim}DIAGNOS${SUBTYPE}(pp${name}(A),B)\n";
        }
    }

    if( $line ){ print $file $line; }
}
