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
#use List::MoreUtils qw(firstidx);
#use List::Util qw(minstr maxstr);
use Data::Dumper;

use classes;



########### VARIABLES TO EXPORT ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(print_f90);
########### VARIABLES TO EXPORT ##########################


########### VARIABLES GLOBAL ##########################
my $VERBOSE = 0;
my $SPACE = "    ";
my $MAX_CHAR_PER_LINE = 76;
my $STRING_INDEX=1;
my %STRING_INDEX_ARRAY=();
my ( $LST_PARAM, $LST_GROUP, $LST_STA, $LST_CONST );
########### VARIABLES GLOBAL ##########################





########### REGULAR EXPRESSIONS ##########################
my $dispatch = {
    #COMMENT LINE FOR PROGRAMMERS
    '\%\!.*' => \&func_COMMENT,
    #DESCRIPTION
    '\%(3|2)d-(diagnos|state)-desc(?:\s|\n)'=> \&func_DESC ,
    '\%(3|2)d-(diaggrp)-desc(?:\s|\n)' => \&func_DESC_DIAGG ,
    #NR
    '([^\%]*)\%(3|2)d-(flux|diagnos|state)-nr(?:\s|\n)'=>    \&func_NR ,
    #ARRAY
    '\%(3|2)d-(diagnos|state)-array(?:\s|\n)'=>    \&func_ARRAY ,
    '\%(3)d-(state)-field-array\s*(\S*)(?:\s|\n)'=>    \&func_ARRAY_FIELD ,
    #PP
    '\%(3|2)d-(diagnos|state)-pp(?:\s|\n)'=>    \&func_PP ,
    '\%(3|2)d-(diaggrp)-pp(?:\s|\n)'=>    \&func_PP_DIAGG ,
    '\%(3|2)d-(diaggrp)-assign-pp(?:\s|\n)'=>    \&func_PP_ASSIGN ,
    #POINTER
    '\%(3|2)d-(diagnos|state)-pointer(?:\s|\n)'=>    \&func_POINT ,
    '\%(3|2)d-(diaggrp)-pointer(?:\s|\n)'=>    \&func_POINT_DIAGG ,
    '\%(3)d-Z-pointer(?:\s|\n)'=>    \&func_POINT_Z ,
    '\%(3)d-(state)-field-pointer\s*(\S*)(?:\s|\n)'=>    \&func_POINT_FIELD ,
    '\%(3|2)d-(diagnos|state)-alloc-pointer(?:\s|\n)' =>    \&func_POINT_ALLOC ,
    '\%(3|2)d-(diaggrp)-alloc-pointer(?:\s|\n)' =>    \&func_POINT_ALLOC_DIAGG ,
    '\%(3)d-(state)-field-alloc-pointer\s*(\S*)(?:\s|\n)' =>    \&func_POINT_ALLOC_FIELD ,
    #HEADER
    '\%(3|2)d-(diagnos|state)-header(?:\s|\n)' =>    \&func_HEADER ,
    '\%(3|2)d-(diaggrp)-header(?:\s|\n)' =>    \&func_HEADER_DIAGG ,
    '\%(3)d-(state)-field-header\s*(\S*)(?:\s|\n)' =>    \&func_HEADER_FIELD ,
    '\%if-exist-header(?:\s|\n)' =>    \&func_HEADER_IF ,
    #STRING
    '\%(3|2)d-(diagnos|state|flux)-string(?:\s|\n)'=>    \&func_STRING ,
    '\%(3|2)d-(diaggrp)-string(?:\s|\n)'=>    \&func_STRING_DIAGG ,
    '\%(3)d-(state)-field-string\s*(\S*)(?:\s|\n)'=>    \&func_STRING_FIELD ,
    '\%dd-string-index(?:\s|\n)'=>    \&func_STRING_INDEX ,
    #ALLOC
    '\%(3|2)d-(diagnos|state)-alloc(?:\s|\n)' => \&func_ALLOC ,
    '\%(3)d-Z-alloc(?:\s|\n)' => \&func_ALLOC_Z ,
    '\%(1|3|2)d-(intvar|variable|variables)-alloc(?:\s|\n)'=>    \&func_ALLOC_INTVAR ,
    #FLUX
    '\%(dd)-(flux)-alloc(?:\s|\n)' => \&func_FLUX_ALLOC ,
    '\%(3|2)d-(flux)-fill(?:\s|\n)' => \&func_FLUX_FILL ,
    #CONSTITUENT
    '\%constituent(?:\s|\n)' => \&func_CONSTITUENT ,
    #GROUP    
    '\%(3|2)d-group-header(?:\s|\n)' => \&func_GROUP_HEADER ,
    '\%(3|2)d-group-parameter(?:\s|\n)' => \&func_GROUP_PARAMETER ,
    '\%(3|2)d-group-function-name(?:\s|\n)' => \&func_GROUP_FUNCTION_NAME ,
    '\%(3|2)d-groupfunctions(?:\s|\n)' => \&func_GROUP_FUNCTIONS ,
    #INTVAR/VARIABLE
    '\%(3|2|1)d-(intvar|variable)(?:\s|\n)' => \&func_INTVAR ,
    #INIT
    '\%(3)d-(state)-Initpp(?:\s|\n)'      => \&func_INIT_PP ,
    '\%(3)d-(state)-InitDefault(?:\s|\n)' => \&func_INIT_DEFAULT ,
    '\%(3)d-(state)-InitSets(?:\s|\n)'    => \&func_INIT_SETS ,
};
########### REGULAR EXPRESSIONS ##########################


########### FUNCTIONS ##########################

sub print_f90 {

    my ($templ_mem, $output, $lst_group, $lst_param, $lst_sta, $lst_const, $verbose) = @_;

    $LST_PARAM = $lst_param;
    $LST_GROUP = $lst_group;
    $LST_STA   = $lst_sta;
    $LST_CONST = $lst_const;
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
                my ($dim, $type, $spec);
                $dim = $1; $type=$2; $spec=$3;
                if($spec){ &{$$dispatch{$key}}($OUT, $dim, $type, $spec); }
                else     { &{$$dispatch{$key}}($OUT, $dim, $type);        }
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
    my ($dim, $type) = @_;

    my $size = 0;

    foreach my $name (keys %$LST_PARAM){
        my $param    = $$LST_PARAM{$name};
        my $dim_var  = $param->getDim();
        my $type_var = $param->getType();
        if( $type_var eq 'diaggrp' ){ $type_var = 'diagnos'; }

        if( $dim_var && $dim_var == $dim && $type_var eq $type ){ 
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
    my ($dim, $type) = @_;

    foreach my $key (keys %$LST_PARAM){
        my $param    = $$LST_PARAM{$key};
        my $dim_var  = $param->getDim();
        my $type_var = $param->getType();

        if( $dim_var && $dim_var == $dim && $type_var eq $type ){ 
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_ALLOC: "; }

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type)  ){ return; }
    
    my $j = "";
    my $TYPE = uc($type);
    if ( $dim == 2 ){ $j = "_XY"; }

    print $file "${SPACE} \n";
    print $file "${SPACE}allocate(D${dim}${TYPE}(1:NO_D${dim}_BOX_${TYPE}S,1:NO_BOXES$j),stat=status)\n";
    print $file "${SPACE}if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}$TYPE\")\n";
    print $file "${SPACE}D${dim}${TYPE} = ZERO\n";
    if ( $type eq "state" ) {
        print $file "#ifdef D1SOURCE\n";
        print $file "${SPACE}  allocate(D${dim}SOURCE(1:NO_D${dim}_BOX_STATES,1:NO_BOXES$j),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SOURCE\")\n";
        print $file "${SPACE}  D${dim}SOURCE = ZERO\n";
        print $file "#else\n";
        print $file "${SPACE}  allocate(D${dim}SOURCE(1:NO_D${dim}_BOX_STATES,1:NO_D${dim}_BOX_STATES,1:NO_BOXES$j),stat=status)\n";
        print $file "${SPACE}  if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SOURCE\")\n";
        print $file "${SPACE}  D${dim}SOURCE = ZERO\n";
        print $file "#ifndef ONESOURCE\n";
        print $file "${SPACE}    allocate(D${dim}SINK(1:NO_D${dim}_BOX_STATES,1:NO_D${dim}_BOX_STATES,1:NO_BOXES$j) ,stat=status)\n";
        print $file "${SPACE}    if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"D${dim}SINK\")\n";
        print $file "${SPACE}    D${dim}SINK = ZERO\n";
        print $file "#endif\n";
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_POINT_ALLOC: "; }

    my $line = "";
    my $TYPE = uc($type);

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        #foreach my $compo (sort keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_PP_ASSIGN: "; }
    
    my $line = "";
    my $TYPE = uc($type);

    #calculate lenght
    my $len=0;
    # foreach my $name (sort keys %$LST_PARAM){
    #     my $param = $$LST_PARAM{$name};
    #     if( $dim == $param->getDim() && $param->getType() eq 'diagnos' ){
    #         $len += sizeWithCompo($param);
    #     } 
    # }

    # if( exists $$LST_STA{"diagnos ${dim}d"} ){
    #     $len = $$LST_STA{"diagnos ${dim}d"};
    # }

    #foreach my $root (sort keys %$LST_PARAM){
    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() && $param->getQuota() ){
            #print $root . " -> " . $param->getQuota() . "\n";
            #foreach my $member (sort keys %$LST_PARAM){
            foreach my $member ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                my $param2 = $$LST_PARAM{$member};
                if( defined($param2->getGroup()) && $param2->getGroup() eq $param->getQuota() ){
                    #$line .= "pp${root}(ii${member})=" . ++$len . "\n";
                    $line .= "pp${root}(ii${member})=" . ++($$LST_STA{"diagnos ${dim}d"}) . "\n";
                }
            }
        }
    }

    if( $line ){ print $file $line; }
}

sub func_POINT_ALLOC_DIAGG  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_PP_ALLOC_DIAG: "; }
    
    my $line = "";

    #foreach my $root (sort keys %$LST_PARAM){
    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() && $param->getQuota() ){
            #print $root . " -> " . $param->getQuota() . "\n";
            my @namesMem = ();
            #foreach my $member (sort keys %$LST_PARAM){
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
    my ( $file, $dim, $type, $spec) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_POINT_ALLOC_FIELD: "; }

    my $line_ini = "";
    my $line_par = "";
    
    my $SPEC       = uc($spec);
    my $spec_short = substr($spec,0,3);

    my $n = 0;
    my $m = sizeDimType($dim, $type);
    if( exists $$LST_STA{"${type} ${dim}d index"} ){
        $n = $$LST_STA{"${type} ${dim}d index"};
        $$LST_STA{"${type} ${dim}d index"} += $m;
    }else{
        $n = $$LST_STA{'diagnos 2d'};
        $$LST_STA{"${type} ${dim}d index"} = $n + $m;
    }
    $$LST_STA{"${type} ${dim}d ${spec} index"} = $n;

    $line_ini .= "${SPACE}PEL${SPEC} => D2DIAGNOS(${n}+1:${n}+${m},:); PEL${SPEC}=ZERO\n";
    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        #foreach my $compo (sort keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_ALLOC_INTVAR: "; }
    
    my $line = "";
    my $j = "";
    my $TYPE = uc($type);
    if ( $dim == 2 ){ $j = "_XY"; }

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            my $mode = 0;
            my $l = "";
            my $groupindex = "";
            if( $param->getDim() > 1 ){ $mode += 2; }
            if( $param->getQuota()   ){ 
                $mode += 1; 
                #my @group = (sort keys %$LST_GROUP);
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_FLUX_ALLOC: "; }

    my $line = "";    
    my ($nn, $mm);
    $nn=0, $mm=0;

    foreach my $dim_tmp ( qw(2d 3d) ){
        if( exists $$LST_STA{"flux $dim_tmp"} && exists $$LST_STA{"select $dim_tmp"} ){
            $nn += $$LST_STA{"flux $dim_tmp"};
            $mm += $$LST_STA{"select $dim_tmp"};
        }
    }

    $line .= "allocate(flx_calc_nr(0:$nn),stat=status)\n";
    $line .= "allocate(flx_CalcIn(1:$nn),stat=status)\n";
    $line .= "allocate(flx_option(1:$nn),stat=status)\n";
    $line .= "allocate(flx_t(1:$mm),stat=status)\n";
    $line .= "allocate(flx_SS(1:$mm),stat=status)\n";
    $line .= "allocate(flx_states(1:$mm),stat=status)\n";
    $line .= "allocate(flx_ostates(1:$mm),stat=status)\n";
    $line .= "flx_calc_nr(0)=0\n";
    $line .= "flx_cal_ben_start=$nn\n";

    if( $line ){ print $file $line; }
}


sub func_FLUX_FILL  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_FLUX_ALLOC: "; }

    my $line = "";
    my $index=1;
    my $jndex=1;

    if( exists $$LST_STA{"flux ${dim}d"} && exists $$LST_STA{"select ${dim}d"} ){
        #foreach my $key (sort keys %$LST_PARAM) { 
        foreach my $key ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
            my $param    = $$LST_PARAM{$key};
            my $function = $param->getFunction();
            if( $param->getType() eq 'flux' ){ 
                $line .= "\n";
                $line .= "! " . $param->getSigla() . "=$$function{xpr}        (normal flux): \n";
                $line .= "flx_calc_nr($index)= " . ($#{$$function{compo1}} + $jndex) . "; ";
                $line .= "flx_CalcIn($index)=iiPel; flx_option($index)=0\n";
                for my $indexC ( 0 .. $#{$$function{compo1}} ) {
                    my $sign1  = $$function{sign2}[$indexC];
                    my $dir    = $$function{dir}[$indexC];
                    my $compo1 = $$function{compo1}[$indexC];
                    my $compo2 = $$function{compo2}[$indexC];
                    if( $compo2 eq '*' ){ $compo2 = $compo1; }
                    $line .= "flx_t($jndex)=${sign1}1.00;flx_SS($jndex)=${dir}; ";
                    $line .= "flx_states($jndex)=pp${compo1};flx_ostates($jndex)=pp${compo2}\n";
                    $jndex++;
                }
                $index++;
            }
        }
    }


    if( $line ){ print $file $line; }
}



sub func_ALLOC_Z  {
    my ( $file, $dim) = @_;
    if ( $VERBOSE ){ print "AllocateMem -> FUNCTION CALLED func_ALLOC_Z: "; }

    my %z_hash = ();

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $param->getZ() ){
            my @z_array_tmp = ( $param->getZ() =~ /\s*(.*)=.*/ );
            #get unique values from array
            $z_hash{"$z_array_tmp[0]"} = 0;
        }
    }

    if( keys %z_hash ){
        print $file "${SPACE}!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
        print $file "${SPACE}!Allocations(s) of and assigning values to alternative  Z-axis\n";
        print $file "${SPACE}\n";
        foreach my $key (keys %z_hash){
            print $file "${SPACE}\n";
            print $file "${SPACE}allocate(${key}(1:NO_BOXES_Z),stat=status)\n";
            print $file "${SPACE}if (status /= 0) call error_msg_prn(ALLOC,\"AllocateMem\", \"${key}\")\n";
            print $file "${SPACE}!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
        }
    }

}




###########################################  MODULE_MEM FUNCTIONS ######################


sub func_DESC  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_DESC: "; }

    my $line     = "";
    my $line_par = "";


    $line = sprintf "! %10s %60s %15s", "${dim}d name", "description", "unit\n";
    $line .= "! ";
    foreach (1..10){ $line .=  "-"; } $line .=  " ";
    foreach (1..60){ $line .=  "-"; } $line .=  " ";
    foreach (1..15){ $line .=  "-"; } $line .=  "\n";
    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        my $comment = $param->getComment();
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        #foreach my $compo (sort keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
                        #my $unit = ${$param->getComponents()}{$compo};
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


sub func_ARRAY  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_ARRAY: "; }

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type)  ){ return; }

    my $line = "${SPACE}real(RLEN),public,pointer,dimension(:,:) :: D${dim}" . uc($type) . "\n";
    if ( $type eq "state" ) {
        $line .= "#ifdef D1SOURCE\n";
        $line .= "${SPACE}real(RLEN),public,pointer,dimension(:,:) :: D${dim}" . "SOURCE\n";
        $line .= "${SPACE}real(RLEN),public,pointer,dimension(:,:) :: D${dim}" . "SINK\n";
        $line .= "#else\n";
        $line .= "${SPACE}real(RLEN),public,pointer,dimension(:,:,:) :: D${dim}" . "SOURCE\n";
        $line .= "${SPACE}real(RLEN),public,pointer,dimension(:,:,:) :: D${dim}" . "SINK\n";
        $line .= "#endif\n";
        $line .= "${SPACE}integer,public,pointer,dimension(:) :: D${dim}" . "STATETYPE\n";
        $line .= "#ifdef BFM_NEMO\n";
        $line .= "${SPACE}integer,public,pointer,dimension(:) :: D${dim}" . "STATEOBC\n";
        $line .= "#endif\n";
    }

    print $file $line;
}

sub func_ARRAY_FIELD  {
    my ( $file, $dim, $type, $spec) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_ARRAY_FIELD: "; }

    #if the variable does not exists => dont print anything
    if( ! checkDimType($dim, $type)  ){ return; }

    #save the spec in stadistics for further calculations of number of variables
    if( ! exists $$LST_STA{'spec num'} ){ $$LST_STA{'spec num'} = 1; }
    else{ $$LST_STA{'spec num'} += 1; }

    my $SPEC = uc($spec);
    
    print $file "${SPACE}real(RLEN),public,pointer,dimension(:,:) :: PEL${SPEC}\n";
}

sub func_NR  {
    my ( $file, $before, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_NR: "; }

    my $title  = "${type} ${dim}d";
    my $number = 0;
    if( exists $$LST_STA{$title} ){
        if( $title eq 'diagnos 2d' ){
            $number = $$LST_STA{"${title}"} + sizeGroup($dim) + ( $$LST_STA{'state 3d'} * $$LST_STA{'spec num'} );
            #if( exists $$LST_STA{'state 3d index'} ){ $number = $$LST_STA{'state 3d index'}; }
            # if( exists $$LST_STA{'state 2d'} ){ $number += $$LST_STA{'state 2d'}; }
            # if( exists $$LST_STA{'state 3d'} ){ $number += $$LST_STA{'state 3d'}; }
            # $number = $number*3 + $$LST_STA{$title};
        }
        elsif( $title eq 'diagnos 3d' ){ 
            $number=sizeDimType($dim, $type); 
        }
        else{ 
            $number=$$LST_STA{$title}; 
        }
    }else{ $number = 0; }
    #print " $title: $number ";

    print $file "${SPACE}$before${number}\n";
}


sub func_PP  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_PP: "; }

    my $line     = "";
    my @line_par = ();


    my $index = 1;
    $line .= "${SPACE}integer,parameter,public :: ";
    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
                        push( @line_par, "pp${nameC}=".$index++ );
                    }
                }
                if( $param->getComponentsEx() && keys(%{$param->getComponentsEx()}) != 0 ){
                    foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                        if( exists ${$param->getComponentsEx()}{$const} ){
                            my $nameC = $name . $const;
                            #foreach my $compo (keys %{$param->getComponentsEx()} ){
                            #my $nameC = $name . $compo;
                            push( @line_par, "pp${nameC}=0" );
                        }
                    }
                }
            }else{
                push( @line_par, "pp${name}=".$index++);
            }
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $line); print $file "\n"; }

}


sub func_POINT  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_POINT: "; }
    
    my $line     = "";
    my @line_par = ();

    $line .= "${SPACE}real(RLEN),public,dimension(:),pointer  :: ";
    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_CONSTITUENT: "; }
    
    my $line     = "";
    my @line_par = ();

    $line .= "${SPACE}integer,parameter,public :: ";
    foreach my $key ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
        push( @line_par, "ii" . uc($key) . "=" . $$LST_CONST{$key} );
    }

    if( $#line_par > 0 ){ printList($file, \@line_par, $line); print $file "\n"; }
}


sub func_GROUP_PARAMETER  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_GROUP_PARAMETER: "; }

    #foreach my $group_name (sort keys %$LST_GROUP){
    foreach my $group_name ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ) {
        my $group = $$LST_GROUP{$group_name};
        if( $dim == $group->getDim() ){
            my @elements = ();
            my $index = 1;
            #foreach my $param_name (sort keys %$LST_PARAM){
            foreach my $param_name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
                my $param = $$LST_PARAM{$param_name};
                if( $param->getGroup() && $group_name eq $param->getGroup() ){ push(@elements, ("ii$param_name=". $index++ ) ); }
            }
            if( $#elements == -1 ){
                print $file "${SPACE}integer,parameter,public     :: ii$group_name=". scalar(@elements) . "\n" ;
            }else{
                print $file "${SPACE}integer,parameter,public     :: ii$group_name=". scalar(@elements) . ", " . join(", ",@elements) . "\n" ;
            }
        }
    }
}


sub func_INTVAR  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_INTVAR: "; }

    my $line_ini = "";
    my $line_par = "";

    if    ( $type eq 'variable') { $line_ini = "${SPACE}real(RLEN),public"; } # :: &\n
    elsif ( $type eq 'intvar'  ) { $line_ini = "${SPACE}integer,public"   ; } # :: &\n
    else                         { print "ERROR: $type not found\n"; exit 1;  }

    #foreach my $param_name (sort keys %$LST_PARAM){
    foreach my $param_name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$param_name};
        if( $param->getDim() == $dim && ($param->getType() eq $type) ){
            my $mode  = 0;
            if( $param->getQuota() ){    $mode++; }
            if( $param->getDim () > 1 ){ $mode++; }
            
            if    ($mode == 0) { $line_par .= "${line_ini}                             :: $param_name ! "; }
            elsif ($mode == 1) { $line_par .= "${line_ini},dimension(:),allocatable    :: $param_name ! "; }
            elsif ($mode == 2) { $line_par .= "${line_ini},dimension(:,:),allocatable  :: $param_name ! "; }

            if( $param->getComment() ){ $line_par .= $param->getComment();           }
            if( $param->getUnit() ){    $line_par .= " (" . $param->getUnit()  . ")";}
            $line_par .= "\n";
        }
    }

    if( $line_par ){
        #$line_par =~ s/(.*)\,\s\&(.*)$/$1$2/; #substitute the last &
        print $file $line_par;
    }
    
}

sub func_DESC_DIAGG  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_DESC_DIAGG: "; }

    #foreach my $root (sort keys %$LST_PARAM){
    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() && $param->getQuota() ){
            #print $root . " -> " . $param->getQuota() . "\n";
            #foreach my $member (sort keys %$LST_PARAM){
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_PP_DIAGG: "; }

    my $line     = "";
    my @line_par = ();

    $line .= "${SPACE}integer,public :: ";
    #foreach my $root (sort keys %$LST_PARAM){
    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() && $param->getQuota() ){
            push(@line_par,"pp${root}(ii". $param->getQuota() .")" );
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $line); print $file "\n"; }
}


sub func_POINT_DIAGG  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_POINT_DIAGG: "; }
    
    my $line     = "";
    my @line_par = ();

    $line .= "${SPACE}real(RLEN),public,dimension(:,:),pointer  :: ";
    #foreach my $root (sort keys %$LST_PARAM){
    foreach my $root ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$root};
        if( $dim == $param->getDim() && $param->getQuota() ){
            push( @line_par, "${root}" );
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $line); print $file "\n"; }

}

sub func_POINT_FIELD  {
    my ( $file, $dim, $type, $spec) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_POINT_FIELD: "; }

    my $line     = "";
    my @line_par = ();

    my $SPEC = substr($spec,0,3);
    $line .= "${SPACE}real(RLEN),public,dimension(:),pointer  :: ";
    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $param->getDim() == $dim && ($param->getType() eq $type) ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #    my $nameC = $name . $compo;
                        push( @line_par, "j${SPEC}${nameC}");
                    }
                }
            }else{
                push( @line_par, "j${SPEC}${name}" );
            }
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $line); print $file "\n"; }
}


sub func_GROUP_FUNCTION_NAME  {
    my ( $file, $dim) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_GROUP_FUNCTION_NAME: "; }

    my $line     = "";
    my @line_par = ();

    $line .= "${SPACE}public ";
    #foreach my $group_name (sort keys %$LST_GROUP){
    foreach my $group_name ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group = $$LST_GROUP{$group_name};
        if( $group->getDim() == $dim ){
            push( @line_par, "pp${group_name}" );
            push( @line_par, "${group_name}" );
        }
    }

    if( $#line_par >= 0 ){ printList($file, \@line_par, $line); }
}


sub func_GROUP_FUNCTIONS  {
    my ( $file, $dim) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_GROUPFUNCTIONS: "; }

    foreach my $pre ("pp", ""){
        #foreach my $groupname (sort keys %$LST_GROUP){
        foreach my $groupname ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
            my $group   = $$LST_GROUP{$groupname};

            if( $group->getDim() == $dim ){
                my @members = ();
                #foreach my $paramname (sort keys %$LST_PARAM){ 
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

sub func_POINT_Z  {
    my ( $file, $dim) = @_;
    if ( $VERBOSE ){ print "ModMem -> FUNCTION CALLED func_POINT_Z: "; }

    my %z_hash = ();

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $param->getZ() ){
            my @z_array_tmp = ( $param->getZ() =~ /\s*(.*)=.*/ );
            #get unique values from array
            $z_hash{"$z_array_tmp[0]"} = 0;
        }
    }

    if( keys %z_hash ){ 
        print $file "${SPACE}!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
        print $file "${SPACE}! Definition(s) of alternative Z-axis\n";
        print $file "${SPACE}real(RLEN),public,dimension(:),pointer  :: " . join(", ", keys %z_hash),"\n";
        print $file "${SPACE}!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
    }
}

###########################################  SET_VAR_INFO_BFM ######################


sub func_STRING  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING: "; }

    my $line  = "";

    $STRING_INDEX_ARRAY{"${type}_${dim}"} = $STRING_INDEX;

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            my $comm  = $param->getComment();
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        my $unitC = ${$param->getComponents()}{$const};
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
                        #my $unitC = ${$param->getComponents()}{$compo};
                        $line .= "${SPACE}var_names($STRING_INDEX)=\"$nameC\"\n";
                        $line .= "${SPACE}var_long($STRING_INDEX)=\"$comm\"\n";
                        $line .= "${SPACE}var_units($STRING_INDEX)=\"$unitC\"\n";
                        $STRING_INDEX++;
                    }
                }
            }else{
                my $unit  = $param->getUnit();
                $line .= "${SPACE}var_names($STRING_INDEX)=\"$name\"\n";
                $line .= "${SPACE}var_long($STRING_INDEX)=\"$comm\"\n";
                $line .= "${SPACE}var_units($STRING_INDEX)=\"$unit\"\n";
                $STRING_INDEX++;
            }
        }
    }

    if( $line ){ print $file $line; }
}


sub func_STRING_DIAGG  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING: "; }

    my $line  = "";
    my $index = 1;

    $STRING_INDEX_ARRAY{"${type}_${dim}"} = $STRING_INDEX;

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $param->getDim() == $dim && $param->getType eq $type ){
            my $group_name = $param->getQuota();
            my $unit       = $param->getUnit();
            my $comm       = $param->getComment();
            #foreach my $name2 (sort keys %$LST_PARAM){
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

    if( $line ){ print $file $line; }
}


sub func_STRING_FIELD  {
    my ( $file, $dim, $type, $spec) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_FIELD: "; }

    my $line;
    
    my $SPEC       = uc($spec);
    my $spec_short = substr($spec,0,3);

    $STRING_INDEX_ARRAY{"${type}_${dim}_${spec}"} = $STRING_INDEX;

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            my $comm  = $param->getComment();
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        my $unitC = decreaseUnit( ${$param->getComponents()}{$const} );
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
                        #my $unitC = decreaseUnit( ${$param->getComponents()}{$compo} );
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

    if( $line ){ print $file $line; }
}


sub func_STRING_INDEX  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "SET_VAR_INFO_BFM -> FUNCTION CALLED func_STRING_INDEX: "; }

    my $line;

    #print Dumper(\%STRING_INDEX_ARRAY) , "\n";
    
    $line .= "${SPACE}stPelStateS=" . $STRING_INDEX_ARRAY{"state_3"}           . "\n";
    $line .= "${SPACE}stPelStateE=" . ( $STRING_INDEX_ARRAY{"diagnos_3"} - 1 ) . "\n";
    $line .= "${SPACE}stPelDiagS="  . $STRING_INDEX_ARRAY{"diagnos_3"}         . "\n";
    $line .= "${SPACE}stPelDiagE="  . ( $STRING_INDEX_ARRAY{"flux_3"} - 1 )    . "\n";
    $line .= "${SPACE}stPelFluxS="  . $STRING_INDEX_ARRAY{"flux_3"}            . "\n";
    $line .= "${SPACE}stPelFluxE="  . ( $STRING_INDEX_ARRAY{"state_2"} - 1 )   ."\n";

    $line .= "${SPACE}stBenStateS=" . $STRING_INDEX_ARRAY{"state_2"}           . "\n";
    $line .= "${SPACE}stBenStateE=" . ( $STRING_INDEX_ARRAY{"diagnos_2"} - 1 ) . "\n";
    $line .= "${SPACE}stBenDiagS="  . $STRING_INDEX_ARRAY{"diagnos_2"}         . "\n";
    $line .= "${SPACE}stBenDiagE="  . ( $STRING_INDEX_ARRAY{"flux_2"} - 1 )    . "\n";
    $line .= "${SPACE}stBenFluxS="  . $STRING_INDEX_ARRAY{"flux_2"}            . "\n";
    $line .= "${SPACE}stBenFluxE="  . ( $STRING_INDEX - 1 )                    . "\n";

    print $file $line;
}


###########################################  INIT_VAR_BFM ######################


sub func_INIT_PP  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_PP: "; }

    my @line_par = ();

    my $line1 = "${SPACE}real(RLEN) :: ";
    my $line2 = "${SPACE}namelist /bfm_init_nml/ ";
    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const . '0';
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo . '0';
                        push( @line_par, ${nameC} );
                    }
                }
            }else{
                push( @line_par, (${name} . '0') );
            }
        }
    }

    if( $#line_par >= 0 ){ 
        printList($file, \@line_par, $line1 ); print $file "\n\n";
        printList($file, \@line_par, $line2 ); print $file "\n";
    }
}


sub func_INIT_DEFAULT  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INIT_DEFAULT: "; }

    my $line = '';

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                #foreach my $compo (keys %{$param->getComponents()} ){
                #$line .= "${SPACE}" . $name . $compo . "0 = _ZERO_\n";
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "INIT_VAR_BFM -> FUNCTION CALLED func_INITS_SET: "; }

    my $line = '';

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        $line .= "${SPACE}  " . $name . $const . " = " . $name . $const . "0\n";
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #$line .= "${SPACE}  " . $name . $compo . " = " . $name . $compo . "0\n";
                    }
                }
            }else{
                $line .= "${SPACE}  " . $name . " = " . $name . "0\n";
            }
        }
    }

    if( $line ){ print $file $line; }
}


###########################################  INCLUDE ######################



sub func_HEADER  {
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_HEADER: "; }

    my $line = '';
    my $TYPE = uc($type);

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
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
    my ( $file, $dim) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_GROUP_HEADER: "; }

    my $line = '';

    #foreach my $name (sort keys %$LST_GROUP){
    foreach my $name ( sort { $$LST_GROUP{$a}->getIndex() cmp $$LST_GROUP{$b}->getIndex() } keys %$LST_GROUP ){
        my $group   = $$LST_GROUP{$name};
        if( $dim == $group->getDim() ){
            $line .= "#define " . $name . "(A,B) D" . $dim . "STATE(pp" . $name . "(A,B),:)\n";
        }
    }

    if( $line ){ print $file $line; }
}


sub func_HEADER_FIELD  {
    my ( $file, $dim, $type, $spec) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_HEADER_FIELD: "; }

    my $line_ini = "";
    my $line_par = "";
    
    my $SPEC       = uc($spec);
    my $spec_short = substr($spec,0,3);

    my $n = 0;
    if( exists $$LST_STA{"${type} ${dim}d ${spec} index"} ){
        $n = $$LST_STA{"${type} ${dim}d ${spec} index"};
    }


    $line_ini .= "#define PEL${SPEC}(A,B) D2DIAGNOS($n+A,B)\n";
    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            if( $param->getComponents() && keys(%{$param->getComponents()}) != 0 ){
                foreach my $const ( sort { $$LST_CONST{$a} cmp $$LST_CONST{$b} } keys %$LST_CONST ){
                    if( exists ${$param->getComponents()}{$const} ){
                        my $nameC = $name . $const;
                        #foreach my $compo (keys %{$param->getComponents()} ){
                        #my $nameC = $name . $compo;
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
    my ( $file, $dim, $type) = @_;
    if ( $VERBOSE ){ print "INCLUDE -> FUNCTION CALLED func_HEADER_DIAGG: "; }

    my $line = '';
    my $TYPE = uc($type);

    #foreach my $name (sort keys %$LST_PARAM){
    foreach my $name ( sort { $$LST_PARAM{$a}->getIndex() cmp $$LST_PARAM{$b}->getIndex() } keys %$LST_PARAM ){
        my $param   = $$LST_PARAM{$name};
        if( $dim == $param->getDim() && $type eq $param->getType() ){
            $line .= "#define " . $name . "(A,B) D" . $dim . "DIAGNOS(pp" . $name . "(A),B)\n";
        }
    }

    if( $line ){ print $file $line; }
}



sub func_HEADER_IF  {}


1;
