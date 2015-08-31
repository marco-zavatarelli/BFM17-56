#!/usr/bin/perl -w
# merge two namelists, a default one that contains all parameters and another one with only the changed parameters 
# usage: namelist.pl [default_namelist] [changed_namelist] [output_namelist]

#Author: Esteban Gutierrez esteban.gutierrez@cmcc.it

use strict;
use warnings;


#my $XPR_START_BLOCK = '^&(\S*)\s*!(.*)';
#my $XPR_VAR_BLOCK = '\s*(\S*)\s*=\s*(\S*)\s*!(.*)';
my $XPR_START_BLOCK = '^&(\S*)\s*(.*)';
my $XPR_END_BLOCK = '^\/';
my $XPR_VAR_BLOCK = '\s*(\S*)\s*=\s*(\S*)\s*(.*)';


my ($nml_def, $nml_chg, $nml_out);
my $blk_name = "";
my %blk_chg  = ();

sub getArguments{
    if($#ARGV != 2){
        print "usage: $0 [default_namelist] [changed_namelist] [output_namelist]\n";
        print "\t- default_namelist: namelist with default parameters\n";
        print "\t- changed_namelist: namelist with parameters to overwrite in default namelist\n";
        print "\t- output_namelist:  output namelist\n";
        exit;
    }else{
        $nml_def=$ARGV[0];
        $nml_chg=$ARGV[1];
        $nml_out=$ARGV[2];
    }
}


#get the arguments
getArguments();

#open the changed namelist
open NML_CHG, "<", "$nml_chg" or die "$nml_chg cannot be opened: $!";
my @lines_chg = <NML_CHG>;
close(NML_CHG);

#iterate through changed namelist getting the new parameters
$blk_name = "";
foreach my $line_chg (@lines_chg){
    if ( my @blk_start = ( $line_chg =~ /$XPR_START_BLOCK/ ) ){
        #start of block name
        $blk_name = $blk_start[0];
        #print "START of block $blk_name\n";
    }else{
        if( my @blk_end = ( $line_chg =~ /$XPR_END_BLOCK/ ) ){
            #end of block name
            #print "END of block $blk_name\n";
            $blk_name="";
        }else{
            if( my @blk_var = ( $line_chg =~ /$XPR_VAR_BLOCK/ ) ){
                my $var_name    = $blk_var[0];
                my $var_value   = $blk_var[1];
                my $var_comment = $blk_var[2];
                if( $blk_name eq ""){
                    print "WARNING: parameter $var_name not belongs to any namelist block\n";
                }else{
                    push @{ $blk_chg{$blk_name,$var_name} }, ( $var_value, $var_comment );
                    #print "\tvar $var_name = $var_value ! $var_comment\n";
                }
            }
        }
    }
    
}




#open the default file
open NML_DEF, "<", "$nml_def" or die "$nml_def cannot be opened: $!";
my @lines_def = <NML_DEF>;
close(NML_DEF);

#open final file to write
open NML_OUT, ">", "$nml_out" or die "$nml_out cannot be opened: $!";

#iterate through default namelist substituting the parameters for new ones
$blk_name = "";
foreach my $line_def (@lines_def){
    if ( my @blk_start = ( $line_def =~ /$XPR_START_BLOCK/ ) ){
        $blk_name = $blk_start[0];
        print NML_OUT $line_def;
    }else{
        if( my @blk_end = ( $line_def =~ /$XPR_END_BLOCK/ ) ){
            $blk_name="";
            print NML_OUT $line_def;
        }else{
            if( my @blk_var = ( $line_def =~ /$XPR_VAR_BLOCK/ ) ){
                my $var_name    = $blk_var[0];
                my $var_value   = $blk_var[1];
                my $var_comment = $blk_var[2];
                if( exists $blk_chg{$blk_name, $var_name} ){
                    my $var_value   = ${$blk_chg{$blk_name, $var_name}}[0];
                    my $var_comment = ${$blk_chg{$blk_name, $var_name}}[1];
#                    print NML_OUT "\t$var_name = $var_value ! $var_comment\n"
                    print NML_OUT "\t$var_name = $var_value  $var_comment\n"
                }else{
                    print NML_OUT $line_def;
                }
            }else{
                print NML_OUT $line_def;
            }
        }
    }
}
close(NML_OUT);


