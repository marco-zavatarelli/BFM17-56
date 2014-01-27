#!/usr/bin/perl -w

# DESCRIPTION
#   Extract and check configuration file and execute tests
#
# AUTHORS
#   Esteban Gutierrez esteban.gutierrez@cmcc.it
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

package bfm_test_modules;

use 5.008_002;

use strict;
use warnings;
use Exporter;
use Data::Dumper;
use File::Path qw(rmtree);
use File::Copy;

use classes;

########### VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(get_configuration generate_test execute_test analyze_test);
########### VARIABLES ##########################

########### FIX VALUES ##########################
my @OPTIONS      = Test->get_options();
my $BFM_OUT_FILE = 'bfm.out';
########### FIX VALUES ##########################

########### DEFAULT ##########################
my $VERBOSE = 0;
########### DEFAULT ##########################

########### FUNCTIONS ##########################
sub get_configuration{
    my ($conf_file, $build_dir, $verbose) = @_;
    my %user_conf;
    $VERBOSE = $verbose;

    get_configuration_test($conf_file, \%user_conf);
    get_configuration_bfm($build_dir , \%user_conf);
    check_conf(\%user_conf);
    return fill_tests(\%user_conf);
}

sub get_configuration_test{
    my ($conf_file, $user_conf) = @_;

    open CONFIG, "<$conf_file" or die "Could not open configuration file for ${conf_file}: $!";
    while (<CONFIG>) {
        chomp;                  # no newline
        s/#.*//;                # no comments
        s/^\s+//;               # no leading white
        s/\s+$//;               # no trailing white
        next unless length;     # anything left?
        my ($var, $value) = split(/\s*=\s*/, $_, 2);
        my @comasep = split(/,/, "$value");
        $$user_conf{$var} = \@comasep;
    }
    close(CONFIG);
}

sub get_configuration_bfm{
    my ($build_dir, $user_conf) = @_;

    my @modes = ();
    my @exes = ();
    foreach my $preset (@{$$user_conf{PRESET}}){
        my $mode = '';
        my $exe = '';
        my $conf_file = "$build_dir/configurations/$preset/configuration";
        if( ! open CONFIG, "<$conf_file" ){ print "\tProblems opening configuration file \"$conf_file\": $!\n"; next; }
        while (<CONFIG>) {
            chomp;                  # no newline
            s/#.*//;                # no comments
            s/^\s+//;               # no leading white
            s/\s+$//;               # no trailing white
            s/\,$//;                # no coma at the end
            s/\'//g;                # no '
            next unless length;     # anything left?
            my ($var, $value) = split(/\s*=\s*/, $_, 2);
            if( $var =~ 'MODE'   ){ $mode  = $value;  }
            if( $var =~ 'BFMEXE' ){ $exe   = $value;  }
        }
        close CONFIG;
        push( @modes, $mode);
        push( @exes,  $exe );
    }
    $$user_conf{MODE} = \@modes;
    $$user_conf{EXE}  = \@exes;
}

sub check_conf{
    my ($user_conf) = @_;

    if( ! exists $$user_conf{NAME} ){ 
        print "ERROR: Configuration file must contain \"NAME\" variable\n"; 
        exit; 
    }
    my $num_tests = scalar(@{$$user_conf{NAME}});
    foreach my $opt (@OPTIONS){
        if( exists $$user_conf{$opt} && scalar(@{$$user_conf{$opt}}) != $num_tests ){
            print "ERROR: \"$opt\" var must have same number of elements as \"NAME\" var ($num_tests elements)\n"; 
            exit;
        }
    }
}

sub fill_tests{
    my($user_conf) = @_;
    my @lst_test;

    my $num_tests = $#{$$user_conf{NAME}};
    #foreach name create a test
    foreach my $test (0..$num_tests){
        my @lst_opt;
        my $name = ${$$user_conf{NAME}}[$test];
        $name =~ s/^\s+|\s+$//g; #remove leading and trailing spaces
        foreach my $opt (@OPTIONS){
            my $value = '';
            #if exists the option in configuration file, replace the empty value
            if( exists $$user_conf{$opt} ){ 
                $value = ${$$user_conf{$opt}}[$test];
                $value =~ s/^\s+|\s+$//g; #remove leading and trailing spaces
            }
            push(@lst_opt, $value);
        }
        #create test object
        push(@lst_test, new Test($name, @lst_opt ) );
    }
    if($VERBOSE){ Test->printAll(\@lst_test); }
    
    return \@lst_test;
}

sub generate_test{
    my ($build_dir, $bfm_exe, $temp_dir, $test) = @_;

    #remove the output folder
    my $out_dir = "${temp_dir}/" . $test->getName();
    rmtree([$out_dir]);

    #execute the test and capture the output
    my $cmd = "export BFMDIR_RUN=$temp_dir; ";
    $cmd   .= "cd ${build_dir}; ";
    $cmd   .= "$bfm_exe -gcd ";
    $cmd   .= $test->generate_opt();
    if($VERBOSE){ print "\tCommand: $cmd\n"; }
    my $out=`$cmd`;
    
    #check for errors and warnings in generation and compilation time
    if($VERBOSE){
        my @out_warning = $out =~ m/WARNING(?::| )+(.*)/ig;
        if(@out_warning){ print "\tWARNING in ". $test->getName() . ":\n\t\t-" . join("\n\t\t-",@out_warning) . "\n"; }
    }
    my @out_error   = $out =~ m/ERROR(?::| )+(.*)/ig;
    if(@out_error)  { print "\tERROR in "  . $test->getName() . ":\n\t\t-" . join("\n\t\t-",@out_error)   . "\n"; return 0; }
    my @out_exit  = $out =~ m/EXITING\.\.\./ig;
    if(@out_exit)   { print "\tERROR in "  . $test->getName() . ": Compiler not exists\n"; return 0; }

    #copy target simlinks to not to depend from NEMO or BFM current compilation
    # because the current compilation could be overwrite by the next test   
    opendir OUTDIR, "$out_dir" or die "ERROR: reading output directory: $out_dir\n";
    while (my $file = readdir(OUTDIR)) {
        my $file_path = "$out_dir/$file";
        if( -f $file_path && -l $file_path ) {
            my $file_target = readlink($file_path);
            unlink($file_path) or die "Cannot remove symbolic link: $!";
            copy( $file_target, $file_path) or die "Copy failed: $!";            
        }
    }
    close OUTDIR;
    return 1;
}

sub execute_test{
    my ($temp_dir, $test) = @_;

    my $test_dir = "$temp_dir/" . $test->getName();
    if( ! -d $test_dir ){ print "\tERROR in " . $test->getName() . " test dir does not exists: $test_dir\n"; return 0; }
    
    # check if executable exists and put with executable attributes
    my $exe_path = "$test_dir/" . $test->getExe();
    if( ! -f $exe_path  ){ print "\tERROR in " . $test->getName() . " executable does not exists: $exe_path\n"; return 0; }
    `chmod u+x $exe_path`;

    #get the script name
    my $script_name = 'runscript_' . $test->getName();
    my $proc = $test->getPreset();

    #execute in background and return the process name
    if( $test->getRun() ){
        #batch job execution (like using LSF) 
        my $cmd = "cd $test_dir; " . $test->getRun() . " < $script_name";
        if($VERBOSE){ print "\tCommand: $cmd\n"; }
        if ( exec($cmd) == -1 ){ print "\tERROR in " . $test->getName() . " batch execution command fails: $!\n"; return 0; }
        return $test->getPreset();

        if($VERBOSE){ print "\tWaiting for Process Name: $proc"; }
        my $count = 0;
        my $time  = 0;
        do{
            #get the process running with the name of the executable
            $count = scalar grep /$proc/, (split /\n/, `bjobs`);
            if( !$time ){ $time = 50; print "\n\t."; }else{ $time--; print "."; }
            sleep 1;
        }while( $count != 0 );
        print "\n";
    }else{
        #standard shell execution
        my $cmd = "cd $test_dir; chmod u+x $script_name; ./$script_name &> $BFM_OUT_FILE";
        if($VERBOSE){ print "\tCommand: $cmd\n"; }
        if ( exec($cmd) == -1 ){ print "\tERROR in " . $test->getName() . " shell execution command fails: $!\n"; return 0; }
    }
    return 1;
}

sub analyze_test{
    my ($temp_dir, $test) = @_;

    my $test_dir = "$temp_dir/" . $test->getName();

    return 1;
}

1;
