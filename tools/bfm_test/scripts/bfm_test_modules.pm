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
########### EXPORT VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(get_configuration generate_test execute_test wait_proc analyze_test);
########### EXPORT VARIABLES ##########################

########### GLOBAL VALUES ##########################
my $BFM_LOG_FILE_EXE = 'bfm.out.exe';
########### GLOBAL VALUES ##########################

########### DEFAULT ##########################
my $VERBOSE = 0;
########### DEFAULT ##########################

########### MODULES ##########################
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
    my $line_next = '';
    foreach my $line (<CONFIG>){
        chomp($line);                    # no newline
        $line =~ s/#.*//;                # no comments
        $line =~ s/^\s+//;               # no leading white
        $line =~ s/\s+$//;               # no trailing white
        if( $line =~ /(.*)\\$/ ){ $line_next .= $1; next; }
        if( $line_next ){ $line = $line_next . $line;  }
        next unless length($line);     # anything left?
        my ($var, $value) = split(/\s*=\s*/, $line, 2);
        my @comasep = split(/\s*,\s*/, "$value");
        $$user_conf{$var} = \@comasep;
        $line_next='';
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
    foreach my $opt ( Test->get_options() ){
        # print "$opt-->" . Dumper($$user_conf{$opt}) . "\n";
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
        foreach my $opt ( Test->get_options()  ){
            my $value = '';
            #if exists the option in configuration file, replace the empty value
            if( exists $$user_conf{$opt} ){ 
                $value = ${$$user_conf{$opt}}[$test];
                $value =~ s/\'//g; # no '
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
    my ($build_dir, $bfm_exe, $overwrite, $temp_dir, $test) = @_;
    my $bfm_log_file_cmp = 'bfm.out.cmp';

    my $test_dir = "${temp_dir}/" . $test->getName();

    #remove the output folder and regenerate if not overwrite
    if( -d $test_dir && !$overwrite ){ $test->setStatcmp(Test->not_exe); return 1; }
    rmtree([$test_dir]);

    #execute the test and capture the output
    my $cmd = "export BFMDIR_RUN=$temp_dir; ";
    if( $test->getPrecmd() ){ $cmd   .= $test->getPrecmd() . "; "; }
    $cmd   .= "cd ${build_dir}; ";
    $cmd   .= "$bfm_exe -gcd ";
    if($VERBOSE){ $cmd .= "-v "; }    
    $cmd   .= $test->generate_opt();
    if($VERBOSE){ print "\tCommand: $cmd\n"; }
    my $out=`$cmd`;

    #save log with compilation in test directory
    if( ! open CMPLOG, ">$test_dir/$bfm_log_file_cmp" ){ print "\tERROR: Problems creating compilation log: $!\n"; $test->setStatcmp(Test->fail); return 0;}
    print CMPLOG $out;
    close CMPLOG;
    if($VERBOSE){ print "\tSee more info in: $test_dir/$bfm_log_file_cmp\n"; }
    
    #check for errors and warnings in generation and compilation time
    if($VERBOSE){
        my @out_warning = $out =~ m/WARNING(?::| )+(.*)/ig;
        if(@out_warning){ print "\tWARNING in ". $test->getName() . ":\n\t\t-" . join("\n\t\t-",@out_warning) . "\n"; }
    }
    my @out_error   = $out =~ m/ERROR(?::| )+(.*)/ig;
    if(@out_error)  { print "\tERROR in "  . $test->getName() . ":\n\t\t-" . join("\n\t\t-",@out_error)   . "\n"; $test->setStatcmp(Test->fail); return 0; }
    my @out_exit  = $out =~ m/EXITING\.\.\./ig;
    if(@out_exit)   { print "\tERROR in "  . $test->getName() . ": Compiler not exists\n"; $test->setStatcmp(Test->fail); return 0; }

    #copy target simlinks to not to depend from NEMO or BFM current compilation
    # because the current compilation could be overwritten by the next test   
    if( ! opendir OUTDIR, "$test_dir" ){ print "\tERROR: reading output directory $test_dir: $!\n"; $test->setStatcmp(Test->fail); return 0; }
    while (my $file = readdir(OUTDIR)) {
        my $file_path = "$test_dir/$file";
        if( -f $file_path && -l $file_path ) {
            my $file_target = readlink($file_path);
            if( ! unlink($file_path) ){ print "\tERROR: Cannot remove symbolic link: $!\n"; $test->setStatcmp(Test->fail); return 0; }
            if( ! copy( $file_target, $file_path) ){ print "\tERROR: Copy failed: $!\n"; $test->setStatcmp(Test->fail); return 0; }
        }
    }

    #if there are forcing dirs, link inside the test directory
    my @forcings = split (';', $test->getForcing());
    foreach my $forcing (@forcings){
        print "\tAdding FORCING dir: $forcing\n";
        my @filelist = glob("$forcing/*");
        foreach my $file (@filelist){
            (my $filename) = $file =~ /.*\/(.*)/;
            if( ! symlink("$file","$test_dir/$filename") ){ print "\tWARNING in ". $test->getName() . ": impossible to link file $file: $!\n"; }
        }
    }

    close OUTDIR;

    $test->setStatcmp(Test->succeed);
    return 1;
}

sub execute_test{
    my ($temp_dir, $test) = @_;

    my $test_dir = "$temp_dir/" . $test->getName();
    if( ! -d $test_dir ){ print "\tERROR in " . $test->getName() . " test dir does not exists: $test_dir\n"; $test->setStatrun(Test->fail); return 0; }
    
    # check if executable exists and put with executable attributes
    my $exe_path = "$test_dir/" . $test->getExe();
    if( ! -f $exe_path  ){ print "\tERROR in " . $test->getName() . " executable does not exists: $exe_path\n"; $test->setStatrun(Test->fail); return 0; }
    `chmod u+x $exe_path`;

    #get the script name
    my $script_name = 'runscript_' . $test->getName();

    #generate the exectuion command
    #generate de change dir
    my $cmd  = "cd $test_dir; ";
    #generate the part to execute before the test command
    if( $test->getPrecmd() ){ $cmd   .= $test->getPrecmd() . "; "; }
    #generate permissions change
    $cmd .= "chmod u+x $script_name; ";
    #generate the test command
    if( $test->getRun() eq 'bsub' ){ #batch jobs needs redirection of input from file 
        $cmd .= $test->getRun() . " < ./$script_name ";
    }else{
        $cmd .= "./$script_name ";
    }
    #generate the loggin part
    $cmd .= " &> $BFM_LOG_FILE_EXE";

    if($VERBOSE){ print "\tCommand: $cmd\n"; }
    if($VERBOSE){ print "\tSee more info in: $test_dir/$BFM_LOG_FILE_EXE\n"; }

    #execute the command
    if ( system($cmd) != 0 ){ print "\tERROR in " . $test->getName() . " cmd execution command fails: $!\n"; $test->setStatrun(Test->fail); return 0; }

    #wait for the command to finish
    wait_proc($test);

    $test->setStatrun(Test->succeed()); 
    return 1;
}

sub wait_proc{
    my ($test) = @_;

    #get command to list process and executable name to wait for
    my $command = 'ps -w';
    my $proc = 'runscript_' . $test->getName();
    if( $test->getRun() eq 'bsub' ){ 
        $command = 'bjobs -w';
        $proc = $test->getPreset();
    }

    my $count = 0;
    # my $time  = 0;
    if($VERBOSE){ print "\tWaiting for Process Name: $proc"; }
    do{
        #get the process running with the name of the executable
        #$count = () = `$command 2> /dev/null` =~ / $proc /;
        $count = scalar grep / $proc /, (split /\n/, `$command 2> /dev/null`);
        # if( !$time ){ $time = 50; print "\n\t."; }else{ $time--; print "."; }
        sleep 1;
    }while( $count != 0 );
    print "\n";
}

sub analyze_test{
    my ($temp_dir, $test) = @_;

    my $timming = '?';

    my $test_dir = "$temp_dir/" . $test->getName();

    #first check if the running finished well
    if($VERBOSE){ print "\tGetting status information\n"; }
    if( $test->getMode() =~ /NEMO/ ){
        #check if time step finished
        my ($out_1, $out_2) = (-1, -2);
        if( open(TIMESTEP, "<", "$test_dir/time.step") ){
            while(<TIMESTEP>){ 
                if( $_ =~ /\s*(\d+)\s*/ ){ $out_1 = $1; }
            }
            close TIMESTEP;
            if($VERBOSE){ print "\t- from file: $test_dir/time.step: $out_1\n"; }
        }
        if( open(NAMELIST, "<", "$test_dir/namelist_cfg") ){
            while(<NAMELIST>){ 
                if( $_ =~ /\s*nn_itend\s*=\s*(\d+)\s*/ ){ $out_2 = $1; } 
            }
            close NAMELIST;
            if($VERBOSE){ print "\t- from file: $test_dir/namelist_cfg: $out_2\n"; }
        }
        if( $out_1 == $out_2){ $test->setStatana(Test->succeed); }
        #check ocean
        if( open(OCEAN, "<", "$test_dir/ocean.output") ){
            if($VERBOSE){ print "\t- from file: $test_dir/ocean.output\n"; }
            while(<OCEAN>){
                if( $_ =~ /"E R R O R"/ ){ $test->setStatana(Test->fail); }
            }
            close OCEAN;
        }
    }else{
        if( open(BFMLOG, "<", "$test_dir/bfm.log") ){ 
            if($VERBOSE){ print "\t- from file: $test_dir/bfm.log\n"; }
            while(<BFMLOG>){ 
                if( $_ =~ /bfm_save : Output saved at the end of the experiment/ ){ $test->setStatana(Test->succeed); }
            }
            close BFMLOG;
        }
    }

    #get timming information
    if($VERBOSE){ print "\tGetting timming information\n"; }
    if( Test->is_succeed($test->getStatana()) ){
        if( $test->getRun() eq 'bsub' ){ # is a batch job
            #get the output file
            my $out_name = "$test_dir/" . $test->getPreset() . "*" . ".out";
            my (@out_files) = glob("$out_name");
            if( $out_files[0] &&  open(JOBEXELOG, "<", "$out_files[0]") ){
                if($VERBOSE){ print "\t- from file: $out_files[0]\n"; }
                #get timming information
                while(<JOBEXELOG>){
                    if( $_ =~ /\s+CPU time\s+:\s+(.+)/ ) { $timming = $1; }
                }
                close(JOBEXELOG);
            }
        }else{ # is a shell execution
            if( open(BFMEXELOG, "<", "$test_dir/$BFM_LOG_FILE_EXE") ){
                if($VERBOSE){ print "\t- from file: $test_dir/$BFM_LOG_FILE_EXE\n"; }
                while(<BFMEXELOG>){
                    if($_ =~ /^real\s*(.*)/ ){ $timming = $1; }
                }
                close BFMEXELOG;
            }
        }
    }

    #compare results if compare dir exists
    if( $test->getCompare ){
        if($VERBOSE){ print "\tComparing results\n"; }
    }

    #check valgring options if valgrind exists
    if( $test->getValgrind ){
        if($VERBOSE){ print "\tComparing results\n"; }
    }

    #generate summary with timing, memmory and results
    my $status = 'N';
    if( Test->is_succeed($test->getStatana()) ){ $status = 'Y'; }
    $test->setResult("SUCCEED: $status TIMMING: $timming");

    return 1;
}
########### MODULES ##########################

1;
