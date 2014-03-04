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
use File::Basename;

use classes;
########### EXPORT VARIABLES ##########################
our @ISA = qw(Exporter);
our @EXPORT= qw(get_configuration generate_test execute_test wait_proc analyze_test);
########### EXPORT VARIABLES ##########################

########### GLOBAL VALUES ##########################
my $BFM_LOG_FILE_CMP = 'bfm.out.cmp';
my $BFM_LOG_FILE_EXE = 'bfm.out.exe';
my $BFM_LOG_FILE_COM = 'bfm.out.com';
my $BFM_LOG_FILE_MEM = 'bfm.out.mem';
my $BFM_LOG_FILE_CAC = 'bfm.out.cac';
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
    
    return \@lst_test;
}

sub generate_test{
    my ($build_dir, $bfm_exe, $overwrite, $temp_dir, $test) = @_;

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
    if( ! open CMPLOG, ">$test_dir/$BFM_LOG_FILE_CMP" ){ print "\tERROR: Problems creating compilation log: $!\n"; $test->setStatcmp(Test->fail); return 0;}
    print CMPLOG $cmd . "\n";
    print CMPLOG $out;
    close CMPLOG;
    if($VERBOSE){ print "\tSee more info in: $test_dir/$BFM_LOG_FILE_CMP\n"; }
    
    #check for errors and warnings in generation and compilation time
    if($VERBOSE){
        my @out_warning = $out =~ m/WARNING(?::| )+(.*)/ig;
        if(@out_warning){ print "\tWARNING in ". $test->getName() . ":\n\t\t-" . join("\n\t\t- ",@out_warning) . "\n"; }
    }
    my @out_error   = $out =~ m/ERROR(?::| )+(.*)/ig;
    if(@out_error)  { print "\tERROR in "  . $test->getName() . ":\n\t\t-" . join("\n\t\t- ",@out_error)   . "\n"; $test->setStatcmp(Test->fail); return 0; }
    my @out_exit  = $out =~ m/EXITING\.\.\./ig;
    if(@out_exit)   { print "\tERROR in "  . $test->getName() . ": Compiler not exists\n"; $test->setStatcmp(Test->fail); return 0; }

    #copy executable simlink to not to depend from NEMO or BFM current compilation
    # because the current compilation could be overwritten by the next test   
    my $file_exe = "$test_dir/" . $test->getExe();
    if( -f $file_exe && -l $file_exe ) {
        my $file_target = readlink($file_exe);
        if( ! unlink($file_exe) ){ print "\tERROR: Cannot remove symbolic link: $!\n"; $test->setStatcmp(Test->fail); return 0; }
        if( ! copy( $file_target, $file_exe) ){ print "\tERROR: Copy failed: $!\n"; $test->setStatcmp(Test->fail); return 0; }
    }

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
    $cmd .= " >>$test_dir/$BFM_LOG_FILE_EXE 2>&1";

    if($VERBOSE){ print "\tCommand: $cmd\n"; }
    #save log with execution in test directory
    if( ! open EXELOG, ">$test_dir/$BFM_LOG_FILE_EXE" ){ print "\tERROR: Problems creating execution log: $!\n"; $test->setStatrun(Test->fail); return 0;}
    print EXELOG $cmd . "\n";
    close EXELOG;
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
        # if( !$time ){ $time = 50; if($VERBOSE){ print "\n\t."; } }else{ $time--; if($VERBOSE){ print "."; } }
        sleep 1;
    }while( $count != 0 );
    if($VERBOSE){ print "\n"; }
}

sub analyze_test{
    my ($temp_dir, $test) = @_;

    my $test_dir = "$temp_dir/" . $test->getName();

    #get status information
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
        my $file_name = "$test_dir/namelist_cfg";
        if( -e "$test_dir/namelist" ){ $file_name = "$test_dir/namelist" } # in NEMO 3.4 this is the name of the namelist
        if( open(NAMELIST, "<", "$file_name") ){
            while(<NAMELIST>){ 
                if( $_ =~ /\s*nn_itend\s*=\s*(\d+)\s*/ ){ $out_2 = $1; } 
            }
            close NAMELIST;
            if($VERBOSE){ print "\t- from file: $file_name: $out_2\n"; }
        }
        if( $out_1 == $out_2){ $test->setStatana(Test->succeed); }
        else{ $test->setStatana(Test->fail); }

        #check ocean
        if( open(OCEAN, "<", "$test_dir/ocean.output") ){
            if($VERBOSE){ print "\t- from file: $test_dir/ocean.output\n"; }
            while(<OCEAN>){
                if( $_ =~ /"E R R O R"/ ){ $test->setStatana(Test->fail); last; }
            }
            close OCEAN;
        }
    }else{
        #by default in shell mode if output is not present result is fail
        $test->setStatana(Test->fail);
        if( $test->getRun() eq 'bsub' ){ # is a batch job
            #get the error output file
            my $err_name = "$test_dir/" . $test->getPreset() . "*" . ".err";
            my (@err_files) = glob("$err_name");
            if( $err_files[0] &&  open(JOBERRLOG, "<", "$err_files[0]") ){
                if($VERBOSE){ print "\t- from file: $err_files[0]\n"; }
                #get end log
                while(<JOBERRLOG>){
                    if( $_ =~ /BFM standalone finished on/ ) { $test->setStatana(Test->succeed); last; }
                }
                close(JOBERRLOG);
            }
        }else{ # is a shell job
            if( open(BFMEXELOG, "<", "$test_dir/$BFM_LOG_FILE_EXE") ){ 
                if($VERBOSE){ print "\t- from file: $test_dir/$BFM_LOG_FILE_EXE\n"; }
                while(<BFMEXELOG>){ 
                    if( $_ =~ /BFM standalone finished on/ ) { $test->setStatana(Test->succeed); last; }
                }
                close BFMEXELOG;
            }
        }
    }
    #return if test fail
    if( Test->is_fail($test->getStatana()) ){ if($VERBOSE){ print "\t- Status FAIL\n"; } return 0; }
    if( Test->is_succeed($test->getStatana()) && $VERBOSE ){ print "\t- Status OK\n"; }

    #get timming information
    if($VERBOSE){ print "\tGetting timming information\n"; }
    if( $test->getRun() eq 'bsub' ){ # is a batch job
        #get the output err file
        my $err_name = "$test_dir/" . $test->getPreset() . "*" . ".err";
        my (@err_files) = glob("$err_name");
        if( $err_files[0] &&  open(JOBEXELOG, "<", "$err_files[0]") ){
            if($VERBOSE){ print "\t- from file: $err_files[0]\n"; }
            #get timming information
            while(<JOBEXELOG>){
                if( $_ =~ /^real\s+(.+)/ ) { $test->setTimming($1); }
            }
            close(JOBEXELOG);
        }
    }else{ # is a shell execution
        if( open(BFMEXELOG, "<", "$test_dir/$BFM_LOG_FILE_EXE") ){
            if($VERBOSE){ print "\t- from file: $test_dir/$BFM_LOG_FILE_EXE\n"; }
            while(<BFMEXELOG>){
                if($_ =~ /^real\s*(.*)/ ){ $test->setTimming($1); }
            }
            close BFMEXELOG;
        }
    }
    if( ! Test->is_timmed($test) ){
        if($VERBOSE){ print "\t- No timming information present\n"; }
    }elsif( $VERBOSE ){ print "\t- Timming: " . $test->getTimming() . "\n"; }

    #get comapre information
    if( $test->getCompare() ){
        my $cmp_dir = $test->getCompare();
        if($VERBOSE){ print "\tComparing results:\n"; }
        if($VERBOSE){ print "\t- $test_dir VS $cmp_dir\n"; }
        if($VERBOSE){ print "\t- See more info in: $test_dir/$BFM_LOG_FILE_COM\n"; }

        #check if command exists
        my $out_cmp = `nccmp -V 2>&1`;
        if( $out_cmp =~ /nccmp \d.\d.\d/){
            my @filelist_test = ();
            my @filelist_cmp  = ();
            my $num_proc = $test->getProc();

            #get all outptut files from BFM in the input dir
            if( open(NAMELIST, "<", "$test_dir/BFM_General.nml") ){
                my $out_bfm = '';
                while(<NAMELIST>){ 
                    if( $_ =~ /\s*out_fname\s*=\s*\"*\'*(\w+)\"*\'*\s*/ ){ $out_bfm = $1; }
                }
                close NAMELIST;
                if( $out_bfm ){
                    if($VERBOSE){ print "\t- from file: $test_dir/BFM_General.nml: $out_bfm\n"; }
                    my @temp_list = glob( "$test_dir/${out_bfm}_*.nc" );
                    foreach my $file ( @temp_list ){ if( $file =~ /${out_bfm}_\d+.nc/ ) {push( @filelist_test, $file ); } }
                }
            }
            #get all outptut files from BFM in the compare dir
            if( open(NAMELIST, "<", "$cmp_dir/BFM_General.nml") ){
                my $out_bfm = '';
                while(<NAMELIST>){ 
                    if( $_ =~ /\s*out_fname\s*=\s*\"*\'*(\w+)\"*\'*\s*/ ){ $out_bfm = $1; }
                }
                close NAMELIST;
                if( $out_bfm ){
                    if($VERBOSE){ print "\t- from file: $cmp_dir/BFM_General.nml: $out_bfm\n"; }
                    my @temp_list = glob( "$cmp_dir/${out_bfm}_*.nc" );
                    foreach my $file ( @temp_list ){ if( $file =~ /${out_bfm}_\d+.nc/ ) {push( @filelist_cmp, $file ); } }
                }
            }

            #check if are the same or not null
            if( ($#filelist_test != $#filelist_cmp) || $#filelist_test < 1 || $#filelist_cmp < 1 ){
                if($VERBOSE){ 
                    print "\tWARNING in ". $test->getName() . ": Comparison impossible, not same or null number of NETCDF files: $#filelist_test <> $#filelist_cmp.\n";
                }
                $test->setStatcom(Test->fail);
            }else{
                #execute nccmp
                my %filelist_test = map {basename($_)=>1} @filelist_test;
                my @only_in_cmp   = grep { !$filelist_test{basename($_)} } @filelist_cmp;
                my %filelist_cmp  = map {basename($_)=>1} @filelist_cmp;
                my @only_in_test  = grep { !$filelist_cmp{basename($_)} } @filelist_test;
                if( ! open COMLOG, ">$test_dir/$BFM_LOG_FILE_COM" ){ print "\tERROR: Problems creating comparison log: $!\n"; $test->setStatcom(Test->fail); return 0;}
                print COMLOG "COMPARISON:\n";
                close COMLOG;
                foreach my $name (keys %filelist_test){
                    my $cmd = "nccmp -dgmf -t 1e-3 $test_dir/$name $cmp_dir/$name";
                    `echo $cmd >>$test_dir/$BFM_LOG_FILE_COM`;
                    `$cmd >>$test_dir/$BFM_LOG_FILE_COM 2>&1`;
                }
                #analyze the output
                if( open(BFMCOMLOG, "<", "$test_dir/$BFM_LOG_FILE_COM") ){
                    my %error_vars;
                    my $name_last = "";
                    my $var_last  = "";
                    my $name_ins  = "";
                    my $var_ins   = "";
                    while(<BFMCOMLOG>){
                        if( $_ =~ /^nccmp \-dgmf \-t 1e-3 (.*) / ){ #get the name of the last file analyzed 
                            $name_last = basename($1);
                        }
                        if( $_ =~ /^DIFFER :[^:]+: ([^:]+) :/ ){ #get the variable which differs
                            $var_last = $1;
                            if( $name_ins ne $name_last || $var_ins ne $var_last ){
                                #if not inserted before, push to the summary hash
                                push( @{$error_vars{$var_last}}, $name_last );
                                $var_ins  = $var_last;
                                $name_ins = $name_last;
                            }
                        }
                    }
                    close BFMCOMLOG;
                    #print Dumper(\%error_vars) . "\n";
                    #check if variables differs (not timestamp)
                    my $count_vars = keys %error_vars;
                    my @list = keys %error_vars;
                    if( exists $error_vars{timeStamp} ){ $count_vars -= 1; }
                    if( $count_vars == 0 ){ $test->setStatcom(Test->succeed); if($VERBOSE){ print "\t- Comparison OK\n";   } }
                    else{                   $test->setStatcom(Test->fail);    if($VERBOSE){ print "\t- Comparison FAIL. Variables differs: @list\n"; } }
                    #insert results at the end of comparison log
                    if( open(BFMCOMLOG, ">>", "$test_dir/$BFM_LOG_FILE_COM") ){
                        print BFMCOMLOG "\n\nSUMMARY:\n\t- $count_vars differs: @list \n" . Dumper(\%error_vars) . "\n";
                        close BFMCOMLOG;
                    }else{
                        if($VERBOSE){ 
                            print "\tWARNING in ". $test->getName() . ": Cannot open comparison file $test_dir/$BFM_LOG_FILE_COM to store summary data: $!.\n";
                        }
                        $test->setStatcom(Test->fail);
                    }
                }else{
                    if($VERBOSE){ 
                        print "\tWARNING in ". $test->getName() . ": Cannot open comparison file $test_dir/$BFM_LOG_FILE_COM: $!.\n";
                    }
                    $test->setStatcom(Test->fail);
                }
            }
        }else{ 
            if($VERBOSE){ 
                print "\tWARNING in ". $test->getName() . ": Comparison not executed because \"nccmp\" command does not exist or is not executable.\n";
            }
            $test->setStatcom(Test->fail);
        }
    }

    #get Valgrind information
    if( $test->getValgrind ){
        if($VERBOSE){ print "\tGetting Valgrind results\n"; }

        #get the massif output
        if( $test->getValgrind =~ /--tool=massif/ ){
            if($VERBOSE){ print "\t- Getting Massif results\n"; }
            my $mass_name = "$test_dir/massif.out.*";
            my (@mass_files) = glob("$mass_name");
            if( $mass_files[0] ){
                if($VERBOSE){ print "\t- memory summary in files: $test_dir/$BFM_LOG_FILE_MEM\n"; }
                my $cmd = "ms_print $mass_files[0] >>$test_dir/$BFM_LOG_FILE_MEM 2>&1";
                if( ! open MEMLOG, ">$test_dir/$BFM_LOG_FILE_MEM" ){ print "\tERROR: Problems creating memory log: $!\n"; $test->setStatval(Test->fail); return 0;}
                print MEMLOG "$cmd\n";
                `$cmd`;
                close MEMLOG;
            }elsif($VERBOSE){ print "\t- Massif FAIL: output not present\n"; }
        }

        #get the cachegrind output
        if( $test->getValgrind =~ /--tool=cachegrind/ ){
            if($VERBOSE){ print "\t- Getting Cachegrind results\n"; }
            my $cache_name  = "$test_dir/cachegrind.out.*";
            my (@cache_files) = glob("$cache_name");
            #if there is output
            if( $cache_files[0] ){
                if($VERBOSE){ print "\t- cache summary in file: $test_dir/$BFM_LOG_FILE_CAC\n"; }
                my $cache_merge = "$test_dir/cachegrind.merge";
                my $cmd_cachemerge = "cg_merge -o $cache_merge " . join(' ',@cache_files) . " >>$test_dir/$BFM_LOG_FILE_CAC 2>&1";
                my $cmd_cacheanno ="cg_annotate --auto=yes $cache_merge >>$test_dir/$BFM_LOG_FILE_CAC 2>&1";
                #sum all outputs
                if( ! open CACLOG, ">$test_dir/$BFM_LOG_FILE_CAC" ){ print "\tERROR: Problems creating cache log for merge: $!\n"; $test->setStatval(Test->fail); return 0;}
                print CACLOG "MERGE:\n$cmd_cachemerge\n";
                close CACLOG;
                `$cmd_cachemerge`;
                #get the summary
                if( ! open CACLOG, ">>$test_dir/$BFM_LOG_FILE_CAC" ){ print "\tERROR: Problems openning cache log for annotate: $!\n"; $test->setStatval(Test->fail); return 0;}
                print CACLOG "\n\nANNOTATE:\n$cmd_cacheanno\n";
                close CACLOG;
                `$cmd_cacheanno`;
            }elsif($VERBOSE){ print "\t- Cachegrind FAIL: output not present\n"; }
        }

        #get the memcheck output
        if( $test->getValgrind =~ /--tool=memcheck/ ){
            if($VERBOSE){ print "\t- Getting Memcheck results\n"; }
            my $val_name = "$test_dir/valgrind_*";
            my (@val_files) = glob("$val_name");
            if( $val_files[0] ){
                if($VERBOSE){ print "\t- memory summary in file: @val_files\n"; }
            }elsif($VERBOSE){ print "\t- Memcheck FAIL: output not present\n"; }
        }

        $test->setStatval(Test->succeed);
    }
    return 1;
}
########### MODULES ##########################

1;
