#!/usr/bin/perl -w


# DESCRIPTION
#   Classes for configuration generator
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

use strict;
use warnings;

###############Class Test
package Test;

my $DEF_RUN       = 'sh';
my $DEF_MODE      = 'STANDALONE';
my $DEF_MODE_NEMO = 'NEMO.*';
my $DEF_EXE_STD   = 'bfm_standalone.x';
my $DEF_EXE_NEMO  = 'nemo.exe';
my $DEF_RUNPROTO  = 'runscript';
my $DEF_ACTIVE    = 'Y';
my $DEF_TIMMING   = '?';
my $DEF_COMPTHR   = '1e-3';

# order is important in options, has to follow same order as used in "new"
my @OPTIONS= ('ACTIVE', 'PRESET', 'COPY', 'ARCH', 'RUN', 'MODE', 'EXE', 'RUNPROTO', 'RUNEXE', 'QUEUE', 'FORCING', 'VALGRIND', 'PRECMD', 'PRERUN', 'COMPARE', 'COMPTHR', 'PROC', 'PAREXE' );
sub get_options{ my ( $self ) = @_; return @OPTIONS; }

#output codes for compilation, running and analysis
my ($FAIL, $SUCCEED, $NOT_EXE) = ('FAIL', 'SUCCEED', 'NOT_EXE');
sub is_fail   { my ( $self, $value ) = @_; if( $value eq $FAIL)   { return 1; }else{ return 0; } }
sub is_succeed{ my ( $self, $value ) = @_; if( $value eq $SUCCEED){ return 1; }else{ return 0; } }
sub is_not_exe{ my ( $self, $value ) = @_; if( $value eq $NOT_EXE){ return 1; }else{ return 0; } }
sub fail    { my ( $self ) = @_; return $FAIL;    }
sub succeed { my ( $self ) = @_; return $SUCCEED; }
sub not_exe { my ( $self ) = @_; return $NOT_EXE; }
sub is_timmed { my ( $self, $value ) = @_; if( $value ne $DEF_TIMMING ){ return 1; }else{ return 0; } }

#CLASS creator and destructor
sub new{
   my $class = shift;
   my $self = {
       _name     => shift,
       _active   => shift,
       _preset   => shift,
       _copy     => shift,
       _arch     => shift,
       _run      => shift,
       _mode     => shift,
       _exe      => shift,
       _runproto => shift,
       _runexe   => shift,
       _queue    => shift,
       _forcing  => shift,
       _valgrind => shift,
       _precmd   => shift,
       _prerun   => shift,
       _compare  => shift,
       _compthr  => shift,
       _proc     => shift,
       _parexe   => shift,
       _statcmp  => $NOT_EXE,
       _statrun  => $NOT_EXE,
       _statana  => $NOT_EXE,
       _statval  => $NOT_EXE,
       _statcom  => $NOT_EXE,
       _timming  => $DEF_TIMMING,
   };
   return bless $self, $class;   
}
sub DESTROY{}

#GET Functions
sub getName    { my( $self ) = @_; return $self->{_name};   }
sub getActive  { my( $self ) = @_; if($self->{_active}){return $self->{_active};}else{ return $DEF_ACTIVE;} }
sub getPreset  { my( $self ) = @_; return $self->{_preset}; }
sub getCopy    { my( $self ) = @_; return $self->{_copy};   }
sub getArch    { my( $self ) = @_; return $self->{_arch};   }
sub getRun     { my( $self ) = @_; if($self->{_run}) {return $self->{_run}; }else{ return $DEF_RUN; } }
sub getMode    { my( $self ) = @_; if($self->{_mode}){return $self->{_mode};}else{ return $DEF_MODE;} }
sub getExe     { 
    my( $self ) = @_; 
    if( $self->{_exe} ){
        return $self->{_exe};
    }elsif( $self->getMode() =~ m/$DEF_MODE_NEMO/ ){
        return $DEF_EXE_NEMO;
    }else{
        return $DEF_EXE_STD;
    }
}
sub getRunproto{ my( $self ) = @_; return $self->{_runproto}; }
sub getRunexe  { my( $self ) = @_; return $self->{_runexe};   }
sub getQueue   { my( $self ) = @_; return $self->{_queue};    }
sub getForcing { my( $self ) = @_; return $self->{_forcing};  }
sub getValgrind{ my( $self ) = @_; return $self->{_valgrind}; }
sub getPrecmd  { my( $self ) = @_; return $self->{_precmd};   }
sub getPrerun  { my( $self ) = @_; return $self->{_prerun};   }
sub getCompare { my( $self ) = @_; return $self->{_compare};  }
sub getCompthr { my( $self ) = @_; if($self->{_compthr}) {return $self->{_compthr}; }else{ return $DEF_COMPTHR; } }
sub getProc    { my( $self ) = @_; return $self->{_proc};     }
sub getParexe  { 
    my( $self ) = @_; 
    if( $self->{_parexe} =~ /\$PROC|\$\{PROC\}/ && $self->getProc() ){
        my $proc = $self->getProc();
        $self->{_parexe} =~ s/\$PROC|\$\{PROC\}/$proc/;
        return $self->{_parexe};
    }else{ 
        return $self->{_parexe};
    }
}
sub getStatcmp { my( $self ) = @_; return $self->{_statcmp};  }
sub getStatrun { my( $self ) = @_; return $self->{_statrun};  }
sub getStatana { my( $self ) = @_; return $self->{_statana};  }
sub getStatval { my( $self ) = @_; return $self->{_statval};  }
sub getStatcom { my( $self ) = @_; return $self->{_statcom};  }
sub getTimming { my( $self ) = @_; return $self->{_timming};  }

#SET functions
sub setName    { my ( $self, $name      ) = @_; $self->{_name}     = $name     if defined($name);     }
sub setActive  { my ( $self, $active    ) = @_; $self->{_active}   = $active   if defined($active);   }
sub setPreset  { my ( $self, $preset    ) = @_; $self->{_preset}   = $preset   if defined($preset);   }
sub setCopy    { my ( $self, $copy      ) = @_; $self->{_copy}     = $copy     if defined($copy);     }
sub setArch    { my ( $self, $arch      ) = @_; $self->{_arch}     = $arch     if defined($arch);     }
sub setRun     { my ( $self, $run       ) = @_; $self->{_run}      = $run      if defined($run);      }
sub setMode    { my ( $self, $mode      ) = @_; $self->{_mode}     = $mode     if defined($mode);     }
sub setExe     { my ( $self, $exe       ) = @_; $self->{_exe}      = $exe      if defined($exe);      }
sub setRunproto{ my ( $self, $runproto  ) = @_; $self->{_runproto} = $runproto if defined($runproto); }
sub setRunexe  { my ( $self, $runexe    ) = @_; $self->{_runexe}   = $runexe   if defined($runexe);   }
sub setQueue   { my ( $self, $queue     ) = @_; $self->{_queue}    = $queue    if defined($queue);    }
sub setForcing { my ( $self, $forcing   ) = @_; $self->{_forcing}  = $forcing  if defined($forcing);  }
sub setValgrind{ my ( $self, $valgrind  ) = @_; $self->{_valgrind} = $valgrind if defined($valgrind); }
sub setPrecmd  { my ( $self, $precmd    ) = @_; $self->{_precmd}   = $precmd   if defined($precmd);   }
sub setPrerun  { my ( $self, $prerun    ) = @_; $self->{_prerun}   = $prerun   if defined($prerun);   }
sub setCompare { my ( $self, $compare   ) = @_; $self->{_compare}  = $compare  if defined($compare);  }
sub setCompthr { my ( $self, $compthr   ) = @_; $self->{_compthr}  = $compthr  if defined($compthr);  }
sub setProc    { my ( $self, $proc      ) = @_; $self->{_proc}     = $proc     if defined($proc);     }
sub setParexe  { my ( $self, $parexe    ) = @_; $self->{_parexe}   = $parexe   if defined($parexe);   }
sub setStatcmp { my ( $self, $statcmp   ) = @_; $self->{_statcmp}  = $statcmp  if defined($statcmp);  }
sub setStatrun { my ( $self, $statrun   ) = @_; $self->{_statrun}  = $statrun  if defined($statrun);  }
sub setStatana { my ( $self, $statana   ) = @_; $self->{_statana}  = $statana  if defined($statana);  }
sub setStatval { my ( $self, $statval   ) = @_; $self->{_statval}  = $statval  if defined($statval);  }
sub setStatcom { my ( $self, $statcom   ) = @_; $self->{_statcom}  = $statcom  if defined($statcom);  }
sub setTimming { my ( $self, $timming   ) = @_; $self->{_timming}  = $timming  if defined($timming);  }

#PRINT functions
sub print{
    my ( $self ) = @_;
    if( $self->getName()     ){ print "\tTest:     " . $self->getName()       . "\n"; }
    if( $self->getActive()   ){ print "\tActive:   " . $self->getActive()     . "\n"; }
    if( $self->getPreset()   ){ print "\tPreset:   " . $self->getPreset()     . "\n"; }
    if( $self->getCopy()     ){ print "\tCopy:     " . $self->getCopy()       . "\n"; }
    if( $self->getArch()     ){ print "\tArch:     " . $self->getArch()       . "\n"; }
    if( $self->getRun()      ){ print "\tRun:      " . $self->getRun()        . "\n"; }
    if( $self->getMode()     ){ print "\tMode:     " . $self->getMode()       . "\n"; }
    if( $self->getExe()      ){ print "\tExe:      " . $self->getExe()        . "\n"; }
    if( $self->getRunproto() ){ print "\tRunproto: " . $self->getRunproto()   . "\n"; }
    if( $self->getRunexe()   ){ print "\tRunexe  : " . $self->getRunexe  ()   . "\n"; }
    if( $self->getQueue()    ){ print "\tQueue:    " . $self->getQueue()      . "\n"; }
    if( $self->getForcing()  ){ print "\tForcing:  " . $self->getForcing()    . "\n"; }
    if( $self->getValgrind() ){ print "\tValgrind: " . $self->getValgrind()   . "\n"; }
    if( $self->getPrecmd()   ){ print "\tPrecmd:   " . $self->getPrecmd()     . "\n"; }
    if( $self->getPrerun()   ){ print "\tPrerun:   " . $self->getPrerun()     . "\n"; }
    if( $self->getCompare()  ){ print "\tCompare:  " . $self->getCompare()    . "\n"; }
    if( $self->getCompthr()  ){ print "\tCompthr:  " . $self->getCompthr()    . "\n"; }
    if( $self->getProc()     ){ print "\tProc:     " . $self->getProc()       . "\n"; }
    if( $self->getParexe()   ){ print "\tParexe:   " . $self->getParexe()     . "\n"; }
    if( $self->getStatcmp()  ){ print "\tStatcmp:  " . $self->getStatcmp()    . "\n"; }
    if( $self->getStatrun()  ){ print "\tStatrun:  " . $self->getStatrun()    . "\n"; }
    if( $self->getStatana()  ){ print "\tStatana:  " . $self->getStatana()    . "\n"; }
    if( $self->getStatcom()  ){ print "\tStatcom:  " . $self->getStatcom()    . "\n"; }
    if( $self->getStatval()  ){ print "\tStatval:  " . $self->getStatval()    . "\n"; }
    if( $self->getTimming()  ){ print "\tTimming:  " . $self->getTimming()    . "\n"; }
}
sub printAll{
    my ($self, $lst_test) = @_;
    print "Test list: \n";
    foreach my $test (@$lst_test){ $test->print(); print "\n"; }
}
sub printSummary{
    my ($self, $lst_test) = @_;
    printf "%-30s%-13s%-5s%-50s\n", "TEST", "COMPILATION", "RUN", "ANALYSIS";
    foreach my $test (@$lst_test){
        my $test_cmp  = $test->getStatcmp();
        my $test_run  = $test->getStatrun();
        my $test_ana  = $test->getStatana();
        my $test_com  = $test->getStatcom();
        my $test_val  = $test->getStatval();
        
        if   ( $test_cmp eq $FAIL    ) { $test_cmp = 'NO' ; }
        elsif( $test_cmp eq $SUCCEED ) { $test_cmp = 'YES'; }
        elsif( $test_cmp eq $NOT_EXE ) { $test_cmp = '-'  ; }

        if   ( $test_run eq $FAIL    ) { $test_run = 'NO' ; }
        elsif( $test_run eq $SUCCEED ) { $test_run = 'YES'; }
        elsif( $test_run eq $NOT_EXE ) { $test_run = '-'  ; }

        if   ( $test_ana eq $FAIL    ) { $test_ana = 'N'; }
        elsif( $test_ana eq $SUCCEED ) { $test_ana = 'Y'; }
        elsif( $test_ana eq $NOT_EXE ) { $test_ana = '?'; }

        if   ( $test_com eq $FAIL    ) { $test_com = 'N'; }
        elsif( $test_com eq $SUCCEED ) { $test_com = 'Y'; }
        elsif( $test_com eq $NOT_EXE ) { $test_com = '?'; }

        if   ( $test_val eq $FAIL    ) { $test_val = 'N'; }
        elsif( $test_val eq $SUCCEED ) { $test_val = 'Y'; }
        elsif( $test_val eq $NOT_EXE ) { $test_val = '?'; }

        my $test_analysis = "SUCCEED: " . $test_ana . "; "
            . "TIMMING: "    . $test->getTimming() ."; "
            . "COMPARISON: " . $test_com ."; "
            . "VALGRIND: "   . $test_val ."; ";

        printf "%-30s", $test->getName();
        printf "%-13s", $test_cmp;
        printf "%-5s", $test_run;
        printf "%-40s", $test_analysis;
        print "\n";
    }
}


#generate options to execute bfm_configure
sub generate_opt{
    my ( $self ) = @_;
    my $cmd    = '';
    if( $self->getName()     ){ $cmd .= "-x "   . $self->getName()     . " "; }
    if( $self->getPreset()   ){ $cmd .= "-p "   . $self->getPreset()   . " "; }
    if( $self->getArch()     ){ $cmd .= "-a "   . $self->getArch()     . " "; }
    if( $self->getProc()     ){ $cmd .= "-r "   . $self->getProc()     . " "; }
    if( $self->getRunexe()   ){ $cmd .= "-R "   . $self->getRunexe()   . " "; }
    if( $self->getQueue()    ){ $cmd .= "-q "   . $self->getQueue()    . " "; }
    if( $self->getForcing()  ){ $cmd .= "-i \"" . $self->getForcing()  . "\" "; }
    if( $self->getParexe()   ){ $cmd .= "-e \"" . $self->getParexe()   . "\" "; }
    if( $self->getValgrind() ){ $cmd .= "-V \"" . $self->getValgrind() . "\" "; }
    return $cmd;
}

sub generate_scriptName{
    my ( $self ) = @_;

    if( $self->getRunexe() ){
        return $self->getRunexe() . '_' . $self->getName();
    }elsif( $self->getRunproto() ){
        return $self->getRunproto() . '_' . $self->getName();
    }else{
        return $DEF_RUNPROTO . '_' . $self->getName();
    }
}

1;
