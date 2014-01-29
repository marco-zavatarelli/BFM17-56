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

use Data::Dumper;

###############Class Test
package Test;

my $DEF_RUN      = 'sh';
my $DEF_MODE     = 'STANDALONE';
my $DEF_EXE_STD  = 'bfm_standalone.x';
my $DEF_EXE_NEMO = 'nemo.exe';
# order is important in options, has to follow same order as used in "new"
my @OPTIONS= ('PRESET', 'ARCH', 'RUN', 'MODE', 'EXE', 'FORCING', 'VALGRIND', 'PRECMD', 'COMPARE', 'PROC', 'PAREXE' );
sub get_options{ my ( $self ) = @_; return @OPTIONS; }

#CLASS creator and destructor
sub new{
   my $class = shift;
   my $self = {
       _name     => shift,
       _preset   => shift,
       _arch     => shift,
       _run      => shift,
       _mode     => shift,
       _exe      => shift,
       _forcing  => shift,
       _valgrind => shift,
       _precmd   => shift,
       _compare  => shift,
       _proc     => shift,
       _parexe   => shift,
   };
   return bless $self, $class;   
}
sub DESTROY{}

#GET Functions
sub getName    { my( $self ) = @_; return $self->{_name};   }
sub getPreset  { my( $self ) = @_; return $self->{_preset}; }
sub getArch    { my( $self ) = @_; return $self->{_arch};   }
sub getRun     { my( $self ) = @_; if($self->{_run}) {return $self->{_run}; }else{ return $DEF_RUN; } }
sub getMode    { my( $self ) = @_; if($self->{_mode}){return $self->{_mode};}else{ return $DEF_MODE;} }
sub getExe     { 
    my( $self ) = @_; 
    if( $self->{_exe} ){
        return $self->{_exe};
    }elsif( $self->{_mode} =~ m/NEMO.*/ ){
        return $DEF_EXE_NEMO;
    }else{
        return $DEF_EXE_STD;
    }
}
sub getForcing { my( $self ) = @_; return $self->{_forcing};  }
sub getValgrind{ my( $self ) = @_; return $self->{_valgrind}; }
sub getPrecmd  { my( $self ) = @_; return $self->{_precmd};   }
sub getCompare { my( $self ) = @_; return $self->{_compare};  }
sub getProc    { my( $self ) = @_; return $self->{_proc};     }
sub getParexe  { 
    my( $self ) = @_; 
    if( $self->{_parexe} =~ /\$PROC|\$\{PROC\}/ && $self->{_proc} ){
        my $proc = $self->{_proc};
        $self->{_parexe} =~ s/\$PROC|\$\{PROC\}/$proc/;
        return $self->{_parexe};
    }else{ 
        return $self->{_parexe};
    }
}
sub getResult  { my( $self ) = @_; return $self->{_result};   }

#SET functions
sub setName    { my ( $self, $name      ) = @_; $self->{_name}     = $name     if defined($name);     }
sub setPreset  { my ( $self, $preset    ) = @_; $self->{_preset}   = $preset   if defined($preset);   }
sub setArch    { my ( $self, $arch      ) = @_; $self->{_arch}     = $arch     if defined($arch);     }
sub setRun     { my ( $self, $run       ) = @_; $self->{_run}      = $run      if defined($run);      }
sub setMode    { my ( $self, $mode      ) = @_; $self->{_mode}     = $mode     if defined($mode);     }
sub setExe     { my ( $self, $exe       ) = @_; $self->{_exe}      = $exe      if defined($exe);      }
sub setForcing { my ( $self, $forcing   ) = @_; $self->{_forcing}  = $forcing  if defined($forcing);  }
sub setValgrind{ my ( $self, $valgrind  ) = @_; $self->{_valgrind} = $valgrind if defined($valgrind); }
sub setPrecmd  { my ( $self, $precmd    ) = @_; $self->{_precmd}   = $precmd   if defined($precmd);   }
sub setCompare { my ( $self, $compare   ) = @_; $self->{_compare}  = $compare  if defined($compare);  }
sub setProc    { my ( $self, $proc      ) = @_; $self->{_proc}     = $proc     if defined($proc);     }
sub setParexe  { my ( $self, $parexe    ) = @_; $self->{_parexe}   = $parexe   if defined($parexe);   }
sub setResult  { my ( $self, $result    ) = @_; $self->{_result}   = $result   if defined($result);   }

#PRINT functions
sub print{
    my ( $self ) = @_;
    if($self->{_name})     { print "\tTest:     " . $self->getName()       . "\n"; }
    if($self->{_preset})   { print "\tPreset:   " . $self->getPreset()     . "\n"; }
    if($self->{_arch})     { print "\tArch:     " . $self->getArch()       . "\n"; }
    if($self->{_run})      { print "\tRun:      " . $self->getRun()        . "\n"; }
    if($self->{_exe})      { print "\tExe:      " . $self->getExe()        . "\n"; }
    if($self->{_mode})     { print "\tMode:     " . $self->getMode()       . "\n"; }
    if($self->{_forcing})  { print "\tForcing:  " . $self->getForcing()    . "\n"; }
    if($self->{_valgrind}) { print "\tValgrind: " . $self->getValgrind()   . "\n"; }
    if($self->{_precmd})   { print "\tPrecmd:   " . $self->getPrecmd()     . "\n"; }
    if($self->{_compare})  { print "\tCompare:  " . $self->getCompare()    . "\n"; }
    if($self->{_proc})     { print "\tProc:     " . $self->getProc()       . "\n"; }
    if($self->{_parexe})   { print "\tParexe:   " . $self->getParexe()     . "\n"; }
    if($self->{_result})   { print "\tResult:   " . $self->getResult()     . "\n"; }
}
sub printAll{
    my ($self, $lst_test) = @_;
    print "Test list: \n";
    foreach my $test (@$lst_test){ $test->print(); print "\n"; }
}

#generate options to execute bfm_configure
sub generate_opt{
    my ( $self ) = @_;
    my $cmd = '';
    if($self->{_name})   { $cmd .= "-x "  . $self->getName()   . " "; }
    if($self->{_preset}) { $cmd .= "-p "  . $self->getPreset() . " "; }
    if($self->{_arch})   { $cmd .= "-a "  . $self->getArch()   . " "; }
    if($self->{_proc})   { $cmd .= "-r "  . $self->getProc()   . " "; }
    if($self->{_parexe}) { $cmd .= "-e "  . $self->getParexe() . " "; }
    return $cmd;
}

1;
