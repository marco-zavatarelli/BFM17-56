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

my $DEF_MODE     = 'STANDALONE';
my $DEF_EXE_STD  = 'bfm_standalone.x';
my $DEF_EXE_NEMO = 'nemo.exe';
# order is important in options, has to follow same order as used in "new"
my @OPTIONS= ('PRESET', 'ARCH', 'RUN', 'MODE', 'EXE', 'FORCING', 'VALGRIND', 'PRECMD');
sub get_options{ my ( $self ) = @_; return @OPTIONS; }

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
   };
   return bless $self, $class;   
}

sub DESTROY{}

sub getName    { my( $self ) = @_; return $self->{_name}    ; }
sub getPreset  { my( $self ) = @_; return $self->{_preset}  ; }
sub getArch    { my( $self ) = @_; return $self->{_arch}    ; }
sub getRun     { my( $self ) = @_; return $self->{_run}     ; }
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
sub getForcing { my( $self ) = @_; return $self->{_forcing} ; }
sub getValgrind{ my( $self ) = @_; return $self->{_valgrind}; }
sub getPrecmd  { my( $self ) = @_; return $self->{_precmd}  ; }

sub setName    { my ( $self, $name      ) = @_; $self->{_name}    = $name      if defined($name)    ; }
sub setPreset  { my ( $self, $preset    ) = @_; $self->{_preset}  = $preset    if defined($preset)  ; }
sub setArch    { my ( $self, $arch      ) = @_; $self->{_arch}    = $arch      if defined($arch)    ; }
sub setRun     { my ( $self, $run       ) = @_; $self->{_run}     = $run       if defined($run)     ; }
sub setMode    { my ( $self, $mode      ) = @_; $self->{_mode}    = $mode      if defined($mode)    ; }
sub setExe     { my ( $self, $exe       ) = @_; $self->{_exe}     = $exe       if defined($exe)     ; }
sub setForcing { my ( $self, $forcing   ) = @_; $self->{_forcing} = $forcing   if defined($forcing) ; }
sub setValgrind{ my ( $self, $valgrind  ) = @_; $self->{_valgrind} = $valgrind if defined($valgrind); }
sub setPrecmd  { my ( $self, $precmd    ) = @_; $self->{_precmd}  = $precmd    if defined($precmd)  ; }


sub print{
    my ( $self ) = @_;
    print "Test ";
    if($self->{_name})     { print "\tName:     " . $self->getName()       . ";\n"; }
    if($self->{_preset})   { print "\tPreset:   " . $self->getPreset()     . ";\n"; }
    if($self->{_arch})     { print "\tArch:     " . $self->getArch()       . ";\n"; }
    if($self->{_run})      { print "\tRun:      " . $self->getRun()        . ";\n"; }
    if($self->{_exe})      { print "\tExe:      " . $self->getExe()        . ";\n"; }
    if($self->{_mode})     { print "\tMode:     " . $self->getMode()       . ";\n"; }
    if($self->{_forcing})  { print "\tForcing:  " . $self->getForcing()    . ";\n"; }
    if($self->{_valgrind}) { print "\tValgrind: " . $self->getValgrind()   . ";\n"; }
    if($self->{_precmd})   { print "\tPrecmd:   " . $self->getPrecmd()     . ";\n"; }
}

#generate options to execute bfm_configure
sub generate_opt{
    my ( $self ) = @_;
    my $cmd = '';
    if($self->{_name})   { $cmd .= " -x "  . $self->getName()      . " "; }
    if($self->{_preset}) { $cmd .= " -p "  . $self->getPreset()    . " "; }
    if($self->{_arch})   { $cmd .= " -a "  . $self->getArch()      . " "; }
    return $cmd;
}

sub printAll{
    my ($self, $lst_test) = @_;
    print "Test list: \n";
    foreach my $test (@$lst_test){ print "\t" ; $test->print(); print "\n"; }
}


1;
