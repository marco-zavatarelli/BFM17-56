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

sub new{
   my $class = shift;
   my $self = {
       _name    => shift,
       _forcing => shift,
   };
   return bless $self, $class;   
}

sub DESTROY{}

sub getName   { my( $self ) = @_; return $self->{_name}   ; }
sub getForcing{ my( $self ) = @_; return $self->{_forcing}; }

sub setName   { my ( $self, $name    ) = @_; $self->{_name}    = $name    if defined($name)   ; }
sub setForcing{ my ( $self, $forcing ) = @_; $self->{_forcing} = $forcing if defined($forcing); }

sub print{
    my ( $self ) = @_;
    print "Test ";
    if($self->{_name})    { print "Name: "    . $self->getName()    . "; "; }
    if($self->{_forcing}) { print "Forcing: " . $self->getForcing() . "; "; }
}

sub printAll{
    my ($self, $lst_test) = @_;
    print "Test list: \n";
    foreach my $test (@$lst_test){ print "\t" ; $test->print(); print "\n"; }
}


1;
