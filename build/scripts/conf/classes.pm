#!/usr/bin/perl -w

#Author: Esteban Gutierrez esteban.gutierrez@cmcc.it

use strict;
use warnings;

use Data::Dumper;

###############Class Parameter
package Parameter;

sub new{
   my $class = shift;
   my $self = {
       _sigla    => shift,
       _dim      => shift,
       _type     => shift,
       _include  => shift,
       _z        => shift,
       _unit     => shift,
       _compo    => shift,
       _compo_ex => shift,
       _comm     => shift,
       _func     => shift,
       _group    => shift,
       _quota    => shift,
   };
   return bless $self, $class;   
}

sub DESTROY{}

sub getSigla{ my( $self ) = @_; return $self->{_sigla}; }
sub getDim{ my( $self ) = @_; return $self->{_dim}; }
sub getType{ my( $self ) = @_; return $self->{_type}; }
sub getInclude{ my( $self ) = @_; return $self->{_include}; }
sub getZ{ my( $self ) = @_; return $self->{_z}; }
sub getUnit{ my( $self ) = @_; return $self->{_unit}; }
sub getComponents{ my( $self ) = @_; return $self->{_compo}; }
sub getComponentsEx{ my( $self ) = @_; return $self->{_compo_ex}; }
sub getComment{ my( $self ) = @_; return $self->{_comm}; }
sub getFunction{ my( $self ) = @_; return $self->{_func}; }
sub getGroup{ my( $self ) = @_; return $self->{_group}; }
sub getQuota{ my( $self ) = @_; return $self->{_quota}; }

sub setSigla{ my ( $self, $sigla ) = @_; $self->{_sigla} = $sigla if defined($sigla);}
sub setDim{ my ( $self, $dim ) = @_; $self->{_dim} = $dim if defined($dim);}
sub setType{ my ( $self, $type ) = @_; $self->{_type} = $type if defined($type);}
sub setInclude{ my ( $self, $include ) = @_; $self->{_include} = $include if defined($include);}
sub setZ{ my ( $self, $z ) = @_; $self->{_z} = $z if defined($z);}
sub setUnit{ my ( $self, $unit ) = @_; $self->{_unit} = $unit if defined($unit);}
sub setComponents{ my ( $self, $compo ) = @_; $self->{_compo} = $compo if defined($compo);}
sub setComponentsEx{ my ( $self, $compo_ex ) = @_; $self->{_compo_ex} = $compo_ex if defined($compo_ex);}
sub setComment{ my ( $self, $comm ) = @_; $self->{_comm} = $comm if defined($comm);}
sub setFunction{ my ( $self, $func ) = @_; $self->{_func} = $func if defined($func);}
sub setGroup{ my ( $self, $group ) = @_; $self->{_group} = $group if defined($group);}
sub setQuota{ my ( $self, $quota ) = @_; $self->{_quota} = $quota if defined($quota);}

sub print{
    my ( $self ) = @_;
    print "PARAMETER => ";
    if($self->{_sigla})    { print              $self->getSigla() . "; ";    }
    if($self->{_dim})      { print "Dim: "   .  $self->getDim() . "; ";      }
    if($self->{_type})     { print "Type: "  .  $self->getType() . "; ";     }
    if($self->{_include})  { print "Incl: "  .  $self->getInclude() . "; ";  }
    if($self->{_z})        { print "-Z: "    .  $self->getZ() . "; ";        }
    if($self->{_unit})     { print "Unit: "  .  $self->getUnit() . "; ";     }
    if($self->{_compo})    { print "Units:" ;   while ( my ($k,$v) = each %{$self->getComponents()}   ) { print " $k=>$v"; } print "; "; }
    if($self->{_compo_ex}) { print "Units Ex:"; while ( my ($k,$v) = each %{$self->getComponentsEx()} ) { print " $k=>$v"; } print "; "; }
    if($self->{_comm})     { print "Comm: "  .  $self->getComment() . "; ";  }
    if($self->{_func})     { 
        print "\n\tFunc\n ";    
        print "\t\txpr:    " . ${$self->getFunction()}{xpr}         . "\n";
        print "\t\tdir:    " . "@{${$self->getFunction()}{dir}}"    . "\n";
        print "\t\tcompo1: " . "@{${$self->getFunction()}{compo1}}" . "\n";
        print "\t\tsign1:  " . "@{${$self->getFunction()}{sign1}}"  . "\n";
        print "\t\tcompo2: " . "@{${$self->getFunction()}{compo2}}" . "\n";
        print "\t\tsign2:  " . "@{${$self->getFunction()}{sign2}}"  . "\n";
    }
    if($self->{_group})    { print "Group: " .  $self->getGroup() . "; ";    }
    if($self->{_quota})    { print "Quota: " .  $self->getQuota() . "; ";    }
    print "\n";
}


###############Class Group
package Group;

sub new{
   my $class = shift;
   my $self = {
       _sigla  => shift,
       _dim    => shift,
       _type    => shift,
       _compo  => shift,
       _include  => shift,
       _z        => shift,
   };
   return bless $self, $class;   
}

sub DESTROY{}

sub getSigla{ my( $self ) = @_; return $self->{_sigla}; }
sub getDim{ my( $self ) = @_; return $self->{_dim}; }
sub getType{ my( $self ) = @_; return $self->{_type}; }
sub getComponents{ my( $self ) = @_; return $self->{_compo}; }
sub getInclude{ my( $self ) = @_; return $self->{_include}; }
sub getZ{ my( $self ) = @_; return $self->{_z}; }

sub setSigla{ my ( $self, $sigla ) = @_; $self->{_sigla} = $sigla if defined($sigla);}
sub setDim{ my ( $self, $dim ) = @_; $self->{_dim} = $dim if defined($dim);}
sub setType{ my ( $self, $type ) = @_; $self->{_type} = $type if defined($type);}
sub setComponents{ my ( $self, $compo ) = @_; $self->{_compo} = $compo if defined($compo);}
sub setInclude{ my ( $self, $include ) = @_; $self->{_include} = $include if defined($include);}
sub setZ{ my ( $self, $z ) = @_; $self->{_z} = $z if defined($z);}

sub print{ 
    my ( $self ) = @_; 
    print "GROUP => ";
    if($self->{_sigla})  { print $self->{_sigla} . "\n";                   }
    if($self->{_dim})    { print "\tDim:   " . $self->getDim() . "\n";     }
    if($self->{_type})   { print "\tType:  " . $self->getType() . "\n";    }
    if($self->{_include}){ print "\tIncl:  " . $self->getInclude() . "\n"; }
    if($self->{_z})      { print "\t-Z:    " . $self->getZ() . "\n";       }
    if($self->{_compo})  { print "\tUnits:" ; while ( my ($k,$v) = each %{$self->getComponents()} ) { print " $k=>$v"; } print "\n"; }
    #print "\n";
}



1;
