#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;

sub new {
  my ($class) = @_;
  my $self = { _name => undef };
  bless $self, $class;
  return $self;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub generateScript {
}

sub getExpectResult {
  my $result = {};
  return $result;
}

1;
