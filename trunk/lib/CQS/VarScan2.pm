#!/usr/bin/perl
package CQS::VarScan2;

use strict;
use warnings;

sub new {
  require "VarScan2/Somatic.pm";
  my $class = "VarScan2::Somatic";
  return $class->new();
}

1;
