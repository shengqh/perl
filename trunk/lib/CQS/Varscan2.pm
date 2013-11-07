#!/usr/bin/perl
package CQS::Varscan2;

use strict;
use warnings;

sub new {
  require "Varscan2/Somatic.pm";
  my $class = "VarScan2::Somatic";
  return $class->new();
}

1;
