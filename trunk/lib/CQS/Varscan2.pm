#!/usr/bin/perl
package CQS::Varscan2;

use strict;
use warnings;

sub new {
  require "Varscan2/SomaticMutation.pm";
  my $class = "VarScan2::SomaticMutation";
  return $class->new();
}

1;
