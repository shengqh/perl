#!/usr/bin/perl
package CQS::ClassFactory;

use strict;
use warnings;
use CQS::Cutadapt;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(new_class)] );

our @classes = ( new CQS::Cutadapt() );

sub new_class {
  my ($className) = @_;
  foreach my $class (@classes) {
    if ( $class->name() eq $className ) {
      return $class->clone();
    }
  }

  die "Cannot find class $className";
}

1;
