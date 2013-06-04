#!/usr/bin/perl
package CQS::ClassFactory;

use strict;
use warnings;
use CQS::Cutadapt;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(new_class)] );

sub new_class {
  my ($className) = @_;
  my @classes = ( new CQS::Cutadapt() );
  foreach my $class (@classes) {
    if ( $class->name() eq $className ) {
      return $class->clone();
    }
  }

  die "Cannot find class $className";
}

1;
