#!/usr/bin/perl
package CQS::ObjectFactory;

use strict;
use warnings;
use CQS::Cutadapt;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(new_object)] );

sub new_object {
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
