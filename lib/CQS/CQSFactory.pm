#!/usr/bin/perl
package CQS::CQSFactory;

use strict;
use warnings;
use CQS::Cutadapt;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(new_object)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

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
