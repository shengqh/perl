#!/usr/bin/perl
package CQS::ClassFactory;

use strict;
use warnings;

require CQS::Cutadapt;
require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(new_class)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

my @classes = ( new CQS::Cutadapt() );

sub new_class {
  my ($className) = @_;
  foreach my $class (@classes) {
    if ( $class->name() eq $className ) {
      return $class;
    }
  }

  die "Cannot find class $className";
}

1;
