#!/usr/bin/perl
package CQS::ClassFactory;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(instantiate)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub instantiate {
  my $requested_type = shift;
  my $location       = "CQS/$requested_type.pm";
  my $class          = "CQS::$requested_type";

  require $location;

  return $class->new(@_);
}

1;
