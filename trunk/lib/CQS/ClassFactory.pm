#!/usr/bin/perl
package CQS::ClassFactory;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(instantiate performTask performConfig)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub instantiate {
  my $requested_type = shift;
  my $location       = "CQS/$requested_type.pm";
  my $class          = "CQS::$requested_type";

  require $location;

  return $class->new(@_);
}

sub performTask {
  my ( $config, $section ) = @_;
  my $classname = $config->{$section}{class} or die "No class defined in section $section.";
  my $obj = instantiate($classname);
  $obj->perform( $config, $section );
}

sub performConfig {
  my ( $config, $pattern, $force ) = @_;
  foreach my $section ( sort keys %{$config} ) {
    if ( !defined $pattern || $section =~ m/$pattern/ ) {
      my $classname = $config->{$section}{class};
      if ( defined $classname ) {
        my $perform;
        if ( defined $force && $force ) {
          $perform = 1;
        }
        else {
          $perform = $config->{$section}{perform};
        }
        
        if ( defined $perform && $perform ) {
          print "Preforming $section by $classname \n";
          my $obj = instantiate($classname);
          $obj->perform( $config, $section );
        }
      }
    }
  }
}

1;
