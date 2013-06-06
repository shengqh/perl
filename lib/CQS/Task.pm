#!/usr/bin/perl
package CQS::Task;

use strict;
use warnings;

sub new {
  my ($class) = @_;
  my $self = { _name => undef };
  bless $self, $class;
  return $self;
}

sub name {
  my ($self) = @_;
  return $self->{_name};
}

sub perform {
}

sub filter {
  my ($self, $resultFiles, $pattern ) = @_;

  if ( !defined $pattern ) {
    return $resultFiles;
  }

  #print $resultFiles . "\n";
  my @filteredFiles = ();

  for my $candidateFile ( @{$resultFiles} ) {
    if ( $candidateFile =~ m/$pattern/ ) {
      push( @filteredFiles, $candidateFile );
    }
  }

  return \@filteredFiles;
}

sub result {
  my $result = {};
  return $result;
}

1;
