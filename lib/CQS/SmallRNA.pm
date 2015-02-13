#!/usr/bin/perl
package CQS::SmallRNA;

use strict;
use warnings;
require CQS::PerformSmallRNA;
require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(performSmallRNA_hg19 performSmallRNATask_hg19 performSmallRNA_hg20 performSmallRNATask_hg20 performSmallRNA_mm10 performSmallRNATask_mm10 performSmallRNA_rn5 performSmallRNATask_rn5 performSmallRNA_cel235 performSmallRNATask_cel235)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub performSmallRNA_hg19 {
  CQS::PerformSmallRNA::performSmallRNA_hg19(@_);
}

sub performSmallRNATask_hg19 {
  CQS::PerformSmallRNA::performSmallRNATask_hg19(@_);
}

sub performSmallRNA_hg20 {
  CQS::PerformSmallRNA::performSmallRNA_hg20(@_);
}

sub performSmallRNATask_hg20 {
  CQS::PerformSmallRNA::performSmallRNATask_hg20(@_);
}

sub performSmallRNA_mm10 {
  CQS::PerformSmallRNA::performSmallRNA_mm10(@_);
}

sub performSmallRNATask_mm10 {
  CQS::PerformSmallRNA::performSmallRNATask_mm10(@_);
}

sub performSmallRNA_rn5 {
  CQS::PerformSmallRNA::performSmallRNA_rn5(@_);
}

sub performSmallRNATask_rn5 {
  CQS::PerformSmallRNA::performSmallRNATask_rn5(@_);
}

sub performSmallRNA_cel235 {
  CQS::PerformSmallRNA::performSmallRNA_cel235(@_);
}

sub performSmallRNATask_cel235 {
  CQS::PerformSmallRNA::performSmallRNATask_cel235(@_);
}

1;
