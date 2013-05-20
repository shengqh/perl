#!/usr/bin/perl
package CQS::NGSCommon;

use strict;
use warnings;

use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(get_sorted_bam)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub get_sorted_bam {
  my $bamFile         = shift;
  my $bamSortedPrefix = change_extension( $bamFile, "_sorted" );
  my $bamSortedFile   = $bamSortedPrefix . ".bam";
  return ( $bamSortedFile, $bamSortedPrefix );
}

1;
