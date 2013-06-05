#!/usr/bin/perl
package CQS::CQSTools;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(mirna_count)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub mirna_count {
  my ( $config, $section ) = @_;
  my $obj = instantiate("MirnaCount");
  $obj->perform( $config, $section );
}

1;
