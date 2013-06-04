#!/usr/bin/perl
package CQS::NGSCommon;

use strict;
use warnings;

use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(get_sorted_bam get_sam2bam_command get_sort_index_command get_stat_command)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub get_sorted_bam {
  my $bamFile         = shift;
  my $bamSortedPrefix = change_extension( $bamFile, "_sorted" );
  my $bamSortedFile   = $bamSortedPrefix . ".bam";
  return ( $bamSortedFile, $bamSortedPrefix );
}

sub get_sam2bam_command {
  my ( $samFile, $bamFile, $indent ) = @_;

  if ( !defined($indent) ) {
    $indent = "";
  }

  my $command = "${indent}if [[ -s $samFile && ! -s $bamFile ]]; then
${indent}  echo sam2bam=`date`
${indent}  samtools view -b -S $samFile -o $bamFile
${indent}fi";

  return ($command);
}

sub get_sort_index_command {
  my ( $bamFile, $bamSortedPrefix, $indent ) = @_;

  if ( !defined($indent) ) {
    $indent = "";
  }

  my $bamSortedFile;
  if ( !defined $bamSortedPrefix ) {
    ( $bamSortedFile, $bamSortedPrefix ) = get_sorted_bam($bamFile);
  }
  else {
    $bamSortedFile = $bamSortedPrefix . ".bam";
  }

  my $command = "${indent}if [[ -s $bamFile && ! -s $bamSortedFile ]]; then
${indent}  echo BamSort=`date` 
${indent}  samtools sort $bamFile $bamSortedPrefix 
${indent}fi

${indent}if [[ -s $bamSortedFile && ! -s ${bamSortedFile}.bai ]]; then
${indent}  echo BamIndex=`date` 
${indent}  samtools index $bamSortedFile
${indent}fi";
  return ($command);
}

sub get_stat_command {
  my ( $bamSortedFile, $indent ) = @_;

  if ( !defined($indent) ) {
    $indent = "";
  }

  my $command = "${indent}if [[ -s $bamSortedFile && ! -s ${bamSortedFile}.stat ]]; then
${indent}  echo bamstat=`date`
${indent}  samtools flagstat $bamSortedFile > ${bamSortedFile}.stat 
${indent}fi";

  return ($command);
}

1;
