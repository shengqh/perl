#!/usr/bin/perl
package CQS::ConfigUtils;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;
use CQS::PBS;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(get_parameter get_param_file get_raw_files get_raw_files2)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_parameter {
  my ( $config, $section ) = @_;

  my $task_name = $config->{general}{task_name} or die "define general::task_name first";
  my $path_file = get_param_file( $config->{general}{path_file}, "path_file", 0 );
  if ( -e $path_file ) {
    $path_file = "source $path_file";
  }
  else {
    $path_file = "";
  }

  my $refPbs     = $config->{$section}{pbs}        or die "define ${section}::pbs parameters first";
  my $target_dir = $config->{$section}{target_dir} or die "define ${section}::target_dir parameters first";
  my ( $logDir, $pbsDir, $resultDir ) = init_dir($target_dir);
  my ($pbsDesc) = get_pbs_desc($refPbs);

  die "define ${section}::option first" if ( !defined $config->{$section}{option} );
  my $option    = $config->{$section}{option};
  my $sh_direct = $config->{$section}{sh_direct};
  if ( !defined $sh_direct ) {
    $sh_direct = 0;
  }

  return ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct );
}

#get parameter which indicates a file. If required, not defined or not exists, die. If defined but not exists, die.
#returned file either undefined or exists.
sub get_param_file {
  my ( $file, $name, $required ) = @_;

  my $result = $file;

  if ($required) {
    if ( !defined $file ) {
      die "$name was not defined!";
    }

    if ( !-e $file ) {
      die "$name $file defined but not exists!";
    }
  }
  else {
    if ( defined($file) ) {
      if ( $file eq "" ) {
        undef($result);
      }
      elsif ( !-e $file ) {
        die "$name $file defined but not exists!";
      }
    }
  }
  return ($result);
}

require CQS::ClassFactory;

sub do_get_raw_files {
  my ( $config, $section, $returnself ) = @_;
  if ( !defined $config->{$section} ) {
    die "section $section was not defined!";
  }

  if ( defined $config->{$section}{unmapped_ref} ) {
    my $alignsection = $config->{$section}{unmapped_ref};
    my $align_dir = $config->{$alignsection}{target_dir} or die "${$alignsection}::target_dir not defined.";
    my ( $logDir, $pbsDir, $resultDir ) = init_dir( $align_dir, 0 );
    my %fqFiles = %{ do_get_raw_files( $config, $alignsection, 0 ) };
    my $result = {};
    for my $sampleName ( keys %fqFiles ) {
      my $fq = "${resultDir}/${sampleName}/${sampleName}_sorted.bam.unmapped.fastq";
      $result->{"${sampleName}_unmapped"} = $fq;
    }
    return ( $result, 0 );
  }

  if ( defined $config->{$section}{source} ) {
    return ( $config->{$section}{source}, 1 );
  }

  if ( defined $config->{$section}{source_ref} ) {
    my $sectionName = $config->{$section}{source_ref};
    my ( $result, $issource ) = do_get_raw_files( $config, $sectionName, 1 );
    return ( $result, 0 );
  }

  if ($returnself) {
    return ( $config->{$section}, 0 );
  }
  else {
    die "define source or source_ref for $section";
  }
}

sub get_raw_files {
  my ( $config, $section ) = @_;
  my ( $result, $issource ) = do_get_raw_files( $config, $section, 0 );
  return $result;
}

#return raw files and if the raw files are extracted from source directly
sub get_raw_files2 {
  my ( $config, $section ) = @_;
  return do_get_raw_files( $config, $section, 0 );
}

1;

