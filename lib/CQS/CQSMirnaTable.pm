#!/usr/bin/perl
package CQS::CQSMirnaTable;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "CQSMirnaTable";
  bless $self, $class;
  return $self;
}

sub get_result {
  my ( $task_name, $option ) = @_;

  my $result;
  if ( $option =~ /-o\s+(\S+)/ ) {
    $result = $1;
  }
  else {
    $result = $task_name . ".count";
    $option = $option . " -o " . $result;
  }
  return ( $result, $option );
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_mt.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "
cd $resultDir
";

  if ( defined $config->{$section}{groups} || defined $config->{$section}{groups_ref} ) {
    my $groups = get_raw_files( $config, $section, "groups" );
    for my $groupName ( sort keys %{$groups} ) {
      my @samples  = @{ $groups->{$groupName} };
      my $filelist = $pbsDir . "/${groupName}_mt.filelist";
      open( FL, ">$filelist" ) or die "Cannot create $filelist";
      for my $sampleName ( sort @samples ) {
        my @countFiles = @{ $rawFiles{$sampleName} };
        my $countFile  = $countFiles[0];
        print FL $sampleName, "\t", $countFile, "\n";
      }
      close(FL);
      my ( $result, $newoption ) = get_result( $groupName, $option );
      print SH "
mono-sgen $cqsFile mirna_table $newoption -l $filelist
";
    }
  }
  else {
    my $filelist = $pbsDir . "/${task_name}_mt.filelist";
    open( FL, ">$filelist" ) or die "Cannot create $filelist";
    for my $sampleName ( sort keys %rawFiles ) {
      my @countFiles = @{ $rawFiles{$sampleName} };
      my $countFile  = $countFiles[0];
      print FL $sampleName, "\t", $countFile, "\n";
    }
    close(FL);
    my ( $result, $newoption ) = get_result( $task_name, $option );

    print SH "
cd $resultDir

mono-sgen $cqsFile mirna_table $newoption -l $filelist
";
  }
  close SH;

  print "!!!shell file $shfile created, you can run this shell file to run CQSMirnaTable task.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = {};
  my ( $resultFile, $newoption ) = get_result( $task_name, $option );
  $resultFile = $resultDir . "/" . $resultFile;
  my $filelist = $resultDir . "/${task_name}.filelist";

  my @resultFiles = ();
  push( @resultFiles, $resultFile );
  push( @resultFiles, $filelist );

  $result->{$task_name} = filter_array( \@resultFiles, $pattern );

  return $result;
}

1;
