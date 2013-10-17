#!/usr/bin/perl
package CQS::CQSSmallRNACategory;

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
  $self->{_name} = "CQSSmallRNACategory";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $cqsFile = get_param_file( $config->{$section}{cqs_tools}, "cqs_tools", 1 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $hasMiRNA = 0;
  my %miRNAFiles = ();
  if ( defined $config->{$section}{"mirna_count"} || defined $config->{$section}{"mirna_count_ref"} ) {
    %miRNAFiles = %{ get_raw_files( $config, $section, "mirna_count" ) };
    $hasMiRNA  =1;
  }
  
  my $ispdf = $config->{$section}{pdfgraph};
  if(!defined $ispdf){
    $ispdf = 0;
  }
  my $pdfoption = "";
  if($ispdf){
    $pdfoption = "--pdf";
  }

  my $shfile = $pbsDir . "/${task_name}_cat.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @bamFiles = @{ $rawFiles{$sampleName} };
    my $bamFile  = $bamFiles[0];
    my $smallrnaoption = "-s \"$bamFile\"";

    my $mirnaoption = "";
    if($hasMiRNA){
      my @mirnas      = @{ $miRNAFiles{$sampleName} };
      my $mirna       = $mirnas[0];
      $mirnaoption = "-m \"$mirna\"";
    }

    my $pbsName = "${sampleName}_cat.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";
    
    my $countFile = $sampleName . ".catcount";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_cat.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

if [ -s $countFile ]; then
  echo job has already been done. if you want to do again, delete $countFile and submit job again.
  exit 0
fi

echo CQSSmallRNACategory=`date` 

mono-sgen $cqsFile smallrna_category $option $mirnaoption $smallrnaoption -o $countFile $pdfoption

echo finished=`date`

exit 1 
";

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all trna_count tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my $countFile = "${resultDir}/${sampleName}.catcount";

    my @resultFiles = ();
    push( @resultFiles, $countFile );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
