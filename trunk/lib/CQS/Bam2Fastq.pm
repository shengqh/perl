#!/usr/bin/perl
package CQS::Bam2Fastq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "Bam2Fastq";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $ispaired = $config->{$section}{ispaired};
  if ( !defined $ispaired ) {
    $ispaired = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}_b2q.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $bamfile     = $sampleFiles[0];

    my $pbsName = "${sampleName}_b2q.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_b2q.log";

    open( OUT, ">$pbsFile" ) or die $!;

    if ($ispaired) {
      my $namesorted = $sampleName . ".namesorted";
      my $fastq     = $sampleName . ".fastq";
      my $fastq1     = $sampleName . ".1.fastq";
      my $fastq2     = $sampleName . ".2.fastq";
      my $finalFile1 = $fastq1 . ".gz";
      my $finalFile2 = $fastq2 . ".gz";
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

echo started=`date`

if [ ! -s $finalFile1 ]; then
  samtools sort -n $bamfile -o $namesorted
  bam2fastx $option -P -o $fastq ${namesorted}.bam
  gzip $fastq1
  gzip $fastq2
fi

echo finished=`date`

exit 1 
";
    }
    else {
      my $fastq     = $sampleName . ".fastq";
      my $finalFile = $fastq . ".gz";
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

cd $resultDir

echo started=`date`

if [ ! -s $finalFile ]; then
  bam2fastx $option -o $fastq $bamfile
  gzip $fastq
fi

echo finished=`date`

exit 1 
";
    }

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Bam2Fastq tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $ispaired = $config->{$section}{ispaired};
  if ( !defined $ispaired ) {
    $ispaired = 0;
  }

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    
    my @resultFiles = ();
    if ($ispaired) {
      my $namesorted = $sampleName . ".namesorted";
      my $fastq     = $sampleName . ".fastq";
      my $fastq1     = $sampleName . ".1.fastq";
      my $fastq2     = $sampleName . ".2.fastq";
      my $finalFile1 = $fastq1 . ".gz";
      my $finalFile2 = $fastq2 . ".gz";
      push( @resultFiles, $finalFile1 );
      push( @resultFiles, $finalFile2 );
    }else{
      my $fastq     = $resultDir . "/" . $sampleName . ".fastq";
      my $finalFile = $fastq . ".gz";
      push( @resultFiles, $finalFile );
    }    

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
