#!/usr/bin/perl
package CQS::QC;

use strict;
use warnings;
use File::Basename;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(fastqc_by_pbs)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub fastqc_by_pbs {
  my ( $rootDir, $seqFile ) = @_;

  my $pathFile = '/home/shengq1/bin/path.txt';

  my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);

  my ($pbsDesc) = get_pbs_desc();

  my $fastqcDir = $resultDir . "/fastqc";

  unless(-e $fastqcDir or mkdir($fastqcDir)){
    die "Cannot create directory $fastqcDir\n";
  }

  my $log = $logDir . "/${sampleName}_fastqc.log";

  my $pbsFile = $pbsDir . "/${sampleName}_fastqc.pbs";

  print "$pbsDir\n";
  open( OUT, ">$pbsFile" ) or die $!;
  print OUT $pbsDesc;
  print OUT "#PBS -o $log\n";
  print OUT "#PBS -j oe\n\n";
  print OUT "source $pathFile\n";
  print OUT "echo fastqc=`date`\n";
  print OUT "fastqc -o $fastqcDir $seqFile\n";
  print OUT "echo finished=`date`\n";
  close OUT;

  `qsub $pbsFile`;
}

1;
