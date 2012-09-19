#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(tophat2_by_pbs_single tophat2_by_pbs_double)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub tophat2_by_pbs_single {
  my ( $faFile, $sampleFile, $rootDir ) = @_;

  my ( $sampleName, $directories, $suffix ) = fileparse($sampleFile);

  $sampleName =~ s/\.(\w+)$//;

  my $pathFile = '/home/shengq1/bin/path.txt';

  my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);

  my ($pbsDesc) = get_pbs_desc();

  my $saiFile       = $sampleName . ".sai";
  my $samFile       = $sampleName . ".sam";
  my $bamFile       = $sampleName . ".bam";
  my $sortedBamFile = $sampleName . "_sort";
  my $log           = $logDir . "/" . $sampleName . ".log";

  if ( -e $resultDir . "/" . $sortedBamFile ) {
    next;
  }

  #my $tag="'\@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:ILLUMINA'";
  my $pbsFile = $pbsDir . "/${sampleName}.pbs";
  print "$pbsDir\n";
  open( OUT, ">$pbsFile" ) or die $!;
  print OUT $pbsDesc;
  print OUT "#PBS -o $log\n";
  print OUT "#PBS -j oe\n\n";
  print OUT "source $pathFile\n";
  print OUT "cd $resultDir\n\n";

  if ( !( -e $resultDir . "/" . $bamFile ) ) {
    if ( !( -e $resultDir . "/" . $samFile ) ) {
      if ( !( -e $resultDir . "/" . $sampleFile ) ) {
        print OUT "echo sai=`date` \n";
        print OUT "bwa aln -q 15 $faFile $sampleFile >$saiFile \n";
      }
      print OUT "echo aln=`date` \n";
      print OUT "bwa sampe -n 3 $faFile $saiFile $sampleFile > $samFile \n";
    }
    print OUT "echo sam2bam=`date`\n";
    print OUT "samtools view -b -S $samFile -o $bamFile\n";
  }
  print OUT "echo sortbam=`date`\n";
  print OUT "samtools sort $bamFile $sortedBamFile\n";
  print OUT "echo bamstat=`date`\n";
  print OUT "samtools flagstat $sortedBamFile.bam > $sortedBamFile.bam.stat\n";
  print OUT "echo finished=`date`\n";
  close OUT;

  `qsub $pbsFile`;
}

sub tophat2_by_pbs_double {
  my ( $dbDir, $sampleFile1, $sampleFile2, $sampleName, $rootDir ) = @_;

  my $pathFile = '/home/shengq1/bin/path.txt';

  my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);

  my $tophatDir = $resultDir . "/tophat2";

  unless(-e $tophatDir or mkdir($tophatDir)){
    die "Cannot create directory $tophatDir\n";
  }

  my ($pbsDesc) = get_pbs_desc();

  my $log = $logDir . "/${sampleName}_tophat2.log";

  my $pbsFile = $pbsDir . "/${sampleName}_tophat2.pbs";

  print "$pbsDir\n";
  open( OUT, ">$pbsFile" ) or die $!;
  print OUT $pbsDesc;
  print OUT "#PBS -o $log\n";
  print OUT "#PBS -j oe\n\n";
  print OUT "source $pathFile\n";
  print OUT "cd $tophatDir\n\n";

  print OUT "echo tophat2=`date` \n";
  print OUT "tophat2 --segment-length 25 -r 0 -p 8 -o $tophatDir $dbDir $sampleFile1 $sampleFile2\n";
  print OUT "echo finished=`date`\n";
  close OUT;

  `qsub $pbsFile`;
}

1;
