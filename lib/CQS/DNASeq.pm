#!/usr/bin/perl
package CQS::DNASeq;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(bwa_by_pbs_single bwa_by_pbs_double)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub bwa_by_pbs_single {
  my ( $faFile, $sampleFile, $rootDir ) = @_;

  my ( $sampleName, $directories, $suffix ) = fileparse($sampleFile);

  $sampleName =~ s/\.(\w+)$//;

  my $pathFile = '/data/cqs/bin/path.txt';

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

  #`qsub $pbsFile`;
}

sub bwa_by_pbs_double {
  my ( $faFile, $sampleFile1, $sampleFile2, $sampleName, $rootDir ) = @_;

  my $pathFile = '/data/cqs/bin/path.txt';

  my ( $logDir, $pbsDir, $resultDir ) = init_dir($rootDir);

  my ($pbsDesc) = get_pbs_desc();

  my ( $sampleName1, $directories1, $suffix1 ) = fileparse($sampleFile1);
  my $saiFile1 = $sampleName1 . ".sai";
  my ( $sampleName2, $directories2, $suffix2 ) = fileparse($sampleFile2);
  my $saiFile2      = $sampleName2 . ".sai";
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
      if ( !( -e $resultDir . "/" . $saiFile1 ) ) {
        print OUT "echo sai=`date` \n";
        print OUT "bwa aln -q 15 $faFile $sampleFile1 >$saiFile1 \n";
      }
      if ( !( -e $resultDir . "/" . $saiFile2 ) ) {
        print OUT "echo sai=`date` \n";
        print OUT "bwa aln -q 15 $faFile $sampleFile2 >$saiFile2 \n";
      }
      print OUT "echo aln=`date` \n";
      print OUT
"bwa sampe -n 3 $faFile $saiFile1 $saiFile2 $sampleFile1 $sampleFile2 > $samFile \n";
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

  #`qsub $pbsFile`;
}

1;
