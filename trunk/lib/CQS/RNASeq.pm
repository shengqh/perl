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
  my ( $dbDir, $rootDir, $sampleName, $sampleFile ) = @_;
  tophat2_by_pbs($dbDir, $rootDir, $sampleName, $sampleFile);
}

sub tophat2_by_pbs_double {
  my ( $dbDir, $rootDir, $sampleName, $sampleFile1, $sampleFile2 ) = @_;
  tophat2_by_pbs($dbDir, $rootDir, $sampleName, $sampleFile1, $sampleFile2);
}

sub tophat2_by_pbs {
  my ( $dbDir, $rootDir, $sampleName, @sampleFiles ) = @_;

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
  print OUT "echo tophat2=`date` \n";
  print OUT "tophat2 --segment-length 25 -r 0 -p 8 -o $tophatDir $dbDir ";
  foreach my $sampleFile in (@sampleFiles){
    print OUT " $sampleFile";
  }
  print OUT "\n";
  print OUT "echo finished=`date`\n";
  close OUT;

  `qsub $pbsFile`;
}

1;
