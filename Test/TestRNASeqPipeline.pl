#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;

my @samples = ("1","3","4","5","10","11","13","16");

my $dbFile = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";
my $rootDir = "/scratch/cqs/shengq1/rnaseq/1769";

unless(-e $rootDir or mkdir($rootDir)){
  die "Unable to create $rootDir\n";
}

foreach my $sample (@samples) {
  my $name = "1769-DPC-" . $sample;	

  my $fastqFile1 = "/scratch/cqs/guoy1/1769/rawdata/" . $name . "_1_sequence.txt";
	
  my $fastqFile2 = "/scratch/cqs/guoy1/1769/rawdata/" . $name . "_2_sequence.txt";
	
  my $curDir = $rootDir . "/" . $name;
	
  unless(-e $curDir or mkdir($curDir)){
    die "Unable to create $curDir\n";
  }

  fastqc_by_pbs($curDir, $name . "_1", $fastqFile1);

  fastqc_by_pbs($curDir, $name . "_2", $fastqFile2);

  #tophat2_by_pbs_double( $dbFile, $curDir, $name, $fastqFile1, $fastqFile2 );
}
