#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;

my @samples = ("1","3","4","5","10","11","13","16");

my $dbFile = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";
my $rootDir = "/scratch/cqs/shengq1/rnaseq";

foreach my $sample (@samples) {
	my $fastqFile1 = "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-" . $sample . "_1_sequence.txt";
	my $fastqFile2 = "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-" . $sample . "_2_sequence.txt";

	tophat2_by_pbs_double( $dbFile, $fastqFile1, $fastqFile2, "1769-DPC-" . $sample, "$rootDir/tophat2_double" );
}
