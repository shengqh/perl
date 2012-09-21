#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;

my @samples = ( "1", "3", "4", "5", "10", "11", "13", "16" );

my $genomeDb     = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";
my $rootDir      = "/scratch/cqs/shengq1/rnaseq/1769_2";
my $gtfFile      = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $gtfIndex     = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";
my $tophat2param = "--segment-length 25 -r 0 -p 6";

my $runNow = 0;

unless ( -e $rootDir or mkdir($rootDir) ) {
	die "Unable to create $rootDir\n";
}

my @sampleNames = ();
my @sampleFiles = ();

foreach my $sample (@samples) {
	my $name = "1769-DPC-" . $sample;

	my $fastqFile1 = "/scratch/cqs/guoy1/1769/rawdata/" . $name . "_1_sequence.txt";

	my $fastqFile2 = "/scratch/cqs/guoy1/1769/rawdata/" . $name . "_2_sequence.txt";

	push( @sampleNames, $name );
	push( @sampleFiles, $fastqFile1 );
	push( @sampleFiles, $fastqFile2 );
}

#tophat2_by_pbs_batch( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $rootDir, "test1769", \@sampleNames, \@sampleFiles );
tophat2_by_pbs_individual( $genomeDb, $gtfFile, $gtfIndex, $tophat2param, $rootDir, \@sampleNames, \@sampleFiles, $runNow );
#tophat2_by_pbs_individual( $genomeDb, "", "", $tophat2param, $rootDir, \@sampleNames, \@sampleFiles );
