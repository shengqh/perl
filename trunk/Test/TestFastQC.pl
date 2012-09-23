#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my @samples = ( "1", "3", "4", "5", "10", "11", "13", "16" );

my $rootDir      = "/scratch/cqs/shengq1/rnaseq/1769";

my $runNow = get_run_now();

create_directory_or_die($rootDir);

my @sampleNames = ();
my @sampleFiles = ();

foreach my $sample (@samples) {
	my $name = "1769-DPC-" . $sample;

	my $fastqFile1 = "/scratch/cqs/guoy1/1769/rawdata/" . $name . "_1_sequence.txt";

	my $fastqFile2 = "/scratch/cqs/guoy1/1769/rawdata/" . $name . "_2_sequence.txt";

	push( @sampleNames, $name . "_1" );
	push( @sampleNames, $name . "_2" );
	push( @sampleFiles, $fastqFile1 );
	push( @sampleFiles, $fastqFile2 );
}

fastqc_by_pbs($rootDir, \@sampleNames, \@sampleFiles, $runNow);
