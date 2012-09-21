#!/usr/bin/perl -w
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::SystemUtils;

my @samples = ( "1", "3", "4", "5", "10", "11", "13", "16" );

my $cufflinkparam = "-p 8";
my $rootDir = "/scratch/cqs/shengq1/rnaseq/1769";

my $runNow = get_run_now();

create_directory_or_die($rootDir);

my @sampleNames = ();
my @sampleFiles = ();
foreach my $sample (@samples) {
	my $name    = "1769-DPC-" . $sample;
	my $bamfile = $rootDir . "/" . $name . "/result/tophat2/accepted_hits.bam";
	push( @sampleNames, $name );
	push( @sampleFiles, $bamfile );
}

cufflinks_by_pbs( $cufflinkparam, $rootDir, \@sampleNames, \@sampleFiles, $runNow );
