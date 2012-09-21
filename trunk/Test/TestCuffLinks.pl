#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::StringUtils;

my @samples = ( "1", "3", "4", "5", "10", "11", "13", "16" );

my $cufflinkparam = "-p 6";
my $rootDir = "/scratch/cqs/shengq1/rnaseq/1769";

my $runNow = 0;
if ($#ARGV > 0){
	my $isRunNow = $ARGV[0]; 
	$runNow = $isRunNow eq "y";
}

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
