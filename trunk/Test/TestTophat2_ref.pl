#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my @samples = ( "1", "3", "4", "5", "10", "11", "13", "16" );

my $tophat2ParamRef = {
   # "root_dir"     => "/scratch/cqs/shengq1/rnaseq/1769_2",
	"genome_db"    => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
	"gtf_file"     => "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf",
	"gtf_index"    => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
	"tophat2_param" => "--segment-length 25 -r 0 -p 8",
	"path_file"    => "/home/shengq1/bin/path.txt"
};

my $pbsParamRef = {
	"email"    => "quanhu.sheng\@vanderbilt.edu",
	"nodes"    => "8",
	"walltime" => "72",
	"mem"      => "15000mb"
};

my $runNow = get_run_now();

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

tophat2_by_pbs_individual2( $tophat2ParamRef, \@sampleNames, \@sampleFiles, $pbsParamRef, $runNow );
