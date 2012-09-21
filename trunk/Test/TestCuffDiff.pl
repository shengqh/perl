#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::StringUtils;

my @samples1 = ( "1", "3", "4", "5" );
my @samples2 = (  "10", "11", "13", "16" );

my $genomeFasta     = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa";
my $gtfFile      = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
#my $gffIndex     = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68.gff";
my $cuffdiffparam = "-p 6 -N";
my $rootDir      = "/scratch/cqs/shengq1/rnaseq/1769";

create_directory_or_die($rootDir);

my @bamFiles1 = ();
my @bamFiles2 = ();

foreach my $sample (@samples1) {
	my $name = "1769-DPC-" . $sample;
	my $bamfile = $rootDir . "/" . $name . "/result/tophat2/accepted_hits.bam";
	push( @bamFiles1, $bamfile );
}

foreach my $sample (@samples2) {
    my $name = "1769-DPC-" . $sample;
    my $bamfile = $rootDir . "/" . $name . "/result/tophat2/accepted_hits.bam";
    push( @bamFiles2, $bamfile );
}

my @files = (merge_string(",", @bamFiles1), merge_string(",", @bamFiles2));

cuffdiff_by_pbs($genomeFasta, $gtfFile, $cuffdiffparam, $rootDir, "test", "G1,G2", @files);
