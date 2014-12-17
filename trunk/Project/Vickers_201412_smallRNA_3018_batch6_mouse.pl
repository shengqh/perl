#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {

	#General options
	task_name  => "3018_b6",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/",
	max_thread => 8,

	#Data
	files => {
		"3018-KCV-6-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-8_ACTTGA_L005_R1_001.fastq.gz"],
		"3018-KCV-6-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-9_GATCAG_L005_R1_001.fastq.gz"],
		"3018-KCV-6-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-10_TAGCTT_L005_R1_001.fastq.gz"],
		"3018-KCV-6-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-11_GGCTAC_L005_R1_001.fastq.gz"],
		"3018-KCV-6-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-12_CTTGTA_L005_R1_001.fastq.gz"],
		"3018-KCV-6-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-13_AGTCAA_L005_R1_001.fastq.gz"],
	},
};

performSmallRNAMouse($def);

1;

