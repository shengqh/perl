#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {

	#General options
	task_name  => "2795",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/",
	max_thread => 8,
	cluster    => "slurm",

	#Data
	files => {
		"2795-KCV-1_RPI12" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI12.fastq.gz"],
		"2795-KCV-1_RPI16" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI16.fastq.gz"],
		"2795-KCV-1_RPI20" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI20.fastq.gz"],
		"2795-KCV-1_RPI24" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI24.fastq.gz"],
		"2795-KCV-1_RPI28" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI28.fastq.gz"],
		"2795-KCV-1_RPI34" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI34.fastq.gz"],
		"2795-KCV-1_RPI38" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI38.fastq.gz"],
		"2795-KCV-1_RPI43" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-1_RPI43.fastq.gz"],
		"2795-KCV-2_RPI11" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI11.fastq.gz"],
		"2795-KCV-2_RPI15" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI15.fastq.gz"],
		"2795-KCV-2_RPI19" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI19.fastq.gz"],
		"2795-KCV-2_RPI23" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI23.fastq.gz"],
		"2795-KCV-2_RPI27" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI27.fastq.gz"],
		"2795-KCV-2_RPI31" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI31.fastq.gz"],
		"2795-KCV-2_RPI36" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI36.fastq.gz"],
		"2795-KCV-2_RPI41" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-2_RPI41.fastq.gz"],
		"2795-KCV-3_RPI13" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI13.fastq.gz"],
		"2795-KCV-3_RPI17" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI17.fastq.gz"],
		"2795-KCV-3_RPI21" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI21.fastq.gz"],
		"2795-KCV-3_RPI25" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI25.fastq.gz"],
		"2795-KCV-3_RPI29" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI29.fastq.gz"],
		"2795-KCV-3_RPI35" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI35.fastq.gz"],
		"2795-KCV-3_RPI39" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI39.fastq.gz"],
		"2795-KCV-3_RPI42" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-3_RPI42.fastq.gz"],
		"2795-KCV-4_RPI14" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-4_RPI14.fastq.gz"],
		"2795-KCV-4_RPI18" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-4_RPI18.fastq.gz"],
		"2795-KCV-4_RPI22" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-4_RPI22.fastq.gz"],
		"2795-KCV-4_RPI26" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-4_RPI26.fastq.gz"],
		"2795-KCV-4_RPI30" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-4_RPI30.fastq.gz"],
		"2795-KCV-4_RPI37" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-4_RPI37.fastq.gz"],
		"2795-KCV-4_RPI40" => ["/gpfs21/scratch/cqs/shengq1/vickers/201501_smallRNA_2795_rat/demultiplexing/result/2795-KCV-4_RPI40.fastq.gz"],
	},
};

performSmallRNARat($def);

1;

