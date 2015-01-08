#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {

	#General options
	task_name  => "2516",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201501_smallRNA_2516_rat/",
	max_thread => 8,
	cluster    => "slurm",

	#Data
	files => {
		"2516-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-01_clipped.fastq.gz"],
		"2516-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-02_clipped.fastq.gz"],
		"2516-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-03_clipped.fastq.gz"],
		"2516-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-04_clipped.fastq.gz"],
		"2516-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-05_clipped.fastq.gz"],
		"2516-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-06_clipped.fastq.gz"],
		"2516-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-07_clipped.fastq.gz"],
		"2516-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-08_clipped.fastq.gz"],
		"2516-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2516_Rat_Islets_diabetes/2516-09_clipped.fastq.gz"],
	},
};

performSmallRNARat($def);

1;

