#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {

	#General options
	task_name  => "2570",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201501_smallRNA_2570_rat/",
	max_thread => 8,
	cluster    => "slurm",

	#Data
	files => {
		"2570-KCV-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-1_1_sequence.txt.gz"],
		"2570-KCV-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-2_1_sequence.txt.gz"],
		"2570-KCV-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-3_1_sequence.txt.gz"],
		"2570-KCV-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-4_1_sequence.txt.gz"],
		"2570-KCV-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-5_1_sequence.txt.gz"],
		"2570-KCV-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-6_1_sequence.txt.gz"],
		"2570-KCV-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-7_1_sequence.txt.gz"],
		"2570-KCV-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-8_1_sequence.txt.gz"],
		"2570-KCV-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-9_1_sequence.txt.gz"],
		"2570-KCV-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-10_1_sequence.txt.gz"],
		"2570-KCV-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-11_1_sequence.txt.gz"],
		"2570-KCV-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-12_1_sequence.txt.gz"],
		"2570-KCV-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-13_1_sequence.txt.gz"],
		"2570-KCV-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-14_1_sequence.txt.gz"],
		"2570-KCV-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-15_1_sequence.txt.gz"],
		"2570-KCV-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-16_1_sequence.txt.gz"],
		"2570-KCV-17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-17_1_sequence.txt.gz"],
		"2570-KCV-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2570_Rat_Liver_diabetes/2570-KCV-18_1_sequence.txt.gz"],
	},
};

performSmallRNARat($def);

1;

