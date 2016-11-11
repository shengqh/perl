#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;

my $def = {

	#General options
	task_name                 => "KCV-56_57_58_cat",
	email                     => "quanhu.sheng\@vanderbilt.edu",
	target_dir                => "/scratch/cqs/shengq1/vickers/20161111_smallRNA_3018-KCV-56_57_58_cat_redo",
	max_thread                => 8,
	cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
	sequencetask_run_time     => 4,

	#Default software parameter (don't change it except you really know it)
	fastq_remove_N        => 1,
	remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
	search_unmapped_reads => 1,
	blast_unmapped_reads  => 0,
    

	#Data
	files => {
		"CatHDL"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-56-i10_S10_R1_001.fastq.gz'],
		"CatLDL"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-57-i26_S10_R1_001.fastq.gz'],
		"Cat1VLDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-58-i41_S11_R1_001.fastq.gz'],
		"Cat2VLDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-58-i42_S12_R1_001.fastq.gz'],
	},
};

performSmallRNA_fca6($def);

1;

