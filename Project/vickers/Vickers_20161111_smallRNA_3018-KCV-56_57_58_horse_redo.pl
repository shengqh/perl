#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;

my $def = {

	#General options
	task_name                 => "KCV-56_57_58_horse",
	email                     => "quanhu.sheng\@vanderbilt.edu",
	target_dir                => "/scratch/cqs/shengq1/vickers/20161111_smallRNA_3018-KCV-56_57_58_horse_redo",
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
        "HorseHDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-56-i4_S4_R1_001.fastq.gz'],
        "HorseLDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-57-i20_S4_R1_001.fastq.gz'],
        "HorseVLDL"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-58-i35_S5_R1_001.fastq.gz'],

	},
};

performSmallRNA_eca2($def);

1;

