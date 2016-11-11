#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;

my $def = {

	#General options
	task_name                 => "KCV-56_57_58_human",
	email                     => "quanhu.sheng\@vanderbilt.edu",
	target_dir                => "/scratch/cqs/shengq1/vickers/20161111_smallRNA_3018-KCV-56_57_58_human_redo",
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
        "HumanHDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-56-i15_S15_R1_001.fastq.gz'],
        "HumanLDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-57-i31_S15_R1_001.fastq.gz'],
        "HumanVLDL"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-58-i46_S16_R1_001.fastq.gz'],

	},
};

performSmallRNA_hg19($def);


1;

