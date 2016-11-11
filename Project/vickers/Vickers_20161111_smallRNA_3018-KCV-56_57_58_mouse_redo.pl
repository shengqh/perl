#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;

my $def = {

	#General options
	task_name                 => "KCV-56_57_58_mouse",
	email                     => "quanhu.sheng\@vanderbilt.edu",
	target_dir                => "/scratch/cqs/shengq1/vickers/20161111_smallRNA_3018-KCV-56_57_58_mouse_redo",
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
        "MouseHDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-56-i13_S13_R1_001.fastq.gz'],
        "MouseLDL"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-57-i29_S13_R1_001.fastq.gz'],
        "MouseVLDL"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-56-58/3018-KCV-58-i44_S14_R1_001.fastq.gz'],

	},
};

performSmallRNA_mm10($def);


1;

