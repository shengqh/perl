#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def_HDL = {

	#General options
	task_name  => "3018_b8_HDL",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/HDL/",
	max_thread => 8,

	#Data
	files => {
		"3018-KCV-8_KetoHDL1"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL1.fastq.gz"],
		"3018-KCV-8_KetoHDL2"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL2.fastq.gz"],
		"3018-KCV-8_KetoHDL3"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL3.fastq.gz"],
		"3018-KCV-8_KetoHDL4"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL4.fastq.gz"],
		"3018-KCV-8_KetoHDL5"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL5.fastq.gz"],
		"3018-KCV-8_KetoHDL6"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL6.fastq.gz"],
		"3018-KCV-8_KetoHDL7"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL7.fastq.gz"],
		"3018-KCV-8_KetoHDL8"    => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoHDL8.fastq.gz"],
	}
};

performSmallRNAHuman($def_HDL);

my $def_Plasma = {

	#General options
	task_name  => "3018_b8_Plasma",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/Plasma/",
	max_thread => 8,

	#Data
	files => {
		"3018-KCV-8_KetoPlasma1" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma1.fastq.gz"],
		"3018-KCV-8_KetoPlasma2" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma2.fastq.gz"],
		"3018-KCV-8_KetoPlasma3" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma3.fastq.gz"],
		"3018-KCV-8_KetoPlasma4" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma4.fastq.gz"],
		"3018-KCV-8_KetoPlasma5" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma5.fastq.gz"],
		"3018-KCV-8_KetoPlasma6" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma6.fastq.gz"],
		"3018-KCV-8_KetoPlasma7" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma7.fastq.gz"],
		"3018-KCV-8_KetoPlasma8" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch8_human/raw/3018-KCV-8_KetoPlasma8.fastq.gz"],
	}
};

performSmallRNAHuman($def_Plasma);

1;

