#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/VANGARD_201411_Vickers_pirmiR_2518/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $samtools = "/home/shengq1/local/bin/samtools/samtools";

my $overlap_count_option = "--gtf_key pri-miR --min_overlap 0.5 -e 2 --not_smallrna";

my $config = {
	general => { "task_name" => "2518", },
	files   => {
		"2518_KCV_01" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/1_Gout/accepted_hits.bam"],
		"2518_KCV_02" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/2_Gout/accepted_hits.bam"],
		"2518_KCV_03" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/3_Gout/accepted_hits.bam"],
		"2518_KCV_04" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/4_Gout/accepted_hits.bam"],
		"2518_KCV_05" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/5_Gout/accepted_hits.bam"],
		"2518_KCV_06" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/6_Gout/accepted_hits.bam"],
		"2518_KCV_07" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/7_Gout/accepted_hits.bam"],
		"2518_KCV_08" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/8_Gout/accepted_hits.bam"],
		"2518_KCV_09" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/9_Gout/accepted_hits.bam"],
		"2518_KCV_10" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/10_Gout/accepted_hits.bam"],
		"2518_KCV_11" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/11_Gout/accepted_hits.bam"],
		"2518_KCV_12" => ["/gpfs21/scratch/cqs/guom1/2518/tophat_G/12_Gout/accepted_hits.bam"],
	},
	primiR_count => {
		class      => "CQSMappedCount",
		perform    => 1,
		target_dir => "${root}/count_overlap",
		option     => $overlap_count_option,
		source_ref => "files",
		cqs_tools  => $cqstools,
		gff_file   => "/scratch/cqs/shengq1/vangard/VANGARD_Vickers/database/20141119_hg19_primiR.gff3",
		samtools   => $samtools,
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "72",
			"mem"      => "20gb"
		},
	},
	primiR_table => {
		class      => "CQSMappedTable",
		perform    => 1,
		target_dir => "${root}/primiR_table",
		option     => "",
		source_ref => [ "primiR_count", ".xml" ],
		cqs_tools  => $cqstools,
		prefix     => "primiR_",
		sh_direct  => 1,
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "10",
			"mem"      => "10gb"
		},
	},

};

performConfig($config);

1;

