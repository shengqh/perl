#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;

my $runNow = get_run_now();

my $root_dir = "/scratch/cqs/shengq1/rnaseq/1769_2";

my $config = {
	general => {
		root_dir      => $root_dir,
		bowtie2_index => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf =>
"/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf",
		path_file => "/home/shengq1/bin/path.txt",
		task_name => "1769-DPC"
	},
	cuffmerge => {
		option    => "-p 8",
		assembies => $root_dir . "/result/cuffmerge/assemblies.txt"
	},
	pbs => {
		"email"    => "quanhu.sheng\@vanderbilt.edu",
		"nodes"    => "8",
		"walltime" => "72",
		"mem"      => "20000mb"
	}
};

cuffdiff_by_pbs( $config, $runNow );
