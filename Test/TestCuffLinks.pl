#!/usr/bin/perl -w
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;

my $runNow = get_run_now();

my $root_dir = "/scratch/cqs/shengq1/rnaseq/1769";

my $config = {
	general => {
		root_dir  => $root_dir,
		path_file => "/home/shengq1/bin/path.txt"
	},
	cufflinks => { option => "-p 8" },
	pbs       => {
		"email"    => "quanhu.sheng\@vanderbilt.edu",
		"nodes"    => "8",
		"walltime" => "72",
		"mem"      => "20000mb"
	},
	tophat2_result => {
		"1769-DPC-1" =>
		  "${root_dir}/result/tophat2/1769-DPC-1/accepted_hits.bam",
		"1769-DPC-3" =>
		  "${root_dir}/result/tophat2/1769-DPC-3/accepted_hits.bam",
		"1769-DPC-4" =>
		  "${root_dir}/result/tophat2/1769-DPC-4/accepted_hits.bam",
		"1769-DPC-5" =>
		  "${root_dir}/result/tophat2/1769-DPC-5/accepted_hits.bam",
		"1769-DPC-10" =>
		  "${root_dir}/result/tophat2/1769-DPC-10/accepted_hits.bam",
		"1769-DPC-11" =>
		  "${root_dir}/result/tophat2/1769-DPC-11/accepted_hits.bam",
		"1769-DPC-13" =>
		  "${root_dir}/result/tophat2/1769-DPC-13/accepted_hits.bam",
		"1769-DPC-16" =>
		  "${root_dir}/result/tophat2/1769-DPC-16/accepted_hits.bam"
	}
};

cufflinks_by_pbs( $config, $runNow );
