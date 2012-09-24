#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;

my $runNow = get_run_now();

my $root_dir = "/scratch/cqs/shengq1/rnaseq/1769_2";

my $config = {
	general => {
		root_dir             => $root_dir,
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		#transcript_gtf       => "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf",
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "1769_2"
	},
	cuffdiff => { option => "-p 8 -N" },
	pbs      => {
		"email"    => "quanhu.sheng\@vanderbilt.edu",
		"nodes"    => "8",
		"walltime" => "72",
		"mem"      => "20000mb"
	},
	files => {
		G1 => [
			$root_dir . "/result/tophat2/1769-DPC-1/accepted_hits.bam",
			$root_dir . "/result/tophat2/1769-DPC-3/accepted_hits.bam",
			$root_dir . "/result/tophat2/1769-DPC-4/accepted_hits.bam",
			$root_dir . "/result/tophat2/1769-DPC-5/accepted_hits.bam"
		],
		G2 => [
			$root_dir . "/result/tophat2/1769-DPC-10/accepted_hits.bam",
			$root_dir . "/result/tophat2/1769-DPC-11/accepted_hits.bam",
			$root_dir . "/result/tophat2/1769-DPC-13/accepted_hits.bam",
			$root_dir . "/result/tophat2/1769-DPC-16/accepted_hits.bam"
		],
	}
};

cuffdiff_by_pbs( $config, $runNow );
