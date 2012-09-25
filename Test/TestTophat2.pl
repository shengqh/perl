#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;

my $runNow = get_run_now();

my $root_dir = "/scratch/cqs/shengq1/rnaseq/1769_Test";

my $config = {
	general => {
		root_dir             => $root_dir,
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf       => "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf",
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
		paired_data          => 1,
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "1769-DPC"
	},
	tophat2 => {
		option    => "--segment-length 25 -r 0 -p 8",
		batchmode => 1
	},
	pbs => {
		"email"    => "quanhu.sheng\@vanderbilt.edu",
		"nodes"    => "8",
		"walltime" => "72",
		"mem"      => "20000mb"
	},
	fastqfiles => {
		"G1" => {
			"1769-DPC-1" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-1_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-1_2_sequence.txt" ],
			"1769-DPC-3" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-3_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-3_2_sequence.txt" ],
			"1769-DPC-4" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-4_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-4_2_sequence.txt" ],
			"1769-DPC-5" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-5_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-5_2_sequence.txt" ],
		},
		"G2" => {
			"1769-DPC-10" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-10_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-10_2_sequence.txt" ],
			"1769-DPC-11" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-11_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-11_2_sequence.txt" ],
			"1769-DPC-13" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-13_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-13_2_sequence.txt" ],
			"1769-DPC-16" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-16_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-16_2_sequence.txt" ],
		}
	}
};

tophat2_by_pbs( $config, $runNow );
