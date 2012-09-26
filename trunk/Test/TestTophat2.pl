#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::SystemUtils;
use XML::Simple;
use Data::Dumper;

my $runNow = get_run_now();

my $raw_dir = "/scratch/cqs/guoy1/1769/rawdata/";

our $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf       => "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf",
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "1769-DPC"
	},
	tophat2 => {
		target_dir  => "/scratch/cqs/shengq1/rnaseq/1769_test/tophat2",
		option      => "--segment-length 25 -r 0 -p 8",
		batchmode   => 0,
		paired_data => 1,
		pbs         => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
		fastqfiles => {
			"G1" => {
				"1769-DPC-1" => [ $raw_dir . "1769-DPC-1_1_sequence.txt", $raw_dir . "1769-DPC-1_2_sequence.txt" ],
				"1769-DPC-3" => [ $raw_dir . "1769-DPC-3_1_sequence.txt", $raw_dir . "1769-DPC-3_2_sequence.txt" ],
				"1769-DPC-4" => [ $raw_dir . "1769-DPC-4_1_sequence.txt", $raw_dir . "1769-DPC-4_2_sequence.txt" ],
				"1769-DPC-5" => [ $raw_dir . "1769-DPC-5_1_sequence.txt", $raw_dir . "1769-DPC-5_2_sequence.txt" ],
			},
			"G2" => {
				"1769-DPC-10" => [ $raw_dir . "1769-DPC-10_1_sequence.txt", $raw_dir . "1769-DPC-10_2_sequence.txt" ],
				"1769-DPC-11" => [ $raw_dir . "1769-DPC-11_1_sequence.txt", $raw_dir . "1769-DPC-11_2_sequence.txt" ],
				"1769-DPC-13" => [ $raw_dir . "1769-DPC-13_1_sequence.txt", $raw_dir . "1769-DPC-13_2_sequence.txt" ],
				"1769-DPC-16" => [ $raw_dir . "1769-DPC-16_1_sequence.txt", $raw_dir . "1769-DPC-16_2_sequence.txt" ],
			}
		}
	},
};

#print Dumper $config;

my $xml = XMLout($config);

print $xml;

#tophat2_by_pbs( $config, $runNow );

1;
