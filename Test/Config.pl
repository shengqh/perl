#!/usr/bin/perl
use strict;
use warnings;

package Task;

my $target_dir             = "/scratch/cqs/shengq1/rnaseq/1769_test";
my $tophat2_dir            = "${target_dir}/tophat2";
my $cufflinks_dir          = "${target_dir}/cufflinks";
my $cuffmerges_dir         = "${target_dir}/cuffmerge";
my $cuffdiff_dir           = "${target_dir}/cuffdiff";
my $cufflinks_cuffdiff_dir = "${target_dir}/cufflinks_cuffdiff";

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";

our $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "1769-DPC"
	},
	fastqfiles => {
		"G1" => {
			"1769-DPC-1" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-1_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-1_2_sequence.txt" ],
			"1769-DPC-3" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-3_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-3_2_sequence.txt" ],
			"1769-DPC-4" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-4_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-4_2_sequence.txt" ],
			"1769-DPC-5" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-5_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-5_2_sequence.txt" ]
		},
		"G2" => {
			"1769-DPC-10" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-10_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-10_2_sequence.txt" ],
			"1769-DPC-11" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-11_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-11_2_sequence.txt" ],
			"1769-DPC-13" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-13_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-13_2_sequence.txt" ],
			"1769-DPC-16" => [ "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-16_1_sequence.txt", "/scratch/cqs/guoy1/1769/rawdata/1769-DPC-16_2_sequence.txt" ]
		},
	},
	tophat2 => {
		target_dir => $tophat2_dir,
		option     => "--segment-length 25 -r 0 -p 8",
		batchmode  => 0,
		source     => "fastqfiles",
		pbs        => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cufflinks => {
		target_dir => "${target_dir}/cufflinks",
		option     => "-p 8",
		source     => "tophat2",
		pbs        => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cufflinks2 => {
		target_dir  => "${target_dir}/cufflinks2",
		option      => "-p 8",
		sourcefiles => {
			"G1" => {
				"1769-DPC-1" => "${tophat2_dir}/1769-DPC-1/accepted_hits.bam",
				"1769-DPC-3" => "${tophat2_dir}/1769-DPC-3/accepted_hits.bam",
				"1769-DPC-4" => "${tophat2_dir}/1769-DPC-4/accepted_hits.bam",
				"1769-DPC-5" => "${tophat2_dir}/1769-DPC-5/accepted_hits.bam",
			},
			"G2" => {
				"1769-DPC-10" => "${tophat2_dir}/1769-DPC-10/accepted_hits.bam",
				"1769-DPC-11" => "${tophat2_dir}/1769-DPC-11/accepted_hits.bam",
				"1769-DPC-13" => "${tophat2_dir}/1769-DPC-13/accepted_hits.bam",
				"1769-DPC-16" => "${tophat2_dir}/1769-DPC-16/accepted_hits.bam",
			},
		},
		pbs => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cuffmerge => {
		target_dir => $cuffmerges_dir,
		option     => "-p 8",
		source     => "cufflinks",
		pbs        => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cuffdiff => {
		target         => $cuffdiff_dir,
		option         => "-p 8 -N",
		transcript_gtf => $transcript_gtf,
		source         => "tophat2",
		pbs            => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cufflinks_cuffdiff => {
		target         => $cufflinks_cuffdiff_dir,
		option         => "-p 8 -N",
		transcript_gtf => "${cuffmerges_dir}/result/merged.gtf",
		source         => "tophat2",
		pbs            => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
};
