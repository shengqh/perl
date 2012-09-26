#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = "/scratch/cqs/shengq1/rnaseq/1769";

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";

my $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		#transcript_gtf       => $transcript_gtf,
		#transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
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
	fastqc => {
		target_dir => "${target_dir}/fastqc",
		option     => "",
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "2",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	tophat2 => {
		target_dir => "${target_dir}/tophat2",
		option     => "--segment-length 25 -r 0 -p 8",
		batchmode  => 0,
		source_ref => "fastqfiles",
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
		source_ref => "tophat2",
		pbs        => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cufflinks2 => {
		target_dir => "${target_dir}/cufflinks2",
		option     => "-p 8",
		source     => {
			"G1" => {
				"1769-DPC-1" => "${target_dir}/tophat2/result/1769-DPC-1/accepted_hits.bam",
				"1769-DPC-3" => "${target_dir}/tophat2/result/1769-DPC-3/accepted_hits.bam",
				"1769-DPC-4" => "${target_dir}/tophat2/result/1769-DPC-4/accepted_hits.bam",
				"1769-DPC-5" => "${target_dir}/tophat2/result/1769-DPC-5/accepted_hits.bam",
			},
			"G2" => {
				"1769-DPC-10" => "${target_dir}/tophat2/result/1769-DPC-10/accepted_hits.bam",
				"1769-DPC-11" => "${target_dir}/tophat2/result/1769-DPC-11/accepted_hits.bam",
				"1769-DPC-13" => "${target_dir}/tophat2/result/1769-DPC-13/accepted_hits.bam",
				"1769-DPC-16" => "${target_dir}/tophat2/result/1769-DPC-16/accepted_hits.bam",
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
		target_dir => "${target_dir}/cuffmerge",
		option     => "-p 8",
		source_ref => "cufflinks",
		pbs        => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cuffmerge2 => {
		target_dir => "${target_dir}/cuffmerge2",
		option     => "-p 8",
		source     => "${target_dir}/cuffmerge2/assemblies.txt",
		pbs        => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cuffdiff => {
		target_dir     => "${target_dir}/cuffdiff",
		option         => "-p 8 -N",
		transcript_gtf => $transcript_gtf,
		source_ref     => "tophat2",
		pbs            => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cuffdiff2 => {
		target_dir     => "${target_dir}/cuffdiff2",
		option         => "-p 8 -N",
		transcript_gtf => $transcript_gtf,
		source         => {
			"G1" => {
				"1769-DPC-1" => "${target_dir}/tophat2/result/1769-DPC-1/accepted_hits.bam",
				"1769-DPC-3" => "${target_dir}/tophat2/result/1769-DPC-3/accepted_hits.bam",
				"1769-DPC-4" => "${target_dir}/tophat2/result/1769-DPC-4/accepted_hits.bam",
				"1769-DPC-5" => "${target_dir}/tophat2/result/1769-DPC-5/accepted_hits.bam",
			},
			"G2" => {
				"1769-DPC-10" => "${target_dir}/tophat2/result/1769-DPC-10/accepted_hits.bam",
				"1769-DPC-11" => "${target_dir}/tophat2/result/1769-DPC-11/accepted_hits.bam",
				"1769-DPC-13" => "${target_dir}/tophat2/result/1769-DPC-13/accepted_hits.bam",
				"1769-DPC-16" => "${target_dir}/tophat2/result/1769-DPC-16/accepted_hits.bam",
			},
		},
		pbs => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
	cufflinks_cuffdiff => {
		target_dir         => "${target_dir}/cufflinks_cuffdiff",
		option             => "-p 8 -N",
		transcript_gtf_ref => "cuffmerge",
		source_ref         => "tophat2",
		pbs                => {
			"email"    => "quanhu.sheng\@vanderbilt.edu",
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20000mb"
		},
	},
};

#fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

#run cuffdiff directly
#cuffdiff_by_pbs( $config, "cuffdiff" );

#run cufflinks-cuffmerge-cuffdiff
#cufflinks_by_pbs( $config, "cufflinks" );

#cuffmerge_by_pbs( $config, "cuffmerge" );

#cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

#define source data directly
#cufflinks_by_pbs( $config, "cufflinks2" );
#cuffmerge_by_pbs( $config, "cuffmerge2" );
#cuffdiff_by_pbs( $config, "cuffdiff2" );

1;
