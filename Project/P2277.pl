#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/P2277");

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "P2277"
	},
	fastqfiles => {
		"G1" => {
			"P2277-01"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-1_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-1_2_sequence.txt" ],
            "P2277-02"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-2_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-2_2_sequence.txt" ],
            "P2277-03"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-3_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-3_2_sequence.txt" ],
            "P2277-04"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-4_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-4_2_sequence.txt" ],
            "P2277-05"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-5_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-5_2_sequence.txt" ],
            "P2277-06"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-6_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-6_2_sequence.txt" ],
            "P2277-07"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-7_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-7_2_sequence.txt" ],
            "P2277-08"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-8_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-8_2_sequence.txt" ],
            "P2277-09"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-9_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-9_2_sequence.txt" ],
            "P2277-10"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-10_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-10_2_sequence.txt" ],
            "P2277-11"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-11_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-11_2_sequence.txt" ],
            "P2277-12"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-12_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-12_2_sequence.txt" ],
            "P2277-13"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-13_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-13_2_sequence.txt" ],
            "P2277-14"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-14_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-14_2_sequence.txt" ],
            "P2277-15"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-15_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-15_2_sequence.txt" ],
            "P2277-16"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-16_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-16_2_sequence.txt" ],
            "P2277-17"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-17_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-17_2_sequence.txt" ],
            "P2277-18"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-18_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-18_2_sequence.txt" ],
            "P2277-19"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-19_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-19_2_sequence.txt" ],
            "P2277-20"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-20_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-20_2_sequence.txt" ],
            "P2277-21"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-21_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-21_2_sequence.txt" ],
            "P2277-22"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-22_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-22_2_sequence.txt" ],
            "P2277-23"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-23_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-23_2_sequence.txt" ],
            "P2277-24"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-24_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-24_2_sequence.txt" ],
            "P2277-25"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-25_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-25_2_sequence.txt" ],
            "P2277-26"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-26_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-26_2_sequence.txt" ],
            "P2277-27"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-27_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-27_2_sequence.txt" ],
            "P2277-28"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-28_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-28_2_sequence.txt" ],
            "P2277-29"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-29_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-29_2_sequence.txt" ],
            "P2277-30"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-30_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-30_2_sequence.txt" ],
            "P2277-31"   => [ "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-31_1_sequence.txt",   "/data/cqs/guom1/2277Rexer_mRNA/2277-BR-31_2_sequence.txt" ],
		},
	},
	fastqc => {
		target_dir => "${target_dir}/fastqc",
		option     => "",
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => $email,
			"nodes"    => "2",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
	tophat2 => {
		target_dir => "${target_dir}/tophat2",
		option     => "--segment-length 25 -r 150 -p 8",
		batchmode  => 0,
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => $email,
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	cufflinks => {
		target_dir => "${target_dir}/cufflinks",
		option     => "-p 8",
		source_ref => "tophat2",
		pbs        => {
			"email"    => $email,
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "20gb"
		},
	},
	cuffmerge => {
		target_dir => "${target_dir}/cuffmerge",
		option     => "-p 8",
		source_ref => "cufflinks",
		pbs        => {
			"email"    => $email,
			"nodes"    => "8",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	cuffdiff => {
		target_dir     => "${target_dir}/cuffdiff",
		option         => "-p 8 -N",
		transcript_gtf => $transcript_gtf,
		source_ref     => "tophat2",
		pbs            => {
			"email"    => $email,
			"nodes"    => "8",
			"walltime" => "240",
			"mem"      => "40gb"
		},
	},
	cufflinks_cuffdiff => {
		target_dir         => "${target_dir}/cufflinks_cuffdiff",
		option             => "-p 8 -N",
		transcript_gtf_ref => "cuffmerge",
		source_ref         => "tophat2",
		pbs                => {
			"email"    => $email,
			"nodes"    => "8",
			"walltime" => "240",
			"mem"      => "40gb"
		},
	},
};

#fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

#run cuffdiff directly
#cuffdiff_by_pbs( $config, "cuffdiff" );

#run cufflinks-cuffmerge-cuffdiff
cufflinks_by_pbs( $config, "cufflinks" );

cuffmerge_by_pbs( $config, "cuffmerge" );

#cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

1;
