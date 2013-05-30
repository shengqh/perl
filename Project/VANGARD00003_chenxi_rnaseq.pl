#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/2379");

my $transcript_gtf = "/data/cqs/guoy1/reference/mm10/mm10_annotation/Mus_musculus.GRCm38.68.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10",
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/mm10_GRCm38_68",
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "2379"
	},
	fastqfiles => {
		"2379_P1" => ["/scratch/cqs/shengq1/rnaseq/2379/rawdata/2379-JRG-1_1_sequence.txt"],
		"2379_P3" => ["/scratch/cqs/shengq1/rnaseq/2379/rawdata/2379-JRG-3_1_sequence.txt"],
		"2379_P5" => ["/scratch/cqs/shengq1/rnaseq/2379/rawdata/2379-JRG-5_1_sequence.txt"],
		"2379_N2" => ["/scratch/cqs/shengq1/rnaseq/2379/rawdata/2379-JRG-2_1_sequence.txt"],
		"2379_N4" => ["/scratch/cqs/shengq1/rnaseq/2379/rawdata/2379-JRG-4_1_sequence.txt"],
		"2379_N6" => ["/scratch/cqs/shengq1/rnaseq/2379/rawdata/2379-JRG-6_1_sequence.txt"],
	},
	groups => {
		"POSITIVE" => [ "2379_P1", "2379_P3", "2379_P5" ],
		"NEGATIVE" => [ "2379_N2", "2379_N4", "2379_N6" ],
	},
	pairs  => { "POSITIVE_vs_NEGATIVE" => [ "POSITIVE", "NEGATIVE" ], },
	fastqc => {
		target_dir => "${target_dir}/fastqc",
		option     => "",
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=2",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
	tophat2 => {
		target_dir => "${target_dir}/tophat2",
		option     => "--segment-length 25 -r 0 -p 6",
		batchmode  => 0,
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=6",
			"walltime" => "72",
			"mem"      => "30gb"
		},
	},
	cufflinks => {
		target_dir     => "${target_dir}/cufflinks",
		option         => "-p 8 -u -N",
		transcript_gtf => $transcript_gtf,
		source_ref     => "tophat2",
		pbs            => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "10gb"
		},
	},
	cuffmerge => {
		target_dir => "${target_dir}/cuffmerge",
		option     => "-p 8",
		source_ref => "cufflinks",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "30gb"
		},
	},
	cufflinks_cuffdiff => {
		target_dir         => "${target_dir}/cufflinks_cuffdiff",
		option             => "-p 8 -u -N",
		transcript_gtf_ref => "cuffmerge",
		source_ref         => "tophat2",
		groups_ref         => "groups",
		pairs_ref          => "pairs",
		pbs                => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "30gb"
		},
	},
};

#fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

#cufflinks_by_pbs( $config, "cufflinks" );

#cuffmerge_by_pbs( $config, "cuffmerge" );

#cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );
1;
