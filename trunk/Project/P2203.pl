#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/P2203");

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "P2203"
	},
	fastqfiles => {
		"P2203-01" => [ "/data/cqs/shengq1/2203/rawdata/2203-WE-1_1_sequence.txt", "/data/cqs/shengq1/2203/rawdata/2203-WE-1_2_sequence.txt" ],
		"P2203-02" => [ "/data/cqs/shengq1/2203/rawdata/2203-WE-2_1_sequence.txt", "/data/cqs/shengq1/2203/rawdata/2203-WE-2_2_sequence.txt" ],
		"P2203-03" => [ "/data/cqs/shengq1/2203/rawdata/2203-WE-3_1_sequence.txt", "/data/cqs/shengq1/2203/rawdata/2203-WE-3_2_sequence.txt" ],
		"P2203-04" => [ "/data/cqs/shengq1/2203/rawdata/2203-WE-4_1_sequence.txt", "/data/cqs/shengq1/2203/rawdata/2203-WE-4_2_sequence.txt" ],
	},
	groups => {
		"FLO-1"       => ["P2203-01"],
		"FLO-1_MLN"   => ["P2203-02"],
		"FLO-1_R"     => ["P2203-03"],
		"FLO-1_R_MLN" => ["P2203-04"],
	},
	pairs => {
		"MLN_vs_None"   => [ "FLO-1_MLN",   "FLO-1" ],
		"R_vs_None"     => [ "FLO-1_R",     "FLO-1" ],
		"R_MLN_vs_None" => [ "FLO-1_R_MLN", "FLO-1" ],
		"R_vs_MLN"      => [ "FLO-1_R",     "FLO-1_MLN" ],
		"R_MLN_vs_MLN"  => [ "FLO-1_R_MLN", "FLO-1_MLN" ],
		"R_MLN_vs_R"    => [ "FLO-1_R_MLN", "FLO-1_R" ]
	},
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
	bwa => {
		target_dir => "${target_dir}/bwa",
		option     => "",
		source_ref => "fastqfiles",
		fasta_file => "/data/cqs/guoy1/reference/hg19/hg19_chr.fa",
		source     => {
			"P2203-01" => [ "/data/cqs/shengq1/2203/rawdata/2203-WE-1_1_sequence.txt", "/data/cqs/shengq1/2203/rawdata/2203-WE-1_2_sequence.txt" ],
		},
		pbs => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "24",
			"mem"      => "20gb"
		},
	},
	tophat2 => {
		target_dir => "${target_dir}/tophat2",
#		option     => "--segment-length 25 -r 150 -p 8",
        option     => "--segment-length 25 -r 5 -p 8",
		batchmode  => 0,
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "240",
			"mem"      => "40gb"
		},
	},
	cufflinks => {
		target_dir     => "${target_dir}/cufflinks",
		option         => "-p 8 -u -N",
		source_ref     => "tophat2",
		transcript_gtf => $transcript_gtf,
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
			"mem"      => "40gb"
		},
	},
	cuffdiff => {
		target_dir     => "${target_dir}/cuffdiff",
		option         => "-p 8 -u -N",
		transcript_gtf => $transcript_gtf,
		source_ref     => "tophat2",
		groups_ref     => "groups",
		pairs_ref      => "pairs",
		pbs            => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "720",
			"mem"      => "40gb"
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
			"walltime" => "720",
			"mem"      => "40gb"
		},
	},
};

bwa_by_pbs_double( $config, "bwa" );

#fastqc_by_pbs( $config, "fastqc" );

#tophat2_by_pbs( $config, "tophat2" );

#cuffdiff_by_pbs( $config, "cuffdiff" );

####run cufflinks-cuffmerge-cuffdiff
#cufflinks_by_pbs( $config, "cufflinks" );

#cuffmerge_by_pbs( $config, "cuffmerge" );

#cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

1;
