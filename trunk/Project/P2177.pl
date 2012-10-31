#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/P2177_2");

my $transcript_gtf = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68",
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "P2177"
	},
	fastqfiles => {
		"P2177-01" => [ "/data/cqs/shengq1/2177/2177-WE-1_1_sequence.txt",  "/data/cqs/shengq1/2177/2177-WE-1_2_sequence.txt" ],
		"P2177-02" => [ "/data/cqs/shengq1/2177/2177-WE-2_1_sequence.txt",  "/data/cqs/shengq1/2177/2177-WE-2_2_sequence.txt" ],
		"P2177-15" => [ "/data/cqs/shengq1/2177/2177-WE-15_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-15_2_sequence.txt" ],
		"P2177-16" => [ "/data/cqs/shengq1/2177/2177-WE-16_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-16_2_sequence.txt" ],
		"P2177-17" => [ "/data/cqs/shengq1/2177/2177-WE-17_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-17_2_sequence.txt" ],
		"P2177-18" => [ "/data/cqs/shengq1/2177/2177-WE-18_1_sequence.txt", "/data/cqs/shengq1/2177/2177-WE-18_2_sequence.txt" ],
	},
	groups => {
		"MKN45"          => [ "P2177-01", "P2177-02" ],
		"AGS_CTRL"       => ["P2177-15"],
		"AGS_DP32"       => ["P2177-16"],
		"MKN45_SC"       => ["P2177-17"],
		"MKN45_DP_SHRNA" => ["P2177-18"],
	},
	pairs  => { "ALL" => [ "MKN45", "AGS_CTRL", "AGS_DP32", "MKN45_SC", "MKN45_DP_SHRNA" ] },
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
		option     => "--segment-length 25 -r 0 -p 8",
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
	cufflinks_NG => {
		target_dir => "${target_dir}/cufflinks_NG",
		option     => "-p 8 -u -N",
		source_ref => "tophat2",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "10gb"
		},
	},
	cuffmerge_NG => {
		target_dir => "${target_dir}/cuffmerge_NG",
		option     => "-p 8",
		source_ref => "cufflinks_NG",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	cufflinks_cuffdiff_NG => {
		target_dir         => "${target_dir}/cufflinks_cuffdiff_NG",
		option             => "-p 8 -u -N",
		transcript_gtf_ref => "cuffmerge_NG",
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
    cufflinks_NG_DEFAULT => {
        target_dir => "${target_dir}/cufflinks_NG_DEFAULT",
        option     => "-p 8",
        source_ref => "tophat2",
        pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=8",
            "walltime" => "72",
            "mem"      => "10gb"
        },
    },
    cuffmerge_NG_DEFAULT => {
        target_dir => "${target_dir}/cuffmerge_NG_DEFAULT",
        option     => "-p 8",
        source_ref => "cufflinks_NG_DEFAULT",
        pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=8",
            "walltime" => "72",
            "mem"      => "40gb"
        },
    },
    cufflinks_cuffdiff_NG_DEFAULT => {
        target_dir         => "${target_dir}/cufflinks_cuffdiff_NG_DEFAULT",
        option             => "-p 8",
        transcript_gtf_ref => "cuffmerge_NG_DEFAULT",
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

#fastqc_by_pbs( $config, "fastqc" );

#tophat2_by_pbs( $config, "tophat2" );

#cuffdiff_by_pbs( $config, "cuffdiff" );

#run cufflinks-cuffmerge-cuffdiff
cufflinks_by_pbs( $config, "cufflinks_NG" );

cuffmerge_by_pbs( $config, "cuffmerge_NG" );

cuffdiff_by_pbs( $config, "cufflinks_cuffdiff_NG" );

cufflinks_by_pbs( $config, "cufflinks_NG_DEFAULT" );

cuffmerge_by_pbs( $config, "cuffmerge_NG_DEFAULT" );

cuffdiff_by_pbs( $config, "cufflinks_cuffdiff_NG_DEFAULT" );

1;
