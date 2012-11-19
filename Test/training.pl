#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/home/shengq1/rnaseq/modules");

my $transcript_gtf =
  "/home/shengq1/rnaseq/references/mm10/annotation/Mus_musculus.GRCm38.68.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $bowtie2_index = "/home/shengq1/rnaseq/references/mm10/bowtie2_index/mm10";

my $bwa_fasta = "/home/shengq1/rnaseq/references/mm10/mm10.fa";

my $config = {
	general => {
		bowtie2_index  => $bowtie2_index,
		transcript_gtf => $transcript_gtf,
		transcript_gtf_index =>
		  "/home/shengq1/rnaseq/references/mm10/annotation/index",
		path_file => "/home/shengq1/bin/path.txt",
		task_name => "Tutorial"
	},
	fastqfiles => {
		"S1" => ["/home/shengq1/rnaseq/rawdata/s1_sequence.txt"],
		"S2" => ["/home/shengq1/rnaseq/rawdata/s2_sequence.txt"],
		"S3" => ["/home/shengq1/rnaseq/rawdata/s3_sequence.txt"],
		"S4" => ["/home/shengq1/rnaseq/rawdata/s4_sequence.txt"],
		"S5" => ["/home/shengq1/rnaseq/rawdata/s5_sequence.txt"],
		"S6" => ["/home/shengq1/rnaseq/rawdata/s6_sequence.txt"],
	},
	groups => {
		"CONTROL" => [ "S1", "S2", "S3" ],
		"SAMPLE"  => [ "S4", "S5", "S6" ],
	},
	pairs  => { "SAMPLE_vs_CONTROL" => [ "SAMPLE", "CONTROL" ], },
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
		target_dir      => "${target_dir}/bwa",
		option          => "-q 15 -t 8",
		option_samse    => "",
		source_ref      => "fastqfiles",
		fasta_file      => $bwa_fasta,
		source_ref      => "fastqfiles",
		pbs             => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "24",
			"mem"      => "20gb"
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
		option         => "-p 8",
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
		option         => "-p 8",
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
		option             => "-p 8",
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

bwa_by_pbs_single( $config, "bwa" );

fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

cuffdiff_by_pbs( $config, "cuffdiff" );

####run cufflinks-cuffmerge-cuffdiff
cufflinks_by_pbs( $config, "cufflinks" );

cuffmerge_by_pbs( $config, "cuffmerge" );

cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

1;
