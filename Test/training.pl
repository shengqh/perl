#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/home/shengq1/rnaseq/modules");

my $email = "quanhu.sheng\@vanderbilt.edu";

my $transcript_gtf = "/home/shengq1/rnaseq/references/mm10/annotation/Mus_musculus.GRCm38.68.gtf";

my $transcript_gtf_index = "/home/shengq1/rnaseq/references/mm10/annotation/index";

my $bowtie2_index = "/home/shengq1/rnaseq/references/mm10/bowtie2_index/mm10";

my $bwa_fasta = "/home/shengq1/rnaseq/references/mm10/mm10.fa";

my $config = {
	#some parameters may be used globally
	general => { 
		bowtie2_index        => $bowtie2_index,
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => $transcript_gtf_index,
		path_file            => "/home/shengq1/bin/path.txt",
		task_name            => "tutorial"
	},
    #raw file definition, the name "fastqfiles" will be used in the following definition
	fastqfiles => {
		"S1" => ["/home/shengq1/rnaseq/rawdata/s1_sequence.txt"],
		"S2" => ["/home/shengq1/rnaseq/rawdata/s2_sequence.txt"],
		"S3" => ["/home/shengq1/rnaseq/rawdata/s3_sequence.txt"],
		"S4" => ["/home/shengq1/rnaseq/rawdata/s4_sequence.txt"],
		"S5" => ["/home/shengq1/rnaseq/rawdata/s5_sequence.txt"],
		"S6" => ["/home/shengq1/rnaseq/rawdata/s6_sequence.txt"],
	},
	#group definition, which will be used at comparison analysis, one raw file can be put into multiple groups based on comparison purpose
	groups => {
		"CONTROL" => [ "S1", "S2", "S3" ],
		"SAMPLE"  => [ "S4", "S5", "S6" ],
	},
	#pair definition, which will be used at comparison analysis
	pairs  => { "SAMPLE_vs_CONTROL" => [ "SAMPLE", "CONTROL" ], },
	#fastqc definition
	fastqc => {
		target_dir => "${target_dir}/fastqc",
		#if you have any additional option...
		option     => "",
		#source files by link
		source_ref => "fastqfiles",
		#pbs parameter
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
	#duplicated fastqc definition, define source file in section itself rather than link 
    fastqc_dup => {
        target_dir => "${target_dir}/fastqc",
        option     => "",
        #source files directly
        source => {
            "S1" => ["/home/shengq1/rnaseq/rawdata/s1_sequence.txt"],
            "S2" => ["/home/shengq1/rnaseq/rawdata/s2_sequence.txt"],
            "S3" => ["/home/shengq1/rnaseq/rawdata/s3_sequence.txt"],
            "S4" => ["/home/shengq1/rnaseq/rawdata/s4_sequence.txt"],
            "S5" => ["/home/shengq1/rnaseq/rawdata/s5_sequence.txt"],
            "S6" => ["/home/shengq1/rnaseq/rawdata/s6_sequence.txt"],
        },
        pbs        => {
            "email"    => $email,
            "nodes"    => "1:ppn=1",
            "walltime" => "2",
            "mem"      => "10gb"
        },
    },
	bwa => {
		target_dir   => "${target_dir}/bwa",
		option       => "-q 15 -t 8",
		option_samse => "",
		source_ref   => "fastqfiles",
		fasta_file   => $bwa_fasta,
		pbs          => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "24",
			"mem"      => "20gb"
		},
	},
	gene_comparison_tophat2 => {
		target_dir           => "${target_dir}/gene_comparison_tophat2",
		option               => "-p 8",
		batchmode            => 0,
		source_ref           => "fastqfiles",
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => $transcript_gtf_index,
		pbs                  => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "240",
			"mem"      => "40gb"
		},
	},
	gene_comparison_cuffdiff => {
		target_dir     => "${target_dir}/gene_comparison_cuffdiff",
		option         => "-p 8",
		transcript_gtf => $transcript_gtf,
		source_ref     => "gene_comparison_tophat2",
		groups_ref     => "groups",
		pairs_ref      => "pairs",
		pbs            => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "720",
			"mem"      => "40gb"
		},
	},
	splicing_comparison_tophat2 => {
		target_dir => "${target_dir}/splicing_comparison_tophat2",
		option     => "-p 8",
		batchmode  => 0,
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "240",
			"mem"      => "40gb"
		},
	},
	splicing_comparison_cufflinks => {
		target_dir     => "${target_dir}/splicing_comparison_cufflinks",
		option         => "-p 8",
		source_ref     => "splicing_comparison_tophat2",
		transcript_gtf => $transcript_gtf,
		pbs            => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "10gb"
		},
	},
	splicing_comparison_cuffmerge => {
		target_dir => "${target_dir}/splicing_comparison_cuffmerge",
		option     => "-p 8",
		source_ref => "splicing_comparison_cufflinks",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	splicing_comparison_cuffdiff => {
		target_dir         => "${target_dir}/splicing_comparison_cuffdiff",
		option             => "-p 8",
		transcript_gtf_ref => "splicing_comparison_cuffmerge",
		source_ref         => "splicing_comparison_tophat2",
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

####comparison
tophat2_by_pbs( $config, "gene_comparison_tophat2" );

cuffdiff_by_pbs( $config, "gene_comparison_cuffdiff" );

####splicing
tophat2_by_pbs( $config, "splicing_comparison_tophat2" );

cufflinks_by_pbs( $config, "splicing_comparison_cufflinks" );

cuffmerge_by_pbs( $config, "splicing_comparison_cuffmerge" );

cuffdiff_by_pbs( $config, "splicing_comparison_cuffdiff" );

1;
