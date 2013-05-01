#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/dnaseq/VANGARD00061_jennifer");
my $bwa_dir = "${target_dir}/bwa_git";

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
	general => {
		bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
		transcript_gtf       => $transcript_gtf,
		transcript_gtf_index => $transcript_gtf_index,
		path_file            => "/home/shengq1/local/bin/path.txt",
		task_name            => "VANGARD00061"
	},
	fastqfiles => {
		"2510-DH-2" => [ "/blue/sequencer/Runs/projects/2510-DH/2013-04-26/2510-DH-2_1.fastq.gz", "/blue/sequencer/Runs/projects/2510-DH/2013-04-26/2510-DH-2_2.fastq.gz" ],
		"2510-DH-3" => [ "/blue/sequencer/Runs/projects/2510-DH/2013-04-26/2510-DH-3_1.fastq.gz", "/blue/sequencer/Runs/projects/2510-DH/2013-04-26/2510-DH-3_2.fastq.gz" ]
	},
	groups => {
		"NORMAL" => ["2510-DH-2"],
		"TUMOR"  => ["2510-DH-3"]
	},
	pairs  => { "TUMOR_vs_NORMAL" => [ "TUMOR", "NORMAL" ] },
	fastqc => {
		target_dir => "${target_dir}/fastqc",
		option     => "",
		source_ref => "fastqfiles",
		pbs        => {
			"email"    => $email,
			"nodes"    => "1:ppn=1",
			"walltime" => "2",
			"mem"      => "10gb"
		},
	},
	bwa => {
		target_dir      => $bwa_dir,
		option          => "-q 15 -t 8",
		option_sampe    => "",
		source_ref      => "fastqfiles",
		fasta_file      => "/data/cqs/shengq1/reference/hg19/hg19_chr.fa",
		estimate_insert => 1,
		source_ref      => "fastqfiles",
		pbs             => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "24",
			"mem"      => "20gb"
		},
	},
	refine => {
		target_dir => $bwa_dir,
		option     => "-Xmx40g",
		thread_count => 8,
		fasta_file => "/data/cqs/guoy1/reference/hg19/hg19_chr.fa",
		vcf_files   => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
		gatk_jar   => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
		markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
		source     => {
			"2510-DH-2" => ["${bwa_dir}/result/2510-DH-2/2510-DH-2_sort.bam"],
			"2510-DH-3" => ["${bwa_dir}/result/2510-DH-3/2510-DH-3_sort.bam"],
		},
		pbs => {
			"email"    => $email,
			"nodes"    => "1:ppn=8",
			"walltime" => "72",
			"mem"      => "40gb"
		},
	},
	snpindel =>{
    target_dir => "${target_dir}/SNPindel",
    option     => "-Xmx40g",
    fasta_file => "/data/cqs/guoy1/reference/hg19/hg19_chr.fa",
    vcf_files   => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
    gatk_jar   => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    gatk_option => "-l INFO -A DepthOfCoverage -A AlleleBalance -G Standard -stand_call_conf 50.0 -stand_emit_conf 10.0 -mbq 20 -deletions 0.05 -dcov 1000 -nct 8",
    source     => {
      "2510-DH" => ["${bwa_dir}/result/2510-DH-2/2510-DH-2_sort.realigned.recal.rmdup.bam", "${bwa_dir}/result/2510-DH-3/2510-DH-3_sort.realigned.recal.rmdup.bam"],
    },
    pbs => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
	}
};

#fastqc_by_pbs( $config, "fastqc" );
bwa_by_pbs_double( $config, "bwa" );
refine_bam_file( $config, "refine" );
gatk_snpindel($config, "snpindel");

1;
