#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::SmallRNA;

my $def = {

	#General options
	task_name  => "3018_b6",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/",
	max_thread => 8,

	#Default software parameter (don't change it except you really know it)
	bowtie1_option_1mm         => "-a -m 100 --best --strata -v 1 -p 8",
	bowtie1_option_pm          => "-a -m 100 --best --strata -v 0 -p 8",
	mirnacount_option          => "-s",                                         #ignore score
	smallrnacount_option       => "-s --min_overlap 0.5 --length --sequence",
	mirna_overlap_count_option => "-s --min_overlap 0.5 --gtf_key miRNA",
	min_read_length            => 16,

	#Software and miRBase database options
	samtools              => "/home/shengq1/local/bin/samtools/samtools",
	cqstools              => "/home/shengq1/cqstools/CQS.Tools.exe",
	mirna_fasta           => "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa",
	mirna_coordinate      => "/data/cqs/shengq1/reference/miRBase20/mmu.gff3",
	bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",
	mirbase_count_option  => "-p mmu",

	#genome database
	trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed",
	trna_fasta          => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/mm10_smallRNA_ucsc_ensembl.bed",
	bowtie1_index       => "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10",

	#Data
	files => {
		"3018-KCV-6-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-8_ACTTGA_L005_R1_001.fastq.gz"],
		"3018-KCV-6-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-9_GATCAG_L005_R1_001.fastq.gz"],
		"3018-KCV-6-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-10_TAGCTT_L005_R1_001.fastq.gz"],
		"3018-KCV-6-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-11_GGCTAC_L005_R1_001.fastq.gz"],
		"3018-KCV-6-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-12_CTTGTA_L005_R1_001.fastq.gz"],
		"3018-KCV-6-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch6_mouse/raw/3018-KCV-6-13_AGTCAA_L005_R1_001.fastq.gz"],
	},
};

performSmallRNA($def);

1;

