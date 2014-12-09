#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::SmallRNA;

my $def = {

	#General options
	task_name  => "3050",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/",
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

	#genome database
	trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed",
	trna_fasta          => "/data/cqs/guoy1/reference/smallrna/mm10_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/mm10_smallRNA_ucsc_ensembl.bed",
	bowtie1_index       => "/data/cqs/shengq1/reference/mm10/bowtie_index/mm10",

	#Data
	files => {
		"3050-KCV-4-25" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-4-25_ACTGAT_L006_R1_001.fastq.gz"],
		"3050-KCV-4-30" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-4-30_CACCGG_L006_R1_001.fastq.gz"],
		"3050-KCV-4-33" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-4-33_CAGGCG_L006_R1_001.fastq.gz"],
		"3050-KCV-4-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-4-36_CCAACA_L006_R1_001.fastq.gz"],
		"3050-KCV-4-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-4-39_CTATAC_L006_R1_001.fastq.gz"],
		"3050-KCV-4-42" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-4-42_TAATCG_L006_R1_001.fastq.gz"],
		"3050-KCV-5-26" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-5-RPI26_stroke70.fastq.gz"],
		"3050-KCV-5-29" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-5-RPI29_stroke73.fastq.gz"],
		"3050-KCV-5-31" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-5-RPI31_stroke75.fastq.gz"],
		"3050-KCV-5-32" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-5-RPI32_stroke76.fastq.gz"],
		"3050-KCV-5-35" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-5-RPI35_stroke79.fastq.gz"],
		"3050-KCV-5-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-5-RPI37_stroke81.fastq.gz"],
		"3050-KCV-5-41" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-5-RPI41_stroke85.fastq.gz"],
		"3050-KCV-6-27" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-6-27_ATTCCT_L008_R1_001.fastq.gz"],
		"3050-KCV-6-28" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-6-28_CAAAAG_L008_R1_001.fastq.gz"],
		"3050-KCV-6-34" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-6-34_CATGGC_L008_R1_001.fastq.gz"],
		"3050-KCV-6-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-6-38_CTAGCT_L008_R1_001.fastq.gz"],
		"3050-KCV-6-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-6-40_CTCAGA_L008_R1_001.fastq.gz"],
		"3050-KCV-6-43" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3050/raw/3050-KCV-6-43_TACAGC_L008_R1_001.fastq.gz"],
	},
};

performSmallRNA($def);

1;

