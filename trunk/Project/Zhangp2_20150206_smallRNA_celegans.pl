#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::SmallRNA;

my $def = {

	#General options
	task_name  => "celegans",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/smallRNA/Zhangp2_20150206_smallRNA_celegans/",
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
	mirna_fasta           => "/data/cqs/shengq1/reference/miRBase21/mature.dna.fa",
	trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
	trna_fasta          => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed",
	bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",
	mirbase_count_option  => "-p cel",

	#genome database
	mirna_coordinate    => "/data/cqs/shengq1/reference/miRBase21/cel.gff3",
	bowtie1_index       => "/scratch/cqs/zhangp2/reference/wormbase/bowtie_index_1.1.0/Caenorhabditis_elegans.WBcel235.dna.toplevel",

	#Data
	files => {
		"3018-KCV-7-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-18_GTCCGC_L006_R1_001.fastq.gz"],
		"3018-KCV-7-19" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-19_GTGAAA_L006_R1_001.fastq.gz"],
		"3018-KCV-7-20" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-20_GTGGCC_L006_R1_001.fastq.gz"],
		"3018-KCV-7-21" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-21_GTTTCG_L006_R1_001.fastq.gz"],
		"3018-KCV-7-22" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-22_CGTACG_L006_R1_001.fastq.gz"],
		"3018-KCV-7-23" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-23_GAGTGG_L006_R1_001.fastq.gz"],
		"3018-KCV-7-24" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-24_GGTAGC_L006_R1_001.fastq.gz"],
		"3018-KCV-7-25" => ["/gpfs21/scratch/cqs/shengq1/vickers/201412_smallRNA_3018_batch7_human/raw/3018-KCV-7-25_ACTGAT_L006_R1_001.fastq.gz"],
	  }
};

performSmallRNA($def);

1;

