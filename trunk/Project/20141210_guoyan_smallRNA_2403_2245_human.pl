#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::SmallRNA;

my $def = {

	#General options
	task_name  => "2403_2245",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/smallRNA/20141210_guoyan_smallRNA_2403_2245_human/",
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
	bowtie1_miRBase_index => "/data/cqs/shengq1/reference/miRBase21/bowtie_index_1.0.1/mature.dna",

	#genome database
	mirna_coordinate    => "/data/cqs/shengq1/reference/miRBase20/hsa.gff3",
	trna_coordinate     => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed",
	trna_fasta          => "/data/cqs/guoy1/reference/smallrna/hg19_tRNA_ucsc_ensembl.bed.fa",
	smallrna_coordinate => "/data/cqs/guoy1/reference/smallrna/hg19_smallRNA_ucsc_ensembl.bed",
	bowtie1_index       => "/data/cqs/guoy1/reference/hg19/bowtie_index_hg19_rCRS_1.0.0/hg19_rCRS",

	#Data
	files => {
		"2403-CRF-01" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-1_1.fastq.gz"],
		"2403-CRF-02" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-2_1.fastq.gz"],
		"2403-CRF-03" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-3_1.fastq.gz"],
		"2403-CRF-04" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-4_1.fastq.gz"],
		"2403-CRF-05" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-5_1.fastq.gz"],
		"2403-CRF-06" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-6_1.fastq.gz"],
		"2403-CRF-07" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-7_1.fastq.gz"],
		"2403-CRF-08" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-8_1.fastq.gz"],
		"2403-CRF-09" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-9_1.fastq.gz"],
		"2403-CRF-10" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-10_1.fastq.gz"],
		"2403-CRF-11" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-11_1.fastq.gz"],
		"2403-CRF-12" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-12_1.fastq.gz"],
		"2403-CRF-13" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-13_1.fastq.gz"],
		"2403-CRF-14" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-14_1.fastq.gz"],
		"2403-CRF-15" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-15_1.fastq.gz"],
		"2403-CRF-18" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-18_1.fastq.gz"],
		"2403-CRF-19" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-19_1.fastq.gz"],
		"2403-CRF-20" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-20_1.fastq.gz"],
		"2403-CRF-21" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-21_1.fastq.gz"],
		"2403-CRF-22" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-22_1.fastq.gz"],
		"2403-CRF-23" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-23_1.fastq.gz"],
		"2403-CRF-24" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-24_1.fastq.gz"],
		"2403-CRF-26" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-26_1.fastq.gz"],
		"2403-CRF-27" => ["/scratch/cqs/guoy1/2403/2013-01-11/2403-CRF-27_1.fastq.gz"],
		"2245-CRF-00" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-0_1_sequence.txt.gz"],
		"2245-CRF-01" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-1_1_sequence.txt.gz"],
		"2245-CRF-02" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-2_1_sequence.txt.gz"],
		"2245-CRF-03" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-3_1_sequence.txt.gz"],
		"2245-CRF-04" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-4_1_sequence.txt.gz"],
		"2245-CRF-05" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-5_1_sequence.txt.gz"],
		"2245-CRF-06" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-6_1_sequence.txt.gz"],
		"2245-CRF-07" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-7_1_sequence.txt.gz"],
		"2245-CRF-08" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-8_1_sequence.txt.gz"],
		"2245-CRF-09" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-9_1_sequence.txt.gz"],
		"2245-CRF-10" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-10_1_sequence.txt.gz"],
		"2245-CRF-11" => ["/scratch/cqs/guoy1/2245/rawdata/2245-CRF-11_1_sequence.txt.gz"],
	  }
};

performSmallRNA($def);

1;

