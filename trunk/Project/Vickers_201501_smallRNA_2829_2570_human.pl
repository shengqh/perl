#!/usr/bin/perl
use strict;
use warnings;

use Pipeline::SmallRNA;

my $def = {

	#General options
	task_name  => "3050",
	email      => "quanhu.sheng\@vanderbilt.edu",
	target_dir => "/scratch/cqs/shengq1/vickers/201501_smallRNA_2829_2570_human/",
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
    "01-18-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-018-Post_CTTGTA_L005_R1_001.fastq.gz"],
    "01-18-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-018-Pre_GGCTAC_L005_R1_001.fastq.gz"],
    "01-28-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-28-Post_CGATGT_L004_R1_001.fastq.gz"],
    "01-28-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-28-Pre_ATCACG_L004_R1_001.fastq.gz"],
    "01-29-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-29-Post_TGACCA_L005_R1_001.fastq.gz"],
    "01-29-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-29-Pre_TTAGGC_L005_R1_001.fastq.gz"],
    "01-31-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-031-Post_GGTAGC_L005_R1_001.fastq.gz"],
    "01-31-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-031-Pre_GAGTGG_L005_R1_001.fastq.gz"],
    "01-36-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-36-Post_GCCAAT_L004_R1_001.fastq.gz"],
    "01-36-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-36-Pre_ACAGTG_L004_R1_001.fastq.gz"],
    "01-61-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-061-Post_ATGAGC_L004_R1_001.fastq.gz"],
    "01-61-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/01-061-Pre_ACTGAT_L004_R1_001.fastq.gz"],
    "03-07-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-007-Post_CAAAAG_L005_R1_001.fastq.gz"],
    "03-07-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-007-Pre_ATTCCT_L005_R1_001.fastq.gz"],
    "03-11-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-011-Post_CACCGG_L004_R1_001.fastq.gz"],
    "03-11-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-011-Pre_CAACTA_L004_R1_001.fastq.gz"],
    "03-15-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-015-Post_CACTCA_L005_R1_001.fastq.gz"],
    "03-15-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-015-Pre_CACGAT_L005_R1_001.fastq.gz"],
    "03-16-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-16-Post_ACTTGA_L005_R1_001.fastq.gz"],
    "03-16-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-16-Pre_CAGATC_L005_R1_001.fastq.gz"],
    "03-17-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-17-Post_TAGCTT_L004_R1_001.fastq.gz"],
    "03-17-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-17-Pre_GATCAG_L004_R1_001.fastq.gz"],
    "03-18-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-018-Post_CGTACG_L004_R1_001.fastq.gz"],
    "03-18-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-018-Pre_GTTTCG_L004_R1_001.fastq.gz"],
    "03-26-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-026-Post_AGTTCC_L004_R1_001.fastq.gz"],
    "03-26-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-026-Pre_AGTCAA_L004_R1_001.fastq.gz"],
    "03-31-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-031-Post_CATGGC_L004_R1_001.fastq.gz"],
    "03-31-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-031-Pre_CAGGCG_L004_R1_001.fastq.gz"],
    "03-33-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-033-Post_CCAACA_L005_R1_001.fastq.gz"],
    "03-33-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-033-Pre_CATTTT_L005_R1_001.fastq.gz"],
    "03-36-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-036-Post_CCGTCC_L005_R1_001.fastq.gz"],
    "03-36-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-036-Pre_ATGTCA_L005_R1_001.fastq.gz"],
    "03-47-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-047-Post_CTAGCT_L004_R1_001.fastq.gz"],
    "03-47-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-047-Pre_CGGAAT_L004_R1_001.fastq.gz"],
    "03-49-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-049-Post_CTCAGA_L005_R1_001.fastq.gz"],
    "03-49-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-049-Pre_CTATAC_L005_R1_001.fastq.gz"],
    "03-63-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-063-Post_GTGGCC_L005_R1_001.fastq.gz"],
    "03-63-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-063-Pre_GTGAAA_L005_R1_001.fastq.gz"],
    "03-65-Post"  => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-065-Post_GTCCGC_L004_R1_001.fastq.gz"],
    "03-65-Pre"   => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2570/03-065-Pre_GTAGAG_L004_R1_001.fastq.gz"],
    "2829-KCV-1A" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1A_ATCACG_L008_R1_001.fastq.gz"],
    "2829-KCV-1B" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1B_CGATGT_L008_R1_001.fastq.gz"],
    "2829-KCV-1C" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1C_TTAGGC_L008_R1_001.fastq.gz"],
    "2829-KCV-1D" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1D_TGACCA_L008_R1_001.fastq.gz"],
    "2829-KCV-1E" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1E_ACAGTG_L008_R1_001.fastq.gz"],
    "2829-KCV-1F" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1F_GCCAAT_L008_R1_001.fastq.gz"],
    "2829-KCV-1G" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1G_CAGATC_L008_R1_001.fastq.gz"],
    "2829-KCV-1H" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1H_ACTTGA_L008_R1_001.fastq.gz"],
    "2829-KCV-1I" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1I_GATCAG_L008_R1_001.fastq.gz"],
    "2829-KCV-1J" => ["/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201501_smallRNA_2829_2570_human/data/2829/2829-KCV-1J_TAGCTT_L008_R1_001.fastq.gz"],
  },
};

performSmallRNA($def);

1;

