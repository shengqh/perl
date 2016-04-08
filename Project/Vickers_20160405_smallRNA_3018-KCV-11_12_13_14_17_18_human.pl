#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

	#General options
	task_name                 => "KCV-11_12_13_14_17_18",
	email                     => "shilin.zhao\@vanderbilt.edu",
	target_dir                => "/scratch/cqs/zhaos/vickers/20160405_smallRNA_3018-KCV-11_12_13_14_17_18_human",
	max_thread                => 8,
	cqstools                  => "/home/shengq1/cqstools/CQS.Tools.exe",
	sequencetask_run_time     => 12,
	table_vis_group_text_size => 14,

	#Default software parameter (don't change it except you really know it)
	fastq_remove_N        => 0,
	remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
	search_unmapped_reads => 1,
	blast_unmapped_reads  => 1,

	#Data
	files => {
		"3018-KCV-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-1_ATCACG_L002_R1_001.fastq.gz"],
		"3018-KCV-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-2_CGATGT_L002_R1_001.fastq.gz"],
		"3018-KCV-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-3_TTAGGC_L002_R1_001.fastq.gz"],
		"3018-KCV-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-4_TGACCA_L002_R1_001.fastq.gz"],
		"3018-KCV-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-5_ACAGTG_L002_R1_001.fastq.gz"],
		"3018-KCV-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-6_GCCAAT_L003_R1_001.fastq.gz"],
		"3018-KCV-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-7_CAGATC_L003_R1_001.fastq.gz"],
		"3018-KCV-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-8_ACTTGA_L003_R1_001.fastq.gz"],
		"3018-KCV-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-9_GATCAG_L003_R1_001.fastq.gz"],
		"3018-KCV-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-10_TAGCTT_L003_R1_001.fastq.gz"],
		"3018-KCV-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-11_GGCTAC_L002_R1_001.fastq.gz"],
		"3018-KCV-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-12_CTTGTA_L002_R1_001.fastq.gz"],
		"3018-KCV-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-13_AGTCAA_L002_R1_001.fastq.gz"],
		"3018-KCV-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-14_AGTTCC_L002_R1_001.fastq.gz"],
		"3018-KCV-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-15_ATGTCA_L002_R1_001.fastq.gz"],
		"3018-KCV-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-16_CCGTCC_L003_R1_001.fastq.gz"],
		"3018-KCV-17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-17_GTAGAG_L003_R1_001.fastq.gz"],
		"3018-KCV-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-18_GTCCGC_L003_R1_001.fastq.gz"],
		"3018-KCV-19" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-19_GTGAAA_L003_R1_001.fastq.gz"],
		"3018-KCV-20" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-20_GTGGCC_L003_R1_001.fastq.gz"],
		"3018-KCV-21" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-21_GTTTCG_L002_R1_001.fastq.gz"],
		"3018-KCV-22" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-22_CGTACG_L003_R1_001.fastq.gz"],
		"3018-KCV-23" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-23_GAGTGG_L002_R1_001.fastq.gz"],
		"3018-KCV-24" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-24_GGTAGC_L003_R1_001.fastq.gz"],

		"3018-KCV-13-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-1_ATCACG_L004_R1_001.fastq.gz"],
		"3018-KCV-13-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-2_CGATGT_L004_R1_001.fastq.gz"],
		"3018-KCV-13-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-3_TTAGGC_L004_R1_001.fastq.gz"],
		"3018-KCV-13-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-4_TGACCA_L004_R1_001.fastq.gz"],
		"3018-KCV-13-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-5_ACAGTG_L004_R1_001.fastq.gz"],
		"3018-KCV-14-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-6_GCCAAT_L005_R1_001.fastq.gz"],
		"3018-KCV-14-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-7_CAGATC_L005_R1_001.fastq.gz"],
		"3018-KCV-14-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-8_ACTTGA_L005_R1_001.fastq.gz"],
		"3018-KCV-14-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-9_GATCAG_L005_R1_001.fastq.gz"],
		"3018-KCV-14-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-10_TAGCTT_L005_R1_001.fastq.gz"],
		"3018-KCV-13-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-11_GGCTAC_L004_R1_001.fastq.gz"],
		"3018-KCV-13-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-12_CTTGTA_L004_R1_001.fastq.gz"],
		"3018-KCV-13-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-13_AGTCAA_L004_R1_001.fastq.gz"],
		"3018-KCV-13-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-14_AGTTCC_L004_R1_001.fastq.gz"],
		"3018-KCV-13-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-15_ATGTCA_L004_R1_001.fastq.gz"],
		"3018-KCV-14-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-16_CCGTCC_L005_R1_001.fastq.gz"],
		"3018-KCV-14-17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-17_GTAGAG_L005_R1_001.fastq.gz"],
		"3018-KCV-14-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-18_GTCCGC_L005_R1_001.fastq.gz"],
		"3018-KCV-14-19" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-19_GTGAAA_L005_R1_001.fastq.gz"],
		"3018-KCV-14-20" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-20_GTGGCC_L005_R1_001.fastq.gz"],

		"3018-KCV-17-21" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-21_GTTTCG_L007_R1_001.fastq.gz"],
		"3018-KCV-17-22" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-22_CGTACG_L007_R1_001.fastq.gz"],
		"3018-KCV-17-23" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-23_GAGTGG_L007_R1_001.fastq.gz"],
		"3018-KCV-17-24" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-24_GGTAGC_L007_R1_001.fastq.gz"],
		"3018-KCV-17-25" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-25_ACTGAT_L007_R1_001.fastq.gz"],
		"3018-KCV-17-31" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-31_CACGAT_L007_R1_001.fastq.gz"],
		"3018-KCV-17-32" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-32_CACTCA_L007_R1_001.fastq.gz"],
		"3018-KCV-17-33" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-33_CAGGCG_L007_R1_001.fastq.gz"],
		"3018-KCV-17-34" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-34_CATGGC_L007_R1_001.fastq.gz"],
		"3018-KCV-17-35" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-35_CATTTT_L007_R1_001.fastq.gz"],
		"3018-KCV-18-26" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-26_ATGAGC_L008_R1_001.fastq.gz"],
		"3018-KCV-18-27" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-27_ATTCCT_L008_R1_001.fastq.gz"],
		"3018-KCV-18-28" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-28_CAAAAG_L008_R1_001.fastq.gz"],
		"3018-KCV-18-29" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-29_CAACTA_L008_R1_001.fastq.gz"],
		"3018-KCV-18-30" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-30_CACCGG_L008_R1_001.fastq.gz"],
		"3018-KCV-18-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-36_CCAACA_L008_R1_001.fastq.gz"],
		"3018-KCV-18-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-37_CGGAAT_L008_R1_001.fastq.gz"],
		"3018-KCV-18-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-38_CTAGCT_L008_R1_001.fastq.gz"],
		"3018-KCV-18-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-39_CTATAC_L008_R1_001.fastq.gz"],
		"3018-KCV-18-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-40_CTCAGA_L008_R1_001.fastq.gz"],

	},
	groups => {
		"HetFH_HDL"   => [ "3018-KCV-01","3018-KCV-02", "3018-KCV-03", "3018-KCV-04", "3018-KCV-05", "3018-KCV-06", "3018-KCV-07", "3018-KCV-08", "3018-KCV-09", "3018-KCV-10" ],
		"Control_HDL" => [ "3018-KCV-13", "3018-KCV-15", "3018-KCV-16", "3018-KCV-17", "3018-KCV-18", "3018-KCV-19", "3018-KCV-20" ],

		"HetFH_VLDL"   => [ "3018-KCV-13-01","3018-KCV-13-02", "3018-KCV-13-03", "3018-KCV-13-04", "3018-KCV-13-05", "3018-KCV-14-06", "3018-KCV-14-07", "3018-KCV-14-08", "3018-KCV-14-09", "3018-KCV-14-10" ],
		"Control_VLDL" => [ "3018-KCV-13-13", "3018-KCV-13-15", "3018-KCV-14-16", "3018-KCV-14-17", "3018-KCV-14-18", "3018-KCV-14-19", "3018-KCV-14-20" ],

		"HetFH_LDL"   => [ "3018-KCV-17-21","3018-KCV-17-22", "3018-KCV-17-23", "3018-KCV-17-24", "3018-KCV-17-25", "3018-KCV-18-26", "3018-KCV-18-27", "3018-KCV-18-28", "3018-KCV-18-29", "3018-KCV-18-30" ],
		"Control_LDL" => [ "3018-KCV-17-33", "3018-KCV-17-35", "3018-KCV-18-36", "3018-KCV-18-37", "3018-KCV-18-38", "3018-KCV-18-39", "3018-KCV-18-40" ],

	},
	pairs => {
		"HetFH_HDL_vs_Control_HDL"     => { groups => [ "Control_HDL",   "HetFH_HDL" ], },
		"HetFH_VLDL_vs_Control_VLDL"  => { groups => [ "Control_VLDL", "HetFH_VLDL" ], },
		"HetFH_LDL_vs_Control_LDL"      => { groups => [ "Control_LDL",   "HetFH_LDL" ], },
	}
};

performSmallRNA_hg19($def);

1;

