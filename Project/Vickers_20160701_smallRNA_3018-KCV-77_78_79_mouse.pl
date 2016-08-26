#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $def = {

	#General options
	task_name                 => "KCV-77_78_79",
	email                     => "quanhu.sheng\@vanderbilt.edu",
	target_dir                => "/scratch/cqs/zhaos/vickers/20160701_smallRNA_3018-KCV-77_78_79_mouse",
	max_thread                => 8,
	cqstools                  => "/home/shengq1/cqstools/CQS.Tools.exe",
	sequencetask_run_time     => 6,
	table_vis_group_text_size => 12,

	#Default software parameter (don't change it except you really know it)
	fastq_remove_N        => 0,
	remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
	search_unmapped_reads => 1,
	blast_unmapped_reads  => 1,
	fastq_remove_random   => 4,                                   #next flex

	#Data
	files => {
		"Liver_WT_1"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i1_S1_R1_001.fastq.gz'],
		"Urine_WT_3"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i10_S8_R1_001.fastq.gz'],
		"Bile_WT_1"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i15_S9_R1_001.fastq.gz'],
		"Bile_SRBIKO_2"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i16_S10_R1_001.fastq.gz'],
		"Bile_WT_3"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i17_S11_R1_001.fastq.gz'],
		"Bile_SRBIKO_4"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i18_S12_R1_001.fastq.gz'],
		"Bile_WT_5"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i19_S13_R1_001.fastq.gz'],
		"Liver_SRBIKO_2"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i2_S2_R1_001.fastq.gz'],
		"HDL_WT_1"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i28_S14_R1_001.fastq.gz'],
		"HDL_SRBIKO_2"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i29_S15_R1_001.fastq.gz'],
		"Liver_WT_3"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i3_S3_R1_001.fastq.gz'],
		"HDL_WT_3"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i30_S16_R1_001.fastq.gz'],
		"HDL_SRBIKO_6"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i31_S17_R1_001.fastq.gz'],
		"HDL_SRBIKO_7"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i32_S18_R1_001.fastq.gz'],
		"Liver_SRBIKO_4"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i4_S4_R1_001.fastq.gz'],
		"APOB_WT_1"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i42_S19_R1_001.fastq.gz'],
		"APOB_SRBIKO_2"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i43_S20_R1_001.fastq.gz'],
		"APOB_WT_3"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i44_S21_R1_001.fastq.gz'],
		"APOB_SRBIKO_4"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i45_S22_R1_001.fastq.gz'],
		"Liver_WT_5"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i5_S5_R1_001.fastq.gz'],
		"Urine_WT_1"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i8_S6_R1_001.fastq.gz'],
		"Urine_SRBIKO_2"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i9_S7_R1_001.fastq.gz'],
		"APOB_WT_8"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i1_S1_R1_001.fastq.gz'],
		"Urine_WT_5"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i11_S7_R1_001.fastq.gz'],
		"Urine_SRBIKO_6"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i12_S8_R1_001.fastq.gz'],
		"Urine_SRBIKO_7"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i13_S9_R1_001.fastq.gz'],
		"Urine_WT_8"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i14_S10_R1_001.fastq.gz'],
		"Urine_SRBIKO_9"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i15_S11_R1_001.fastq.gz'],
		"APOB_SRBIKO_9"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i2_S2_R1_001.fastq.gz'],
		"Bile_SRBIKO_6"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i20_S12_R1_001.fastq.gz'],
		"Bile_SRBIKO_7"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i21_S13_R1_001.fastq.gz'],
		"Bile_WT_8"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i22_S14_R1_001.fastq.gz'],
		"Bile_SRBIKO_9"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i23_S15_R1_001.fastq.gz'],
		"HDL_WT_8"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i33_S16_R1_001.fastq.gz'],
		"HDL_SRBIKO_9"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i34_S17_R1_001.fastq.gz'],
		"HDL_WT_10"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i35_S18_R1_001.fastq.gz'],
		"HDL_SRBIKO_11"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i36_S19_R1_001.fastq.gz'],
		"APOB_WT_5"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i46_S20_R1_001.fastq.gz'],
		"APOB_SRBIKO_6"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i47_S21_R1_001.fastq.gz'],
		"APOB_SRBIKO_7"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i48_S22_R1_001.fastq.gz'],
		"Liver_SRBIKO_6"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i6_S3_R1_001.fastq.gz'],
		"Liver_SRBIKO_7"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i7_S4_R1_001.fastq.gz'],
		"Liver_WT_8"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i8_S5_R1_001.fastq.gz'],
		"Liver_SRBIKO_9"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i9_S6_R1_001.fastq.gz'],
		"Liver_WT_10"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i10_S6_R1_001.fastq.gz'],
		"Liver_SRBIKO_11" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i11_S7_R1_001.fastq.gz'],
		"Liver_WT_12"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i12_S8_R1_001.fastq.gz'],
		"Liver_SRBIKO_13" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i13_S9_R1_001.fastq.gz'],
		"Liver_WT_14"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i14_S10_R1_001.fastq.gz'],
		"Urine_SRBIKO_11" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i16_S11_R1_001.fastq.gz'],
		"Urine_SRBIKO_13" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i17_S12_R1_001.fastq.gz'],
		"Urine_WT_14"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i18_S13_R1_001.fastq.gz'],
		"Bile_WT_10"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i24_S14_R1_001.fastq.gz'],
		"Bile_WT_12"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i25_S15_R1_001.fastq.gz'],
		"Bile_SRBIKO_13"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i26_S16_R1_001.fastq.gz'],
		"Bile_WT_14"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i27_S17_R1_001.fastq.gz'],
		"APOB_WT_10"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i3_S1_R1_001.fastq.gz'],
		"HDL_SRBIKO_4"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i37_S18_R1_001.fastq.gz'],
		"HDL_WT_5"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i38_S19_R1_001.fastq.gz'],
		"HDL_WT_12"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i39_S20_R1_001.fastq.gz'],
		"APOB_SRBIKO_11"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i4_S2_R1_001.fastq.gz'],
		"HDL_SRBIKO_13"   => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i40_S21_R1_001.fastq.gz'],
		"HDL_WT_14"       => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i41_S22_R1_001.fastq.gz'],
		"APOB_WT_12"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i5_S3_R1_001.fastq.gz'],
		"APOB_SRBIKO_13"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i6_S4_R1_001.fastq.gz'],
		"APOB_WT_14"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i7_S5_R1_001.fastq.gz'],

	},
	groups => {
		"APOB_SRBIKO"  => [ "APOB_SRBIKO_2",  "APOB_SRBIKO_4",  "APOB_SRBIKO_9",  "APOB_SRBIKO_6",  "APOB_SRBIKO_7",   "APOB_SRBIKO_11",  "APOB_SRBIKO_13" ],
		"APOB_WT"      => [ "APOB_WT_1",      "APOB_WT_3",      "APOB_WT_8",      "APOB_WT_5",      "APOB_WT_10",      "APOB_WT_12",      "APOB_WT_14" ],
		"Bile_SRBIKO"  => [ "Bile_SRBIKO_2",  "Bile_SRBIKO_4",  "Bile_SRBIKO_6",  "Bile_SRBIKO_7",  "Bile_SRBIKO_9",   "Bile_SRBIKO_13" ],
		"Bile_WT"      => [ "Bile_WT_1",      "Bile_WT_3",      "Bile_WT_5",      "Bile_WT_8",      "Bile_WT_10",      "Bile_WT_12",      "Bile_WT_14" ],
		"HDL_SRBIKO"   => [ "HDL_SRBIKO_2",   "HDL_SRBIKO_6",   "HDL_SRBIKO_7",   "HDL_SRBIKO_9",   "HDL_SRBIKO_11",   "HDL_SRBIKO_4",    "HDL_SRBIKO_13" ],
		"HDL_WT"       => [ "HDL_WT_1",       "HDL_WT_3",       "HDL_WT_8",       "HDL_WT_10",      "HDL_WT_5",        "HDL_WT_12",       "HDL_WT_14" ],
		"Liver_SRBIKO" => [ "Liver_SRBIKO_2", "Liver_SRBIKO_4", "Liver_SRBIKO_6", "Liver_SRBIKO_7", "Liver_SRBIKO_9",  "Liver_SRBIKO_11", "Liver_SRBIKO_13" ],
		"Liver_WT"     => [ "Liver_WT_1",     "Liver_WT_3",     "Liver_WT_5",     "Liver_WT_8",     "Liver_WT_10",     "Liver_WT_12",     "Liver_WT_14" ],
		"Urine_SRBIKO" => [ "Urine_SRBIKO_2", "Urine_SRBIKO_6", "Urine_SRBIKO_7", "Urine_SRBIKO_9", "Urine_SRBIKO_11", "Urine_SRBIKO_13" ],
		"Urine_WT"     => [ "Urine_WT_3",     "Urine_WT_1",     "Urine_WT_5",     "Urine_WT_8",     "Urine_WT_14" ],

	},
	pairs => {
		"APOB_SRBIKO_vs_WT" => {
			groups => [ "APOB_WT", "APOB_SRBIKO" ],
			Batch  => [
				'3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-77',
				'3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79',
			],
		},
		"Bile_SRBIKO_vs_WT" => {
			groups => [ "Bile_WT", "Bile_SRBIKO" ],
			Batch  => [
				'3018-KCV-77', '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-77',
				'3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79',
			],
		},
		"HDL_SRBIKO_vs_WT" => {
			groups => [ "HDL_WT", "HDL_SRBIKO" ],
			Batch  => [
				'3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-77',
				'3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79',
			],
		},
		"Liver_SRBIKO_vs_WT" => {
			groups => [ "Liver_WT", "Liver_SRBIKO" ],
			Batch  => [
				'3018-KCV-77', '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-77',
				'3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79',
			],
		},
		"Urine_SRBIKO_vs_WT" => {
			groups => [ "Urine_WT",    "Urine_SRBIKO" ],
			Batch  => [ '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', ],
		},
	},
	groups_vis_layout => {
        "Row_Group" =>
          [ "HDL", "HDL", "APOB", "APOB", "LIVER", "LIVER", "BILE", "BILE","URINE","URINE" ],
        "Col_Group" => [
            "WT", "SR-BI KO", "WT", "SR-BI KO","WT", "SR-BI KO",
            "WT", "SR-BI KO", "WT", "SR-BI KO"
        ],
        "Groups" => [
            "HDL_WT",   "HDL_SRBIKO",
            "APOB_WT",  "APOB_SRBIKO",
            "Liver_WT", "Liver_SRBIKO",
            "Bile_WT",  "Bile_SRBIKO",
            "Urine_WT",  "Urine_SRBIKO"
        ],
    },
    pairs_host_deseq2_vis_layout => {
        "Col_Group" => [
            "HDL",  "HDL",  "HDL",      "APOB",  "APOB",
            "APOB",  "LIVER", "LIVER", "LIVER", 
            "BILE", "BILE", "BILE","URINE", "URINE", "URINE"
        ],
        "Row_Group" => [
            "miRNA",           "tRNA",
            "Other small RNA", 
            "miRNA",           "tRNA",
            "Other small RNA", 
            "miRNA",           "tRNA",
            "Other small RNA", 
            "miRNA",           "tRNA",
            "Other small RNA",      "miRNA",           "tRNA",
            "Other small RNA"
        ],
        "Groups" => [
            "deseq2_miRNA_HDL_SRBIKO_vs_WT",
            "deseq2_tRNA_HDL_SRBIKO_vs_WT",
            "deseq2_otherSmallRNA_HDL_SRBIKO_vs_WT",
            "deseq2_miRNA_APOB_SRBIKO_vs_WT",
            "deseq2_tRNA_APOB_SRBIKO_vs_WT",
            "deseq2_otherSmallRNA_APOB_SRBIKO_vs_WT",
            "deseq2_miRNA_Liver_SRBIKO_vs_WT",
            "deseq2_tRNA_Liver_SRBIKO_vs_WT",
            "deseq2_otherSmallRNA_Liver_SRBIKO_vs_WT",
            "deseq2_miRNA_Bile_SRBIKO_vs_WT",
            "deseq2_tRNA_Bile_SRBIKO_vs_WT",
            "deseq2_otherSmallRNA_Bile_SRBIKO_vs_WT",
            "deseq2_miRNA_Urine_SRBIKO_vs_WT",
            "deseq2_tRNA_Urine_SRBIKO_vs_WT",
            "deseq2_otherSmallRNA_Urine_SRBIKO_vs_WT"
        ],
    },
    pairs_nonHostGroups_deseq2_vis_layout => {
        "Col_Group" => [
            "HDL",  "HDL",  "HDL",   
            "APOB",  "APOB","APOB",  
            "LIVER", "LIVER", "LIVER",
            "BILE", "BILE", "BILE",
            "URINE", "URINE", "URINE"
        ],
        "Row_Group" => [
            "Microbiome","Environment","Fungus",
            "Microbiome","Environment","Fungus",
            "Microbiome","Environment","Fungus",
            "Microbiome","Environment","Fungus",
            "Microbiome","Environment","Fungus"
        ],
        "Groups" => [
            "deseq2_bacteria_group1_HDL_Knockout_VS_WildType",
            "deseq2_bacteria_group2_HDL_Knockout_VS_WildType",
            "deseq2_fungus_group4_HDL_Knockout_VS_WildType",
            "deseq2_bacteria_group1_APOB_Knockout_VS_WildType",
            "deseq2_bacteria_group2_APOB_Knockout_VS_WildType",
            "deseq2_fungus_group4_APOB_Knockout_VS_WildType",
            "deseq2_bacteria_group1_Liver_Knockout_VS_WildType",
            "deseq2_bacteria_group2_Liver_Knockout_VS_WildType",
            "deseq2_fungus_group4_Liver_Knockout_VS_WildType",
            "deseq2_bacteria_group1_Bile_Knockout_VS_WildType",
            "deseq2_bacteria_group2_Bile_Knockout_VS_WildType",
            "deseq2_fungus_group4_Bile_Knockout_VS_WildType",
            "deseq2_bacteria_group1_Urine_Knockout_VS_WildType",
            "deseq2_bacteria_group2_Urine_Knockout_VS_WildType",
            "deseq2_fungus_group4_Urine_Knockout_VS_WildType"
        ],
    },
    pairs_nonHosttRNArRNA_deseq2_vis_layout => {
        "Col_Group" => [
            "HDL",  "HDL",  "HDL",   
            "APOB",  "APOB","APOB",  
            "LIVER", "LIVER", "LIVER",
            "BILE", "BILE", "BILE",
            "URINE", "URINE", "URINE"
        ],
        "Row_Group" => [
            "tRNA","rRNAL","rRNAS",
            "tRNA","rRNAL","rRNAS",
            "tRNA","rRNAL","rRNAS",
            "tRNA","rRNAL","rRNAS",
            "tRNA","rRNAL","rRNAS"
        ],
        "Groups" => [
            "deseq2_nonHost_tRna_HDL_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAL_HDL_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAS_HDL_Knockout_VS_WildType",
            "deseq2_nonHost_tRna_APOB_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAL_APOB_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAS_APOB_Knockout_VS_WildType",
            "deseq2_nonHost_tRna_Liver_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAL_Liver_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAS_Liver_Knockout_VS_WildType",
            "deseq2_nonHost_tRna_Bile_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAL_Bile_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAS_Bile_Knockout_VS_WildType",
            "deseq2_nonHost_tRna_Urine_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAL_Urine_Knockout_VS_WildType",
            "deseq2_nonHost_rRNAS_Urine_Knockout_VS_WildType"
        ],
    },
};

my $config = performSmallRNA_mm10($def, 0);

$config->{identical_sequence_count_table}{target_dir} = "/scratch/cqs/shengq1/vickers/20160701_smallRNA_3018-KCV-77_78_79_mouse/class_independent/identical_sequence_count_table";
$config->{deseq2_top100Reads}{target_dir} = "/scratch/cqs/shengq1/vickers/20160701_smallRNA_3018-KCV-77_78_79_mouse/class_independent/deseq2_top100Reads";
$config->{deseq2_top100Contigs}{target_dir} = "/scratch/cqs/shengq1/vickers/20160701_smallRNA_3018-KCV-77_78_79_mouse/class_independent/deseq2_top100Contigs";

#performTask($config, "identical_sequence_count_table");
#performTask($config, "deseq2_top100Reads");
performTask($config, "deseq2_top100Contigs");

1;

