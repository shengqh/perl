#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name  => "3018-KCV-23-25",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150427_smallRNA_3018-KCV-23-25_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "Healthy_1_B06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-1_ATCACG_L006_R1_001.fastq.gz"],
    "Healthy_1_B07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-3_TTAGGC_L006_R1_001.fastq.gz"],
    "Healthy_1_B08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-4_TGACCA_L006_R1_001.fastq.gz"],
    "Healthy_1_B09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-6_GCCAAT_L006_R1_001.fastq.gz"],
    "Healthy_1_B10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-17_GTAGAG_L006_R1_001.fastq.gz"],
    "Healthy_1_B11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-19_GTGAAA_L006_R1_001.fastq.gz"],
    "Healthy_1_B12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-20_GTGGCC_L006_R1_001.fastq.gz"],
    "Healthy_2_B06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-22_CGTACG_L006_R1_001.fastq.gz"],
    "Healthy_2_B07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-26_ATGAGC_L006_R1_001.fastq.gz"],
    "Healthy_2_B08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-27_ATTCCT_L006_R1_001.fastq.gz"],
    "Healthy_2_B09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-30_CACCGG_L006_R1_001.fastq.gz"],
    "Healthy_2_B10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-37_CGGAAT_L006_R1_001.fastq.gz"],
    "Healthy_2_B11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-40_CTCAGA_L006_R1_001.fastq.gz"],
    "Healthy_2_B12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-23/3018-KCV-23-41_GACGAC_L006_R1_001.fastq.gz"],
    "Healthy_3_B06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-5_ACAGTG_L007_R1_001.fastq.gz"],
    "Healthy_3_B07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-7_CAGATC_L007_R1_001.fastq.gz"],
    "Healthy_3_B08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-8_ACTTGA_L007_R1_001.fastq.gz"],
    "Healthy_3_B09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-10_TAGCTT_L007_R1_001.fastq.gz"],
    "Healthy_3_B10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-14_AGTTCC_L007_R1_001.fastq.gz"],
    "Healthy_3_B11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-21_GTTTCG_L007_R1_001.fastq.gz"],
    "Healthy_3_B12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-23_GAGTGG_L007_R1_001.fastq.gz"],
    "HetFH_1_B06"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-24_GGTAGC_L007_R1_001.fastq.gz"],
    "HetFH_1_B07"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-25_ACTGAT_L007_R1_001.fastq.gz"],
    "HetFH_1_B08"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-28_CAAAAG_L007_R1_001.fastq.gz"],
    "HetFH_1_B09"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-29_CAACTA_L007_R1_001.fastq.gz"],
    "HetFH_1_B10"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-32_CACTCA_L007_R1_001.fastq.gz"],
    "HetFH_1_B11"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-33_CAGGCG_L007_R1_001.fastq.gz"],
    "HetFH_1_B12"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-24/3018-KCV-24-36_CCAACA_L007_R1_001.fastq.gz"],
    "HetFH_2_B06"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-2_CGATGT_L008_R1_001.fastq.gz"],
    "HetFH_2_B07"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-9_GATCAG_L008_R1_001.fastq.gz"],
    "HetFH_2_B08"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-11_GGCTAC_L008_R1_001.fastq.gz"],
    "HetFH_2_B09"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-12_CTTGTA_L008_R1_001.fastq.gz"],
    "HetFH_2_B10"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-13_AGTCAA_L008_R1_001.fastq.gz"],
    "HetFH_2_B11"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-15_ATGTCA_L008_R1_001.fastq.gz"],
    "HetFH_2_B12"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-16_CCGTCC_L008_R1_001.fastq.gz"],
    "HetFH_3_B06"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-18_GTCCGC_L008_R1_001.fastq.gz"],
    "HetFH_3_B07"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-31_CACGAT_L008_R1_001.fastq.gz"],
    "HetFH_3_B08"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-34_CATGGC_L008_R1_001.fastq.gz"],
    "HetFH_3_B09"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-35_CATTTT_L008_R1_001.fastq.gz"],
    "HetFH_3_B10"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-38_CTAGCT_L008_R1_001.fastq.gz"],
    "HetFH_3_B11"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-39_CTATAC_L008_R1_001.fastq.gz"],
    "HetFH_3_B12"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-25/3018-KCV-25-41_GACGAC_L008_R1_001.fastq.gz"],
  }
};

my $config = performSmallRNA_hg19($def, 0);
performTask($config, "identical_sequence_count_table");

1;

