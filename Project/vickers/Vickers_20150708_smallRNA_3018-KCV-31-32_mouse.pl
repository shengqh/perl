#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-31-32",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150708_smallRNA_3018-KCV-31-32_mouse",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "LDL_01_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i1_ATCACG_L001_R1_001.fastq.gz"],
    "LDL_02_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i2_CGATGT_L001_R1_001.fastq.gz"],
    "LDL_03_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i3_TTAGGC_L001_R1_001.fastq.gz"],
    "LDL_04_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i4_TGACCA_L001_R1_001.fastq.gz"],
    "LDL_05_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i5_ACAGTG_L001_R1_001.fastq.gz"],
    "LDL_06_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i6_GCCAAT_L001_R1_001.fastq.gz"],
    "LDL_07_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i7_CAGATC_L001_R1_001.fastq.gz"],
    "LDL_08_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i8_ACTTGA_L001_R1_001.fastq.gz"],
    "LDL_09_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i9_GATCAG_L001_R1_001.fastq.gz"],
    "LDL_10_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i10_TAGCTT_L001_R1_001.fastq.gz"],
    "LDL_11_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i11_GGCTAC_L001_R1_001.fastq.gz"],
    "LDL_12_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-31-i12_CTTGTA_L001_R1_001.fastq.gz"],
    "HDL_01_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i13_AGTCAA_L002_R1_001.fastq.gz"],
    "HDL_02_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i14_AGTTCC_L002_R1_001.fastq.gz"],
    "HDL_03_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i15_ATGTCA_L002_R1_001.fastq.gz"],
    "HDL_04_CETPAPOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i16_CCGTCC_L002_R1_001.fastq.gz"],
    "HDL_05_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i17_GTAGAG_L002_R1_001.fastq.gz"],
    "HDL_06_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i18_GTCCGC_L002_R1_001.fastq.gz"],
    "HDL_07_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i19_GTGAAA_L002_R1_001.fastq.gz"],
    "HDL_08_CETP"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i20_GTGGCC_L002_R1_001.fastq.gz"],
    "HDL_09_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i21_GTTTCG_L002_R1_001.fastq.gz"],
    "HDL_10_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i22_CGTACG_L002_R1_001.fastq.gz"],
    "HDL_11_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i23_GAGTGG_L002_R1_001.fastq.gz"],
    "HDL_12_WT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-32-i24_GGTAGC_L002_R1_001.fastq.gz"],
  }
};

performSmallRNA_mm10($def);

1;
