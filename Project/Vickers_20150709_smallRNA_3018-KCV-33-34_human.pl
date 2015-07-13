#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-33-34",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-33-34_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "Prog_034c_00Baseline" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i7_CAGATC_L003_R1_001.fastq.gz"],
    "Prog_034c_06mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i8_ACTTGA_L003_R1_001.fastq.gz"],
    "Prog_034c_12mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i9_GATCAG_L003_R1_001.fastq.gz"],
    "Prog_034c_24mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i10_TAGCTT_L003_R1_001.fastq.gz"],
    "Prog_034c_36mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i11_GGCTAC_L003_R1_001.fastq.gz"],
    "Prog_038c_00Baseline" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i1_ATCACG_L003_R1_001.fastq.gz"],
    "Prog_038c_06mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i2_CGATGT_L003_R1_001.fastq.gz"],
    "Prog_038c_12mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i3_TTAGGC_L003_R1_001.fastq.gz"],
    "Prog_038c_24mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i45_TCATTC_L003_R1_001.fastq.gz"],
    "Prog_038c_36mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i5_ACAGTG_L003_R1_001.fastq.gz"],
    "Prog_038c_48mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-33-i6_GCCAAT_L003_R1_001.fastq.gz"],
    "Prog_045c_00Baseline" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i12_CTTGTA_L004_R1_001.fastq.gz"],
    "Prog_045c_06mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i13_AGTCAA_L004_R1_001.fastq.gz"],
    "Prog_045c_12mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i46_TCCCGA_L004_R1_001.fastq.gz"],
    "Prog_045c_24mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i15_ATGTCA_L004_R1_001.fastq.gz"],
    "Prog_045c_36mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i16_CCGTCC_L004_R1_001.fastq.gz"],
    "Prog_045c_48mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i17_GTAGAG_L004_R1_001.fastq.gz"],
    "Prog_059c_06mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i47_TCGAAG_L004_R1_001.fastq.gz"],
    "Prog_059c_12mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i48_TCGGCA_L004_R1_001.fastq.gz"],
    "Prog_059c_24mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i19_GTGAAA_L004_R1_001.fastq.gz"],
    "Prog_059c_36mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i20_GTGGCC_L004_R1_001.fastq.gz"],
    "Prog_059c_48mo"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-31-38/3018-KCV-34-i21_GTTTCG_L004_R1_001.fastq.gz"],
  },
};

performSmallRNA_hg19($def);

1;

