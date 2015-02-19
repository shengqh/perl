#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-9",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150219_smallRNA_3018-KCV-9_human",
  max_thread => 8,

  #Data
  files => {
    "3018-KCV-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-1_ATCACG_L005_R1_001.fastq.gz"],
    "3018-KCV-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-2_CGATGT_L005_R1_001.fastq.gz"],
    "3018-KCV-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-4_TGACCA_L005_R1_001.fastq.gz"],
    "3018-KCV-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-6_GCCAAT_L005_R1_001.fastq.gz"],
    "3018-KCV-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-8_ACTTGA_L005_R1_001.fastq.gz"],
    "3018-KCV-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-10_TAGCTT_L005_R1_001.fastq.gz"],
    "3018-KCV-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-11_GGCTAC_L005_R1_001.fastq.gz"],
    "3018-KCV-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-12_CTTGTA_L005_R1_001.fastq.gz"],
    "3018-KCV-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-14_AGTTCC_L005_R1_001.fastq.gz"],
    "3018-KCV-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-15_ATGTCA_L005_R1_001.fastq.gz"],
    "3018-KCV-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-16_CCGTCC_L005_R1_001.fastq.gz"],
    "3018-KCV-17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-17_GTAGAG_L005_R1_001.fastq.gz"],
    "3018-KCV-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-18_GTCCGC_L005_R1_001.fastq.gz"],
    "3018-KCV-19" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-9_Human/3018-KCV-19_GTGAAA_L005_R1_001.fastq.gz"]
  }
};

performSmallRNAHuman($def);

1;

