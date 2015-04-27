#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-21-22",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150427_smallRNA_3018-KCV-21-22_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "3018-KCV-21-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-1_ATCACG_L004_R1_001.fastq.gz"],
    "3018-KCV-21-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-2_CGATGT_L004_R1_001.fastq.gz"],
    "3018-KCV-21-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-3_TTAGGC_L004_R1_001.fastq.gz"],
    "3018-KCV-21-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-4_TGACCA_L004_R1_001.fastq.gz"],
    "3018-KCV-21-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-5_ACAGTG_L004_R1_001.fastq.gz"],
    "3018-KCV-21-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-6_GCCAAT_L004_R1_001.fastq.gz"],
    "3018-KCV-21-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-7_CAGATC_L004_R1_001.fastq.gz"],
    "3018-KCV-21-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-8_ACTTGA_L004_R1_001.fastq.gz"],
    "3018-KCV-21-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-21/3018-KCV-21-9_GATCAG_L004_R1_001.fastq.gz"],
    "3018-KCV-22-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-10_TAGCTT_L005_R1_001.fastq.gz"],
    "3018-KCV-22-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-11_GGCTAC_L005_R1_001.fastq.gz"],
    "3018-KCV-22-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-12_CTTGTA_L005_R1_001.fastq.gz"],
    "3018-KCV-22-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-13_AGTCAA_L005_R1_001.fastq.gz"],
    "3018-KCV-22-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-14_AGTTCC_L005_R1_001.fastq.gz"],
    "3018-KCV-22-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-15_ATGTCA_L005_R1_001.fastq.gz"],
    "3018-KCV-22-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-16_CCGTCC_L005_R1_001.fastq.gz"],
    "3018-KCV-22-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-22/3018-KCV-22-18_GTCCGC_L005_R1_001.fastq.gz"],
  }
};

performSmallRNA_hg19($def);

1;

