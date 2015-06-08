#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-26",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150608_smallRNA_3018-KCV-26_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "3018-KCV-26-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-1_ATCACG_L008_R1_001.fastq.gz"],
    "3018-KCV-26-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-2_CGATGT_L008_R1_001.fastq.gz"],
    "3018-KCV-26-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-3_TTAGGC_L008_R1_001.fastq.gz"],
    "3018-KCV-26-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-4_TGACCA_L008_R1_001.fastq.gz"],
    "3018-KCV-26-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-5_ACAGTG_L008_R1_001.fastq.gz"],
    "3018-KCV-26-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-6_GCCAAT_L008_R1_001.fastq.gz"],
    "3018-KCV-26-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-7_CAGATC_L008_R1_001.fastq.gz"],
    "3018-KCV-26-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-8_ACTTGA_L008_R1_001.fastq.gz"],
    "3018-KCV-26-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-9_GATCAG_L008_R1_001.fastq.gz"],
    "3018-KCV-26-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-10_TAGCTT_L008_R1_001.fastq.gz"],
    "3018-KCV-26-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-15_ATGTCA_L008_R1_001.fastq.gz"],
    "3018-KCV-26-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-16_CCGTCC_L008_R1_001.fastq.gz"],
    "3018-KCV-26-17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-26/3018-KCV-26-17_GTAGAG_L008_R1_001.fastq.gz"],
  },
};

performSmallRNA_hg19($def);

1;

