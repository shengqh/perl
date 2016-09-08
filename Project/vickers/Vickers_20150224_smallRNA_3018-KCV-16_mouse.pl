#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-16",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-16_mouse",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 1,

  #Data
  files => {
    "3018-KCV-16-1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-1_ATCACG_L005_R1_001.fastq.gz"],
    "3018-KCV-16-2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-2_CGATGT_L005_R1_001.fastq.gz"],
    "3018-KCV-16-3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-3_TTAGGC_L005_R1_001.fastq.gz"],
    "3018-KCV-16-4" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-4_TGACCA_L005_R1_001.fastq.gz"],
    "3018-KCV-16-5" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-5_ACAGTG_L005_R1_001.fastq.gz"],
    "3018-KCV-16-6" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-6_GCCAAT_L005_R1_001.fastq.gz"],
    "3018-KCV-16-7" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-7_CAGATC_L005_R1_001.fastq.gz"],
    "3018-KCV-16-8" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-8_ACTTGA_L005_R1_001.fastq.gz"],
    "3018-KCV-16-9" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-16_Mouse/3018-KCV-16-9_GATCAG_L005_R1_001.fastq.gz"],
  }
};

performSmallRNA_mm10($def);

1;
