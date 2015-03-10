#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-13-14",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150310_smallRNA_3018-KCV-13-14_VLDL",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 1,

  #Data
  files => {
  "3018-KCV-1314-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-1_ATCACG_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-2_CGATGT_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-3_TTAGGC_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-4_TGACCA_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-5_ACAGTG_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-6_GCCAAT_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-7_CAGATC_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-8_ACTTGA_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-9_GATCAG_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-10_TAGCTT_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-11_GGCTAC_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-12_CTTGTA_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-13_AGTCAA_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-14_AGTTCC_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-13_VLDL_RNA/3018-KCV-13-15_ATGTCA_L004_R1_001.fastq.gz"],
  "3018-KCV-1314-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-16_CCGTCC_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-17_GTAGAG_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-18_GTCCGC_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-19" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-19_GTGAAA_L005_R1_001.fastq.gz"],
  "3018-KCV-1314-20" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-14_VLDL_RNA/3018-KCV-14-20_GTGGCC_L005_R1_001.fastq.gz"],
    }
};

performSmallRNA_hg19($def);

1;

