#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-19-20",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150409_smallRNA_3018-KCV-19-20_mouse",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "3018-19-20-Citrate-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-1_ATCACG_L004_R1_001.fastq.gz"],
    "3018-19-20-Citrate-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-2_CGATGT_L005_R1_001.fastq.gz"],
    "3018-19-20-Citrate-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-3_TTAGGC_L004_R1_001.fastq.gz"],
    "3018-19-20-Citrate-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-4_TGACCA_L005_R1_001.fastq.gz"],
    "3018-19-20-Citrate-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-5_ACAGTG_L004_R1_001.fastq.gz"],
    "3018-19-20-Citrate-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-6_GCCAAT_L005_R1_001.fastq.gz"],
    "3018-19-20-Citrate-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-7_CAGATC_L004_R1_001.fastq.gz"],
    "3018-19-20-Citrate-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-8_ACTTGA_L005_R1_001.fastq.gz"],
    "3018-19-20-Citrate-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-9_GATCAG_L004_R1_001.fastq.gz"],
    "3018-19-20-Citrate-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-10_TAGCTT_L005_R1_001.fastq.gz"],
    "3018-19-20-Citrate-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-11_GGCTAC_L004_R1_001.fastq.gz"],
    "3018-19-20-Citrate-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-12_CTTGTA_L005_R1_001.fastq.gz"],
    "3018-19-20-Citrate-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-13_AGTCAA_L004_R1_001.fastq.gz"],
    "3018-19-20-Citrate-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-14_AGTTCC_L004_R1_001.fastq.gz"],
    "3018-19-20-STZ-15"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-15_ATGTCA_L004_R1_001.fastq.gz"],
    "3018-19-20-STZ-16"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-16_CCGTCC_L005_R1_001.fastq.gz"],
    "3018-19-20-STZ-17"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-17_GTAGAG_L004_R1_001.fastq.gz"],
    "3018-19-20-STZ-19"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-19_GTGAAA_L005_R1_001.fastq.gz"],
    "3018-19-20-STZ-20"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-20_GTGGCC_L005_R1_001.fastq.gz"],
    "3018-19-20-STZ-22"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-22_CGTACG_L005_R1_001.fastq.gz"],
    "3018-19-20-STZ-23"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-23_GAGTGG_L004_R1_001.fastq.gz"],
    "3018-19-20-STZ-24"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-24_GGTAGC_L005_R1_001.fastq.gz"],
    "3018-19-20-STZ-25"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-19_Mouse_HDL/3018-KCV-19-25_ACTGAT_L004_R1_001.fastq.gz"],
    "3018-19-20-STZ-26"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-20_Mouse_HDL/3018-KCV-20-26_ATGAGC_L005_R1_001.fastq.gz"],
  }
};

performSmallRNA_mm10($def);

1;
