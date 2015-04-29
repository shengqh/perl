#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-11-12",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150224_smallRNA_3018-KCV-11-12_HDL",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 1,

  #Data
  files => {
    "3018-KCV-01" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-1_ATCACG_L002_R1_001.fastq.gz"],
    "3018-KCV-02" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-2_CGATGT_L002_R1_001.fastq.gz"],
    "3018-KCV-03" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-3_TTAGGC_L002_R1_001.fastq.gz"],
    "3018-KCV-04" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-4_TGACCA_L002_R1_001.fastq.gz"],
    "3018-KCV-05" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-5_ACAGTG_L002_R1_001.fastq.gz"],
    "3018-KCV-06" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-6_GCCAAT_L003_R1_001.fastq.gz"],
    "3018-KCV-07" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-7_CAGATC_L003_R1_001.fastq.gz"],
    "3018-KCV-08" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-8_ACTTGA_L003_R1_001.fastq.gz"],
    "3018-KCV-09" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-9_GATCAG_L003_R1_001.fastq.gz"],
    "3018-KCV-10" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-10_TAGCTT_L003_R1_001.fastq.gz"],
    "3018-KCV-11" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-11_GGCTAC_L002_R1_001.fastq.gz"],
    "3018-KCV-12" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-12_CTTGTA_L002_R1_001.fastq.gz"],
    "3018-KCV-13" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-13_AGTCAA_L002_R1_001.fastq.gz"],
    "3018-KCV-14" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-14_AGTTCC_L002_R1_001.fastq.gz"],
    "3018-KCV-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-15_ATGTCA_L002_R1_001.fastq.gz"],
    "3018-KCV-16" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-16_CCGTCC_L003_R1_001.fastq.gz"],
    "3018-KCV-17" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-17_GTAGAG_L003_R1_001.fastq.gz"],
    "3018-KCV-18" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-18_GTCCGC_L003_R1_001.fastq.gz"],
    "3018-KCV-19" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-19_GTGAAA_L003_R1_001.fastq.gz"],
    "3018-KCV-20" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-20_GTGGCC_L003_R1_001.fastq.gz"],
    "3018-KCV-21" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-21_GTTTCG_L002_R1_001.fastq.gz"],
    "3018-KCV-22" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-22_CGTACG_L003_R1_001.fastq.gz"],
    "3018-KCV-23" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-11_HDL/3018-KCV-11-23_GAGTGG_L002_R1_001.fastq.gz"],
    "3018-KCV-24" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-12_HDL/3018-KCV-12-24_GGTAGC_L003_R1_001.fastq.gz"],
    }
};

performSmallRNATask_hg19($def, "identical_sequence_table");


1;

