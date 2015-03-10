#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-17-18",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150310_smallRNA_3018-KCV-17-18_LDL",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 1,

  #Data
  files => {
  "3018-KCV-17-21" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-21_GTTTCG_L007_R1_001.fastq.gz"],
  "3018-KCV-17-22" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-22_CGTACG_L007_R1_001.fastq.gz"],
  "3018-KCV-17-23" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-23_GAGTGG_L007_R1_001.fastq.gz"],
  "3018-KCV-17-24" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-24_GGTAGC_L007_R1_001.fastq.gz"],
  "3018-KCV-17-25" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-25_ACTGAT_L007_R1_001.fastq.gz"],
  "3018-KCV-17-31" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-31_CACGAT_L007_R1_001.fastq.gz"],
  "3018-KCV-17-32" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-32_CACTCA_L007_R1_001.fastq.gz"],
  "3018-KCV-17-33" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-33_CAGGCG_L007_R1_001.fastq.gz"],
  "3018-KCV-17-34" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-34_CATGGC_L007_R1_001.fastq.gz"],
  "3018-KCV-17-35" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-17_LDL_RNA/3018-KCV-17-35_CATTTT_L007_R1_001.fastq.gz"],
  "3018-KCV-18-26" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-26_ATGAGC_L008_R1_001.fastq.gz"],
  "3018-KCV-18-27" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-27_ATTCCT_L008_R1_001.fastq.gz"],
  "3018-KCV-18-28" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-28_CAAAAG_L008_R1_001.fastq.gz"],
  "3018-KCV-18-29" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-29_CAACTA_L008_R1_001.fastq.gz"],
  "3018-KCV-18-30" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-30_CACCGG_L008_R1_001.fastq.gz"],
  "3018-KCV-18-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-36_CCAACA_L008_R1_001.fastq.gz"],
  "3018-KCV-18-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-37_CGGAAT_L008_R1_001.fastq.gz"],
  "3018-KCV-18-38" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-38_CTAGCT_L008_R1_001.fastq.gz"],
  "3018-KCV-18-39" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-39_CTATAC_L008_R1_001.fastq.gz"],
  "3018-KCV-18-40" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-18_LDL_RNA/3018-KCV-18-40_CTCAGA_L008_R1_001.fastq.gz"],
    }
};

performSmallRNA_hg19($def);

1;

