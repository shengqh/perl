#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name                 => "mouse_spikein",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/vickers/20161121_smallRNA_3018_85_spikein_run2/mouse",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 6,
  table_vis_group_text_size => 12,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 0,
  top_read_number       => 100,
  blast_top_reads       => 0,
  blast_localdb         => "/scratch/cqs/shengq1/references/blastdb",
  
  special_sequence_file => "/scratch/cqs/shengq1/vickers/20161121_smallRNA_3018_85_spikein_run2/spikein.txt",

  #next flex
  fastq_remove_random => 4,

  #Data
  files => {
    "WT_HDL_23" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i23_S1_R1_001.fastq.gz"],
    "SRBIKO_HDL_24" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i24_S2_R1_001.fastq.gz"],
    "WT_HDL_25" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i25_S3_R1_001.fastq.gz"],
    "SRBIKO_HDL_26" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i26_S4_R1_001.fastq.gz"],
    "WT_HDL_27" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i27_S5_R1_001.fastq.gz"],
    "SRBIKO_HDL_28" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i28_S6_R1_001.fastq.gz"],
    "SRBIKO_HDL_29" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i29_S7_R1_001.fastq.gz"],
    "WT_HDL_30" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i30_S8_R1_001.fastq.gz"],
    "SRBIKO_HDL_31" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i31_S9_R1_001.fastq.gz"],
    "WT_HDL_32" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i32_S10_R1_001.fastq.gz"],
    "SRBIKO_HDL_33" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i33_S11_R1_001.fastq.gz"],
    "WT_HDL_34" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i34_S12_R1_001.fastq.gz"],
    "SRBIKO_HDL_35" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i35_S13_R1_001.fastq.gz"],
    "WT_HDL_36" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i36_S14_R1_001.fastq.gz"],
  },

  groups => {
    "WT"  => [ "WT_HDL_23",  "WT_HDL_25",  "WT_HDL_27",  "WT_HDL_30",  "WT_HDL_32",  "WT_HDL_34", "WT_HDL_36" ],
    "SRBIKO" => [ "SRBIKO_HDL_24", "SRBIKO_HDL_26", "SRBIKO_HDL_28", "SRBIKO_HDL_29", "SRBIKO_HDL_31", "SRBIKO_HDL_33", "SRBIKO_HDL_35" ],
  },
  pairs => {
    "SRBIKO_vs_WT" => [ "WT", "SRBIKO" ],
  },
};

performSmallRNA_mm10($def);

1;

