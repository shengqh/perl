#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name                 => "smallrnaseq",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/pipelines/smallrnaseq2",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 6,
  table_vis_group_text_size => 12,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,

  #nextflex kit
  fastq_remove_random => 4,

  #time cost task
  blast_top_reads         => 1,
  blast_unmapped_reads    => 0,
  perform_contig_analysis => 0,

  #Data
  files => {
    "Liver_WT_1"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i1_S1_R1_001.fastq.gz'],
    "Liver_WT_3"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i3_S3_R1_001.fastq.gz'],
    "Liver_WT_5"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i5_S5_R1_001.fastq.gz'],
    "Liver_SRBIKO_2" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i2_S2_R1_001.fastq.gz'],
    "Liver_SRBIKO_4" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i4_S4_R1_001.fastq.gz'],
    "Liver_SRBIKO_6" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i6_S3_R1_001.fastq.gz'],
  },
  groups => {
    "Liver_SRBIKO" => [ "Liver_SRBIKO_2", "Liver_SRBIKO_4", "Liver_SRBIKO_6" ],
    "Liver_WT"     => [ "Liver_WT_1",     "Liver_WT_3",     "Liver_WT_5" ],
  },
  pairs => {
    "Liver_SRBIKO_vs_WT" => [ "Liver_WT", "Liver_SRBIKO" ],
  },
};

performSmallRNA_mm10($def);

1;

