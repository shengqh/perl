#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name                 => "human_spikein",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/vickers/20161121_smallRNA_3018_85_spikein_run2/human",
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

  #NTA consideration
  consider_tRNA_NTA=>1,

  #next flex
  fastq_remove_random => 4,
  #Data
  files => {
    "BM23APOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i37_S15_R1_001.fastq.gz"],
    "BM23HDL"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i38_S16_R1_001.fastq.gz"],
    "BM23FP"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i39_S17_R1_001.fastq.gz"],
    "LM55APOB" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i40_S18_R1_001.fastq.gz"],
    "LM55HDL"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i41_S19_R1_001.fastq.gz"],
    "LM55FP"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i42_S20_R1_001.fastq.gz"],
    "RNASpike" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i43_S21_R1_001.fastq.gz"],
    "SeqSpike" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-85_2nd_run/3018/3018-KCV-85-i44_S22_R1_001.fastq.gz"],
  },
};

my $config = performSmallRNA_hg19($def, 0);
performTask($config, "identical_NTA");

1;

