#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA3;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name                 => "kcv77",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/vickers/Vickers_20161222_smallRNA_v3_test",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 6,
  table_vis_group_text_size => 12,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 1,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_host_genome    => 1,
  search_miRBase        => 0,
  search_unmapped_reads => 0,
  blast_unmapped_reads  => 0,
  top_read_number       => 100,
  blast_top_reads       => 0,
  blast_localdb         => "/scratch/cqs/shengq1/references/blastdb",

  #NTA consideration
  consider_tRNA_NTA => 1,

  #next flex
  fastq_remove_random => 4,

  #Data
  files => {
    "Liver_WT_1"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i1_S1_R1_001.fastq.gz'],
  },
};

my $config = performSmallRNA_mm10( $def, 1 );

#performTask($config, "identical_NTA");

1;

