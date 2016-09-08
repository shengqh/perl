#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::FileUtils;
use Pipeline::SmallRNAUtils;
use Pipeline::SmallRNA;

my $genome = hg19_genome();

my $root = create_directory_or_die("/scratch/cqs/shengq1/vickers/test_20160721_smallRNA_3018-KCV-59_human_repeat_10/");

my $userdef = {

  #General options
  email                     => "quanhu.sheng\@vanderbilt.edu",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 12,
  table_vis_group_text_size => 14,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 1,

  #General options
  task_name  => "repeat10",
  target_dir => $root,

  #Data
  files => {
    "i01" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i02" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i03" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i04" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i05" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i06" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i07" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i08" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i09" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "i10" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
  },
};

my $defWithContig = getSmallRNADefinition( $userdef, $genome );
performSmallRNA($defWithContig);

1;

