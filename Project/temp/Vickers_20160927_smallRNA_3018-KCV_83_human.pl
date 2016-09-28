#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name             => "3018-KCV-83",
  email                 => "shilin.zhao\@vanderbilt.edu",
  emailType             => "FAIL",
  target_dir            => "/scratch/cqs/shengq1/temp/20160927_smallRNA_3018-KCV_83_human",
  max_thread            => 8,
  cqstools              => "/home/shengq1/cqstools/CQS.Tools.exe",
  sequencetask_run_time => 7,
  
  #time cost task
  blast_top_reads         => 1,
  blast_unmapped_reads    => 1,
  perform_contig_analysis => 1,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,

  fastq_remove_random   => 4,                                   #next flex
  #Data
  files => {
    "CtrStim1"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i1_S1_R1_001.fastq.gz"],
    "CtrStimHDL3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i10_S9_R1_001.fastq.gz"],
    "SLEStim3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i11_S10_R1_001.fastq.gz"],
    "SLEStimHDL3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i12_S11_R1_001.fastq.gz"],
    "CtrStim3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i13_S12_R1_001.fastq.gz"],
    "CtrStimHDL4" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i14_S13_R1_001.fastq.gz"],
    "SLEStim4" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i15_S14_R1_001.fastq.gz"],
    "SLEStimHDL4" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i16_S15_R1_001.fastq.gz"],
    "CtrStimHDL1"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i2_S2_R1_001.fastq.gz"],
    "SLEStim1"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i3_S3_R1_001.fastq.gz"],
    "SLEStimHDL1"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i4_S4_R1_001.fastq.gz"],
    "CtrStim2"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i5_S5_R1_001.fastq.gz"],
    "CtrStimHDL2"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i6_S6_R1_001.fastq.gz"],
    "SLEStim2"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i7_S7_R1_001.fastq.gz"],
    "SLEStimHDL2"  => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-83/3018/3018-KCV-83-i8_S8_R1_001.fastq.gz"],

  },
  groups => {
    "CtrStim"    => [ "CtrStim1", "CtrStim2", "CtrStim3"],
    "SLEStim" => [ "SLEStim1", "SLEStim2", "SLEStim3", "SLEStim4" ],
    "CtrStimHDL" => [ "CtrStimHDL1", "CtrStimHDL2", "CtrStimHDL3", "CtrStimHDL4" ],
    "SLEStimHDL" => [ "SLEStimHDL1", "SLEStimHDL2", "SLEStimHDL3","SLEStimHDL4"],
  },
  pairs => {
    "SLEStim_VS_CtrStim"   => { groups => [ "CtrStim",  "SLEStim" ], },
    "SLEStimHDL_VS_CtrStimHDL" => { groups => [ "CtrStimHDL", "SLEStimHDL" ], },
    "CtrStimHDL_VS_CtrStim" => { groups => [ "CtrStim", "CtrStimHDL" ], },
    "SLEStimHDL_VS_SLEStim"    => { groups => [ "SLEStim", "SLEStimHDL" ], },
  },
};

my $config = performSmallRNA_hg19($def, 1);
#performTask($config, "reads_mapping_summary");

1;

