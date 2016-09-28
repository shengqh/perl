#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name             => "3018-KCV-84",
  email                 => "shilin.zhao\@vanderbilt.edu",
  emailType             => "FAIL",
  target_dir            => "/scratch/cqs/shengq1/temp/20160928_smallRNA_3018-KCV_84_human",
  max_thread            => 8,
  cqstools              => "/home/shengq1/cqstools/CQS.Tools.exe",
  sequencetask_run_time => 7,
  
  #time cost task
  blast_top_reads         => 0,
  blast_unmapped_reads    => 0,
  perform_contig_analysis => 0,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  fastq_remove_random   => 4,      #next flex
  #Data
  files => {
  "CtrMinusHDL1" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i17_S1_R1_001.fastq.gz"],
  "CtrPlusHDL1" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i18_S2_R1_001.fastq.gz"],
  "CtrMinusHDL2" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i19_S3_R1_001.fastq.gz"],
  "CtrMinusHDL3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i21_S4_R1_001.fastq.gz"],
  "CtrPlusHDL2" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i22_S5_R1_001.fastq.gz"],
  "CtrMinusHDL4" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i23_S6_R1_001.fastq.gz"],
  "CtrPlusHDL3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i24_S7_R1_001.fastq.gz"],
  "SLEMinusHDL1" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i25_S8_R1_001.fastq.gz"],
  "SLEPlusHDL1" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i26_S9_R1_001.fastq.gz"],
  "SLEMinusHDL2" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i27_S10_R1_001.fastq.gz"],
  "SLEPlusHDL2" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i28_S11_R1_001.fastq.gz"],
  "SLEMinusHDL3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i29_S12_R1_001.fastq.gz"],
  "SLEPlusHDL3" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i30_S13_R1_001.fastq.gz"],
  "SLEMinusHDL4" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i31_S14_R1_001.fastq.gz"],
  "SLEPlusHDL4" => ["/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-84/3018/3018-KCV-84-i32_S15_R1_001.fastq.gz"],

  },
  groups => {
    "CtrPlusHDL"    => [ "CtrPlusHDL1", "CtrPlusHDL2", "CtrPlusHDL3"],
    "CtrMinusHDL" => [ "CtrMinusHDL1", "CtrMinusHDL2", "CtrMinusHDL13", "CtrMinusHDL14" ],
    "SLEPlusHDL" => [ "SLEPlusHDL1", "SLEPlusHDL2", "SLEPlusHDL3", "SLEPlusHDL4" ],
    "SLEMinusHDL" => [ "SLEMinusHDL1", "SLEMinusHDL2", "SLEMinusHDL3","SLEMinusHDL4"],
  },
  pairs => {
    "CtrPlusHDL_VS_CtrMinusHDL"   => { groups => [ "CtrMinusHDL",  "CtrPlusHDL" ], },
    "SLEMinusHDL_VS_CtrMinusHDL" => { groups => [ "CtrMinusHDL", "SLEMinusHDL" ], },

    "SLEPlusHDL_VS_SLEMinusHDL" => { groups => [ "SLEMinusHDL", "SLEPlusHDL" ], },
    "SLEPlusHDL_VS_CtrPlusHDL"    => { groups => [ "CtrPlusHDL", "SLEPlusHDL" ], },

  },
};

my $config = performSmallRNA_hg19($def, 0);
performTask($config, "identical_sequence_count_table");

1;

