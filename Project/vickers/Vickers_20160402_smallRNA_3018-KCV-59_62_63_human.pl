#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name                 => "KCV-59_62_63",
  email                     => "quanhu.sheng\@vanderbilt.edu",
  target_dir                => "/scratch/cqs/shengq1/vickers/20160402_smallRNA_3018-KCV-59_62_63_human",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/CQS.Tools.exe",
  sequencetask_run_time     => 12,
  table_vis_group_text_size => 14,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 1,

  #Data
  files => {
    "1010_hem_propranolol_ShortTerm"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i1_S1_R1_001.fastq.gz'],
    "1097_hem_Baseline"                      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i10_S10_R1_001.fastq.gz'],
    "1114_hem_Baseline"                      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i13_S11_R1_001.fastq.gz'],
    "1033_hem_propranolol_LongTerm"          => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i14_S12_R1_001.fastq.gz'],
    "1051_hem_propranolol_ShortTerm_Sample2" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i15_S13_R1_001.fastq.gz'],
    "1053_hem_propranolol_LongTerm"          => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i16_S14_R1_001.fastq.gz'],
    "1077_hem_propranolol_ShortTerm"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i17_S15_R1_001.fastq.gz'],
    "1024_hem_propranolol_LongTerm"          => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i2_S2_R1_001.fastq.gz'],
    "1033_hem_propranolol_ShortTerm"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i3_S3_R1_001.fastq.gz'],
    "1051_hem_propranolol_ShortTerm_Sample1" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i4_S4_R1_001.fastq.gz'],
    "1053_hem_propranolol_ShortTerm"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i5_S5_R1_001.fastq.gz'],
    "1072_hem_propranolol_ShortTerm"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i6_S6_R1_001.fastq.gz'],
    "1077_hem_Baseline"                      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i7_S7_R1_001.fastq.gz'],
    "1078_hem_propranolol_ShortTerm"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i8_S8_R1_001.fastq.gz'],
    "1093_hem_propranolol_ShortTerm"         => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-59/3018-KCV-59-i9_S9_R1_001.fastq.gz'],

    "1010_hem_Serum_NoDrug"          => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i1_S1_R1_001.fastq.gz'],
    "1119_hem_Baseline"              => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i10_S10_R1_001.fastq.gz'],
    "1146_hem_Baseline_Sample1"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i11_S11_R1_001.fastq.gz'],
    "1024_hem_Serum_NoDrug"          => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i2_S2_R1_001.fastq.gz'],
    "1051_hem_propranolol_LongTerm"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i3_S3_R1_001.fastq.gz'],
    "1069_control_NoDrug"            => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i4_S4_R1_001.fastq.gz'],
    "1072_hem_propranolol_LongTerm"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i5_S5_R1_001.fastq.gz'],
    "1078_hem_propranolol_LongTerm"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i6_S6_R1_001.fastq.gz'],
    "1097_hem_propranolol_ShortTerm" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i7_S7_R1_001.fastq.gz'],
    "1114_hem_propranolol_LongTerm"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i8_S8_R1_001.fastq.gz'],
    "1118_hem_Baseline"              => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-62/3018-KCV-62-i9_S9_R1_001.fastq.gz'],

    "1114_hem_propranolol_ShortTerm" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i10_S9_R1_001.fastq.gz'],
    "1146_hem_Baseline_Sample2"      => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i11_S10_R1_001.fastq.gz'],
    "1147_control_NoDrug"            => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i12_S11_R1_001.fastq.gz'],
    "1010_hem_propranolol_LongTerm"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i2_S1_R1_001.fastq.gz'],
    "1024_hem_propranolol_ShortTerm" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i3_S2_R1_001.fastq.gz'],
    "1048_control_NoDrug"            => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i4_S3_R1_001.fastq.gz'],
    "1053_hem_Baseline"              => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i5_S4_R1_001.fastq.gz'],
    "1068_control_NoDrug"            => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i6_S5_R1_001.fastq.gz'],
    "1072_hem_Baseline"              => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i7_S6_R1_001.fastq.gz'],
    "1093_hem_Baseline"              => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i8_S7_R1_001.fastq.gz'],
    "1097_hem_propranolol_LongTerm"  => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-63/3018-KCV-63-i9_S8_R1_001.fastq.gz'],
  },
  groups => {
    "Control_NoDrug" => [ "1069_control_NoDrug", "1147_control_NoDrug", "1048_control_NoDrug", "1068_control_NoDrug" ],
    "Hem_Baseline"   => [
      "1097_hem_Baseline",         "1114_hem_Baseline", "1077_hem_Baseline", "1119_hem_Baseline", "1146_hem_Baseline_Sample1", "1118_hem_Baseline",
      "1146_hem_Baseline_Sample2", "1053_hem_Baseline", "1093_hem_Baseline"
    ],
    "Hem_propranolol_LongTerm" => [
      "1033_hem_propranolol_LongTerm", "1053_hem_propranolol_LongTerm", "1024_hem_propranolol_LongTerm", "1051_hem_propranolol_LongTerm",
      "1072_hem_propranolol_LongTerm", "1078_hem_propranolol_LongTerm", "1114_hem_propranolol_LongTerm", "1010_hem_propranolol_LongTerm",
      "1097_hem_propranolol_LongTerm"
    ],
    "Hem_propranolol_ShortTerm" => [
      "1010_hem_propranolol_ShortTerm",         "1051_hem_propranolol_ShortTerm_Sample2", "1077_hem_propranolol_ShortTerm", "1033_hem_propranolol_ShortTerm",
      "1051_hem_propranolol_ShortTerm_Sample1", "1053_hem_propranolol_ShortTerm",         "1072_hem_propranolol_ShortTerm", "1078_hem_propranolol_ShortTerm",
      "1093_hem_propranolol_ShortTerm",         "1097_hem_propranolol_ShortTerm",         "1114_hem_propranolol_ShortTerm", "1024_hem_propranolol_ShortTerm"
    ],
    "Hem_Serum_NoDrug" => [ "1010_hem_Serum_NoDrug", "1024_hem_Serum_NoDrug" ],

  },
  pairs => {
    "Hem_Baseline_VS_Control_NoDrug"                        => { groups => [ "Control_NoDrug",            "Hem_Baseline" ], },
    "Hem_propranolol_LongTerm_VS_Hem_Baseline"              => { groups => [ "Hem_Baseline",              "Hem_propranolol_LongTerm" ], },
    "Hem_propranolol_ShortTerm_VS_Hem_Baseline"             => { groups => [ "Hem_Baseline",              "Hem_propranolol_ShortTerm" ], },
    "Hem_propranolol_LongTerm_VS_Hem_propranolol_ShortTerm" => { groups => [ "Hem_propranolol_ShortTerm", "Hem_propranolol_LongTerm" ], },
    "Hem_Serum_NoDrug_VS_Hem_Baseline"                      => { groups => [ "Hem_Baseline",              "Hem_Serum_NoDrug" ], },
  }
};

performSmallRNA_hg19($def);

1;

