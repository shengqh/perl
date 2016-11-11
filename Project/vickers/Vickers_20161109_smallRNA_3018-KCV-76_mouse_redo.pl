#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA2;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name             => "3018-KCV_76_mouse",
  email                 => "quanhu.sheng\@vanderbilt.edu",
  target_dir            => "/scratch/cqs/shengq1/vickers/20161109_smallRNA_3018-KCV-76_mouse_redo",
  max_thread            => 8,
  cqstools              => "/home/shengq1/cqstools/CQS.Tools.exe",
  sequencetask_run_time => 12,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 1,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 1,
  fastq_remove_random   => 4,

  #Data
  files => {
    "IP_WT_Chow_01"     => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i25_S1_R1_001.fastq.gz' ],
    "IP_WT_Chow_02"     => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i26_S2_R1_001.fastq.gz' ],
    "IP_WT_Chow_03"     => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i27_S3_R1_001.fastq.gz' ],
    "IP_WT_Chow_04"     => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i28_S4_R1_001.fastq.gz' ],
    "IP_ApoE_Chow_01"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i29_S5_R1_001.fastq.gz' ],
    "IP_ApoE_Chow_02"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i30_S6_R1_001.fastq.gz' ],
    "IP_ApoE_Chow_03"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i31_S7_R1_001.fastq.gz' ],
    "IP_ApoE_Chow_04"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i32_S8_R1_001.fastq.gz' ],
    "IP_ApoE_HFD_01"    => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i33_S9_R1_001.fastq.gz' ],
    "IP_ApoE_HFD_02"    => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i34_S10_R1_001.fastq.gz' ],
    "IP_ApoE_HFD_03"    => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i35_S11_R1_001.fastq.gz' ],
    "IP_ApoE_HFD_04"    => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i36_S12_R1_001.fastq.gz' ],
    "FPLC_WT_Chow_01"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i37_S13_R1_001.fastq.gz' ],
    "FPLC_WT_Chow_02"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i38_S14_R1_001.fastq.gz' ],
    "FPLC_WT_Chow_03"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i39_S15_R1_001.fastq.gz' ],
    "FPLC_WT_Chow_04"   => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i40_S16_R1_001.fastq.gz' ],
    "FPLC_ApoE_Chow_01" => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i41_S17_R1_001.fastq.gz' ],
    "FPLC_ApoE_Chow_02" => [
      '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i42_S18_R1_001.fastq.gz'
    ],
    "FPLC_ApoE_Chow_03" => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i43_S19_R1_001.fastq.gz' ],
    "FPLC_ApoE_Chow_04" => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i44_S20_R1_001.fastq.gz' ],
    "FPLC_ApoE_HFD_01"  => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i45_S21_R1_001.fastq.gz' ],
    "FPLC_ApoE_HFD_03"  => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i46_S22_R1_001.fastq.gz' ],
    "FPLC_ApoE_HFD_04"  => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i47_S23_R1_001.fastq.gz' ],
    "H2O"               => [ '/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-76/3018/3018-KCV-76-i48_S24_R1_001.fastq.gz' ],

  },
  groups => {
    "FPLC_ApoE_Chow" => [ "FPLC_ApoE_Chow_01", "FPLC_ApoE_Chow_02", "FPLC_ApoE_Chow_03", "FPLC_ApoE_Chow_04" ],
    "FPLC_ApoE_HFD"  => [ "FPLC_ApoE_HFD_01",  "FPLC_ApoE_HFD_03",  "FPLC_ApoE_HFD_04" ],
    "FPLC_WT_Chow"   => [ "FPLC_WT_Chow_01",   "FPLC_WT_Chow_02",   "FPLC_WT_Chow_03",   "FPLC_WT_Chow_04" ],
    "IP_ApoE_Chow"   => [ "IP_ApoE_Chow_01",   "IP_ApoE_Chow_02",   "IP_ApoE_Chow_03",   "IP_ApoE_Chow_04" ],
    "IP_ApoE_HFD"    => [ "IP_ApoE_HFD_01",    "IP_ApoE_HFD_02",    "IP_ApoE_HFD_03",    "IP_ApoE_HFD_04" ],
    "IP_WT_Chow"     => [ "IP_WT_Chow_01",     "IP_WT_Chow_02",     "IP_WT_Chow_03",     "IP_WT_Chow_04" ],
  },
  pairs => {
    "IP_ApoE_Chow_VS_IP_WT_Chow"      => { groups => [ "IP_WT_Chow",     "IP_ApoE_Chow" ], },
    "IP_ApoE_HFD_VS_IP_ApoE_Chow"     => { groups => [ "IP_ApoE_Chow",   "IP_ApoE_HFD" ], },
    "IP_ApoE_HFD_VS_IP_WT_Chow"       => { groups => [ "IP_WT_Chow",     "IP_ApoE_HFD" ], },
    "FPLC_ApoE_Chow_VS_FPLC_WT_Chow"  => { groups => [ "FPLC_WT_Chow",   "FPLC_ApoE_Chow" ], },
    "FPLC_ApoE_HFD_VS_FPLC_ApoE_Chow" => { groups => [ "FPLC_ApoE_Chow", "FPLC_ApoE_HFD" ], },
    "FPLC_ApoE_HFD_VS_FPLC_WT_Chow"   => { groups => [ "FPLC_WT_Chow",   "FPLC_ApoE_HFD" ], },

  }
};

performSmallRNA_mm10($def);
1;

