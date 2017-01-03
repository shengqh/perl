#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA3;
use CQS::ClassFactory;

my $email = "quanhu.sheng\@vanderbilt.edu";

my $def = {

  #General options
  task_name                 => "KCV_3018_77_78_79_86",
  email                     => $email,
  target_dir                => "/scratch/cqs/shengq1/vickers/20161223_smallRNA_3018-KCV-77_78_79_86_v3",
  max_thread                => 8,
  cqstools                  => "/home/shengq1/cqstools/cqstools.exe",
  sequencetask_run_time     => 6,
  table_vis_group_text_size => 12,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 1,
  remove_sequences      => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
  search_unmapped_reads => 1,
  blast_unmapped_reads  => 0,
  top_read_number       => 100,
  blast_top_reads       => 0,
  blast_localdb         => "/scratch/cqs/shengq1/references/blastdb",

  #next flex
  fastq_remove_random => 4,

  #Data
  files => {

    "APOB_WT_01"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i42_S19_R1_001.fastq.gz'],
    "APOB_WT_03"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i44_S21_R1_001.fastq.gz'],
    "APOB_WT_05"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i46_S20_R1_001.fastq.gz'],
    "APOB_WT_08"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i1_S1_R1_001.fastq.gz'],
    "APOB_WT_10"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i3_S1_R1_001.fastq.gz'],
    "APOB_WT_12"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i5_S3_R1_001.fastq.gz'],
    "APOB_WT_14"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i7_S5_R1_001.fastq.gz'],
    "APOB_WT_86_04"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i4_S4_R1_001.fastq.gz'],
    "APOB_WT_86_05"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i5_S5_R1_001.fastq.gz'],
    "APOB_WT_86_06"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i6_S6_R1_001.fastq.gz'],
    "APOB_SRBIKO_02"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i43_S20_R1_001.fastq.gz'],
    "APOB_SRBIKO_04"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i45_S22_R1_001.fastq.gz'],
    "APOB_SRBIKO_06"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i47_S21_R1_001.fastq.gz'],
    "APOB_SRBIKO_07"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i48_S22_R1_001.fastq.gz'],
    "APOB_SRBIKO_09"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i2_S2_R1_001.fastq.gz'],
    "APOB_SRBIKO_11"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i4_S2_R1_001.fastq.gz'],
    "APOB_SRBIKO_13"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i6_S4_R1_001.fastq.gz'],
    "APOB_SRBIKO_86_01" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i1_S1_R1_001.fastq.gz'],
    "APOB_SRBIKO_86_03" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i3_S3_R1_001.fastq.gz'],
    "APOB_SRBIKO_86_20" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i2_S2_R1_001.fastq.gz'],

    "Bile_WT_01"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i15_S9_R1_001.fastq.gz'],
    "Bile_WT_03"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i17_S11_R1_001.fastq.gz'],
    "Bile_WT_05"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i19_S13_R1_001.fastq.gz'],
    "Bile_WT_08"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i22_S14_R1_001.fastq.gz'],
    "Bile_WT_10"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i24_S14_R1_001.fastq.gz'],
    "Bile_WT_12"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i25_S15_R1_001.fastq.gz'],
    "Bile_WT_14"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i27_S17_R1_001.fastq.gz'],
    "Bile_WT_86_04"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i16_S16_R1_001.fastq.gz'],
    "Bile_WT_86_05"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i17_S17_R1_001.fastq.gz'],
    "Bile_WT_86_06"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i18_S18_R1_001.fastq.gz'],
    "Bile_SRBIKO_02"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i16_S10_R1_001.fastq.gz'],
    "Bile_SRBIKO_04"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i18_S12_R1_001.fastq.gz'],
    "Bile_SRBIKO_06"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i20_S12_R1_001.fastq.gz'],
    "Bile_SRBIKO_07"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i21_S13_R1_001.fastq.gz'],
    "Bile_SRBIKO_09"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i23_S15_R1_001.fastq.gz'],
    "Bile_SRBIKO_13"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i26_S16_R1_001.fastq.gz'],
    "Bile_SRBIKO_86_01" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i13_S13_R1_001.fastq.gz'],
    "Bile_SRBIKO_86_03" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i15_S15_R1_001.fastq.gz'],
    "Bile_SRBIKO_86_20" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i14_S14_R1_001.fastq.gz'],

    "HDL_WT_01"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i28_S14_R1_001.fastq.gz'],
    "HDL_WT_03"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i30_S16_R1_001.fastq.gz'],
    "HDL_WT_08"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i33_S16_R1_001.fastq.gz'],
    "HDL_WT_10"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i35_S18_R1_001.fastq.gz'],
    "HDL_WT_05"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i38_S19_R1_001.fastq.gz'],
    "HDL_WT_12"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i39_S20_R1_001.fastq.gz'],
    "HDL_WT_14"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i41_S22_R1_001.fastq.gz'],
    "HDL_WT_86_04"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i10_S10_R1_001.fastq.gz'],
    "HDL_WT_86_05"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i11_S11_R1_001.fastq.gz'],
    "HDL_WT_86_06"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i12_S12_R1_001.fastq.gz'],
    "HDL_SRBIKO_02"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i29_S15_R1_001.fastq.gz'],
    "HDL_SRBIKO_06"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i31_S17_R1_001.fastq.gz'],
    "HDL_SRBIKO_07"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i32_S18_R1_001.fastq.gz'],
    "HDL_SRBIKO_09"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i34_S17_R1_001.fastq.gz'],
    "HDL_SRBIKO_11"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i36_S19_R1_001.fastq.gz'],
    "HDL_SRBIKO_04"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i37_S18_R1_001.fastq.gz'],
    "HDL_SRBIKO_13"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i40_S21_R1_001.fastq.gz'],
    "HDL_SRBIKO_86_01" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i7_S7_R1_001.fastq.gz'],
    "HDL_SRBIKO_86_03" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i9_S9_R1_001.fastq.gz'],
    "HDL_SRBIKO_86_20" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i8_S8_R1_001.fastq.gz'],

    "Liver_WT_01"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i1_S1_R1_001.fastq.gz'],
    "Liver_WT_03"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i3_S3_R1_001.fastq.gz'],
    "Liver_WT_05"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i5_S5_R1_001.fastq.gz'],
    "Liver_WT_08"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i8_S5_R1_001.fastq.gz'],
    "Liver_WT_10"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i10_S6_R1_001.fastq.gz'],
    "Liver_WT_12"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i12_S8_R1_001.fastq.gz'],
    "Liver_WT_14"        => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i14_S10_R1_001.fastq.gz'],
    "Liver_WT_86_04"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i22_S22_R1_001.fastq.gz'],
    "Liver_WT_86_05"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i23_S23_R1_001.fastq.gz'],
    "Liver_WT_86_06"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i24_S24_R1_001.fastq.gz'],
    "Liver_SRBIKO_02"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i2_S2_R1_001.fastq.gz'],
    "Liver_SRBIKO_04"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i4_S4_R1_001.fastq.gz'],
    "Liver_SRBIKO_06"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i6_S3_R1_001.fastq.gz'],
    "Liver_SRBIKO_07"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i7_S4_R1_001.fastq.gz'],
    "Liver_SRBIKO_09"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i9_S6_R1_001.fastq.gz'],
    "Liver_SRBIKO_11"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i11_S7_R1_001.fastq.gz'],
    "Liver_SRBIKO_13"    => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i13_S9_R1_001.fastq.gz'],
    "Liver_SRBIKO_86_01" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i19_S19_R1_001.fastq.gz'],
    "Liver_SRBIKO_86_03" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i21_S21_R1_001.fastq.gz'],
    "Liver_SRBIKO_86_20" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-86/3018-KCV-86-i20_S20_R1_001.fastq.gz'],

    "Urine_WT_01"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i8_S6_R1_001.fastq.gz'],
    "Urine_WT_03"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i10_S8_R1_001.fastq.gz'],
    "Urine_WT_05"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i11_S7_R1_001.fastq.gz'],
    "Urine_WT_08"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i14_S10_R1_001.fastq.gz'],
    "Urine_WT_14"     => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i18_S13_R1_001.fastq.gz'],
    "Urine_SRBIKO_02" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-77/3018/3018-KCV-77-i9_S7_R1_001.fastq.gz'],
    "Urine_SRBIKO_06" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i12_S8_R1_001.fastq.gz'],
    "Urine_SRBIKO_07" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i13_S9_R1_001.fastq.gz'],
    "Urine_SRBIKO_09" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-78/3018/3018-KCV-78-i15_S11_R1_001.fastq.gz'],
    "Urine_SRBIKO_11" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i16_S11_R1_001.fastq.gz'],
    "Urine_SRBIKO_13" => ['/gpfs21/scratch/cqs/zhaos/vickers/data/3018/3018-KCV-79/3018/3018-KCV-79-i17_S12_R1_001.fastq.gz'],
  },
  groups => {
    "APOB_WT" => [ "APOB_WT_01", "APOB_WT_03", "APOB_WT_05", "APOB_WT_08", "APOB_WT_10", "APOB_WT_12", "APOB_WT_14", "APOB_WT_86_04", "APOB_WT_86_05", "APOB_WT_86_06" ],
    "APOB_SRBIKO" =>
      [ "APOB_SRBIKO_02", "APOB_SRBIKO_04", "APOB_SRBIKO_06", "APOB_SRBIKO_07", "APOB_SRBIKO_09", "APOB_SRBIKO_11", "APOB_SRBIKO_13", "APOB_SRBIKO_86_01", "APOB_SRBIKO_86_03", "APOB_SRBIKO_86_20" ],
    "Bile_WT" => [ "Bile_WT_01", "Bile_WT_03", "Bile_WT_05", "Bile_WT_08", "Bile_WT_10", "Bile_WT_12", "Bile_WT_14", "Bile_WT_86_04", "Bile_WT_86_05", "Bile_WT_86_06" ],
    "Bile_SRBIKO" => [ "Bile_SRBIKO_02", "Bile_SRBIKO_04", "Bile_SRBIKO_06", "Bile_SRBIKO_07", "Bile_SRBIKO_09", "Bile_SRBIKO_13", "Bile_SRBIKO_86_01", "Bile_SRBIKO_86_03", "Bile_SRBIKO_86_20" ],
    "HDL_WT" => [ "HDL_WT_01", "HDL_WT_03", "HDL_WT_08", "HDL_WT_10", "HDL_WT_05", "HDL_WT_12", "HDL_WT_14", "HDL_WT_86_04", "HDL_WT_86_05", "HDL_WT_86_06" ],
    "HDL_SRBIKO" =>
      [ "HDL_SRBIKO_02", "HDL_SRBIKO_06", "HDL_SRBIKO_07", "HDL_SRBIKO_09", "HDL_SRBIKO_11", "HDL_SRBIKO_04", "HDL_SRBIKO_13", "HDL_SRBIKO_86_01", "HDL_SRBIKO_86_03", "HDL_SRBIKO_86_20" ],
    "Liver_WT"     => [ "Liver_WT_01", "Liver_WT_03", "Liver_WT_05", "Liver_WT_08", "Liver_WT_10", "Liver_WT_12", "Liver_WT_14", "Liver_WT_86_04", "Liver_WT_86_05", "Liver_WT_86_06" ],
    "Liver_SRBIKO" => [
      "Liver_SRBIKO_02", "Liver_SRBIKO_04",    "Liver_SRBIKO_06",    "Liver_SRBIKO_07", "Liver_SRBIKO_09", "Liver_SRBIKO_11",
      "Liver_SRBIKO_13", "Liver_SRBIKO_86_01", "Liver_SRBIKO_86_03", "Liver_SRBIKO_86_20"
    ],
    "Urine_WT"     => [ "Urine_WT_01",     "Urine_WT_03",     "Urine_WT_05",     "Urine_WT_08",     "Urine_WT_14" ],
    "Urine_SRBIKO" => [ "Urine_SRBIKO_02", "Urine_SRBIKO_06", "Urine_SRBIKO_07", "Urine_SRBIKO_09", "Urine_SRBIKO_11", "Urine_SRBIKO_13" ],
  },
  pairs => {
    "APOB_SRBIKO_vs_WT" => {
      groups => [ "APOB_WT", "APOB_SRBIKO" ],
      Batch  => [
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86',
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86'
      ],
    },

    "Bile_SRBIKO_vs_WT" => {
      groups => [ "Bile_WT", "Bile_SRBIKO" ],
      Batch  => [
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86',
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86'
      ],
    },
    "HDL_SRBIKO_vs_WT" => {
      groups => [ "HDL_WT", "HDL_SRBIKO" ],
      Batch  => [
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86',
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86'
      ],
    },
    "Liver_SRBIKO_vs_WT" => {
      groups => [ "Liver_WT", "Liver_SRBIKO" ],
      Batch  => [
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86',
        '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79', '3018-KCV-86', '3018-KCV-86', '3018-KCV-86'
      ],
    },
    "Urine_SRBIKO_vs_WT" => {
      groups => [ "Urine_WT",    "Urine_SRBIKO" ],
      Batch  => [ '3018-KCV-77', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-77', '3018-KCV-78', '3018-KCV-78', '3018-KCV-78', '3018-KCV-79', '3018-KCV-79' ],
    },
  },
  groups_vis_layout => {
    "Row_Group" => [ "HDL",    "HDL",        "APOB",    "APOB",        "LIVER",    "LIVER",        "BILE",    "BILE",        "URINE",    "URINE" ],
    "Col_Group" => [ "WT",     "SR-BI KO",   "WT",      "SR-BI KO",    "WT",       "SR-BI KO",     "WT",      "SR-BI KO",    "WT",       "SR-BI KO" ],
    "Groups"    => [ "HDL_WT", "HDL_SRBIKO", "APOB_WT", "APOB_SRBIKO", "Liver_WT", "Liver_SRBIKO", "Bile_WT", "Bile_SRBIKO", "Urine_WT", "Urine_SRBIKO" ],
  },
  pairs_host_deseq2_vis_layout => {
    "Col_Group" => [ "HDL", "HDL", "HDL", "APOB", "APOB", "APOB", "LIVER", "LIVER", "LIVER", "BILE", "BILE", "BILE", "URINE", "URINE", "URINE" ],
    "Row_Group" =>
      [ "miRNA", "tRNA", "Other small RNA", "miRNA", "tRNA", "Other small RNA", "miRNA", "tRNA", "Other small RNA", "miRNA", "tRNA", "Other small RNA", "miRNA", "tRNA", "Other small RNA" ],
    "Groups" => [
      "deseq2_miRNA_HDL_SRBIKO_vs_WT",           "deseq2_tRNA_HDL_SRBIKO_vs_WT",           "deseq2_otherSmallRNA_HDL_SRBIKO_vs_WT", "deseq2_miRNA_APOB_SRBIKO_vs_WT",
      "deseq2_tRNA_APOB_SRBIKO_vs_WT",           "deseq2_otherSmallRNA_APOB_SRBIKO_vs_WT", "deseq2_miRNA_Liver_SRBIKO_vs_WT",       "deseq2_tRNA_Liver_SRBIKO_vs_WT",
      "deseq2_otherSmallRNA_Liver_SRBIKO_vs_WT", "deseq2_miRNA_Bile_SRBIKO_vs_WT",         "deseq2_tRNA_Bile_SRBIKO_vs_WT",         "deseq2_otherSmallRNA_Bile_SRBIKO_vs_WT",
      "deseq2_miRNA_Urine_SRBIKO_vs_WT",         "deseq2_tRNA_Urine_SRBIKO_vs_WT",         "deseq2_otherSmallRNA_Urine_SRBIKO_vs_WT"
    ],
  },
  pairs_nonHostGroups_deseq2_vis_layout => {
    "Col_Group" => [ "HDL", "HDL", "HDL", "APOB", "APOB", "APOB", "LIVER", "LIVER", "LIVER", "BILE", "BILE", "BILE", "URINE", "URINE", "URINE" ],
    "Row_Group" => [
      "Microbiome", "Environment", "Fungus",      "Microbiome", "Environment", "Fungus",      "Microbiome", "Environment",
      "Fungus",     "Microbiome",  "Environment", "Fungus",     "Microbiome",  "Environment", "Fungus"
    ],
    "Groups" => [
      "deseq2_bacteria_group1_HDL_Knockout_VS_WildType",   "deseq2_bacteria_group2_HDL_Knockout_VS_WildType",
      "deseq2_fungus_group4_HDL_Knockout_VS_WildType",     "deseq2_bacteria_group1_APOB_Knockout_VS_WildType",
      "deseq2_bacteria_group2_APOB_Knockout_VS_WildType",  "deseq2_fungus_group4_APOB_Knockout_VS_WildType",
      "deseq2_bacteria_group1_Liver_Knockout_VS_WildType", "deseq2_bacteria_group2_Liver_Knockout_VS_WildType",
      "deseq2_fungus_group4_Liver_Knockout_VS_WildType",   "deseq2_bacteria_group1_Bile_Knockout_VS_WildType",
      "deseq2_bacteria_group2_Bile_Knockout_VS_WildType",  "deseq2_fungus_group4_Bile_Knockout_VS_WildType",
      "deseq2_bacteria_group1_Urine_Knockout_VS_WildType", "deseq2_bacteria_group2_Urine_Knockout_VS_WildType",
      "deseq2_fungus_group4_Urine_Knockout_VS_WildType"
    ],
  },
  pairs_nonHosttRNArRNA_deseq2_vis_layout => {
    "Col_Group" => [ "HDL",  "HDL",   "HDL",   "APOB", "APOB",  "APOB",  "LIVER", "LIVER", "LIVER", "BILE", "BILE",  "BILE",  "URINE", "URINE", "URINE" ],
    "Row_Group" => [ "tRNA", "rRNAL", "rRNAS", "tRNA", "rRNAL", "rRNAS", "tRNA",  "rRNAL", "rRNAS", "tRNA", "rRNAL", "rRNAS", "tRNA",  "rRNAL", "rRNAS" ],
    "Groups"    => [
      "deseq2_nonHost_tRna_HDL_Knockout_VS_WildType",    "deseq2_nonHost_rRNAL_HDL_Knockout_VS_WildType",
      "deseq2_nonHost_rRNAS_HDL_Knockout_VS_WildType",   "deseq2_nonHost_tRna_APOB_Knockout_VS_WildType",
      "deseq2_nonHost_rRNAL_APOB_Knockout_VS_WildType",  "deseq2_nonHost_rRNAS_APOB_Knockout_VS_WildType",
      "deseq2_nonHost_tRna_Liver_Knockout_VS_WildType",  "deseq2_nonHost_rRNAL_Liver_Knockout_VS_WildType",
      "deseq2_nonHost_rRNAS_Liver_Knockout_VS_WildType", "deseq2_nonHost_tRna_Bile_Knockout_VS_WildType",
      "deseq2_nonHost_rRNAL_Bile_Knockout_VS_WildType",  "deseq2_nonHost_rRNAS_Bile_Knockout_VS_WildType",
      "deseq2_nonHost_tRna_Urine_Knockout_VS_WildType",  "deseq2_nonHost_rRNAL_Urine_Knockout_VS_WildType",
      "deseq2_nonHost_rRNAS_Urine_Knockout_VS_WildType"
    ],
  },
};

my $config = performSmallRNA_mm10( $def, 0 );
performTask($config, "count_table_correlation");

1;

