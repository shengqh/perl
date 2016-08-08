#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use Pipeline::SmallRNAUtils;
use Pipeline::ParclipSmallRNA;
use CQS::PerformSmallRNA;
use Data::Dumper;

use Hash::Merge qw( merge );

my $userdef = merge(
  {

    #General options
    task_name  => "3050-KCV-08",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160721_parclip_3050-KCV-08_mouse_mouse/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N      => 0,
    remove_sequences    => "'CCACGTTCCCGTGG;ACAGTCCGACGATC'",
    fastq_remove_random => 4,                                       #for nextflex kit
    adapter             => "TGGAATTCTCGGGTGCCAAGG",
    cqstools            => "/home/shengq1/cqstools/cqstools.exe",

    #Data
    files => {
      "UPRT_MIP_25_2"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i26_S1_R1_001.fastq.gz"],
      "UPRT_MIP_25_3"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i27_S2_R1_001.fastq.gz"],
      "UPRT_28_1"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i28_S3_R1_001.fastq.gz"],
      "UPRT_28_2"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i29_S4_R1_001.fastq.gz"],
      "UPRT_MIP_28_4"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i30_S5_R1_001.fastq.gz"],
      "UPRT_28_5"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i31_S6_R1_001.fastq.gz"],
      "UPRT_MIP_29_2"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i32_S7_R1_001.fastq.gz"],
      "UPRT_MIP_29_3"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i33_S8_R1_001.fastq.gz"],
      "MIP_31_1"        => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i34_S9_R1_001.fastq.gz"],
      "MIP_31_2"        => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i35_S10_R1_001.fastq.gz"],
      "UPRT_MIP_31_3"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i36_S11_R1_001.fastq.gz"],
      "UPRT_MIP_25_2_X" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i38_S12_R1_001.fastq.gz"],
      "UPRT_MIP_25_3_X" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i39_S13_R1_001.fastq.gz"],
      "UPRT_28_1_X"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i40_S14_R1_001.fastq.gz"],
      "UPRT_28_2_X"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i41_S15_R1_001.fastq.gz"],
      "UPRT_MIP_28_4_X" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i42_S16_R1_001.fastq.gz"],
      "UPRT_28_5_X"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i43_S17_R1_001.fastq.gz"],
      "UPRT_MIP_29_2_X" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i44_S18_R1_001.fastq.gz"],
      "UPRT_MIP_29_3_X" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i45_S19_R1_001.fastq.gz"],
      "MIP_31_1_X"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i46_S20_R1_001.fastq.gz"],
      "MIP_31_2_X"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i47_S21_R1_001.fastq.gz"],
      "UPRT_MIP_31_3_X" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-08/3050/3050-KCV-8-i48_S22_R1_001.fastq.gz"],
    },
  },
  mm10_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
