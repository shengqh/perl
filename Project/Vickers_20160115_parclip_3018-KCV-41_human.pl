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
    task_name  => "3018-KCV-15",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160115_parclip_3018-KCV-41_human/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    adapter        => "TGGAATTCTCGGGTGCCAAGG",
    cqstools       => "/home/shengq1/cqstools/CQS.Tools.exe",

    #Data
    files => {
      "S6_T"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i11_S1_R1_001.fastq.gz"],
      "S6_B"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i12_S2_R1_001.fastq.gz"],
      "S6_M"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i13_S3_R1_001.fastq.gz"],
      "S7_T"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i21_S4_R1_001.fastq.gz"],
      "S7_B"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i22_S5_R1_001.fastq.gz"],
      "S7_M"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i23_S6_R1_001.fastq.gz"],
      "S7_T_no_4SU"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i24_S7_R1_001.fastq.gz"],
      "S7_B_no_4SU"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i25_S8_R1_001.fastq.gz"],
      "S7_M_no_4SU"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i26_S9_R1_001.fastq.gz"],
      "S8_T" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i31_S10_R1_001.fastq.gz"],
      "S8_B" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i32_S11_R1_001.fastq.gz"],
      "S8_M" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i33_S12_R1_001.fastq.gz"],
      "S8_T_no_4SU" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i34_S13_R1_001.fastq.gz"],
      "S8_B_no_4SU" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i35_S14_R1_001.fastq.gz"],
      "S8_M_no_4SU" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-51/3018/3018-KCV-51-i36_S15_R1_001.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, hg19_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
