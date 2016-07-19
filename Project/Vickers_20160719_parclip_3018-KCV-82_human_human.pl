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
    task_name  => "3018-KCV-82-H-H",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160719_parclip_3018-KCV-82_human_human/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N      => 0,
    fastq_remove_random => 0,                                       #for trueseq kit
    adapter             => "TGGAATTCTCGGGTGCCAAGG",
    cqstools            => "/home/shengq1/cqstools/cqstools.exe",

    #Data
    files => {
      "T1_T_T"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i1_S1_R1_001.fastq.gz"],
      "B1_T_B"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i2_S2_R1_001.fastq.gz"],
      "M1_T_M"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i3_S3_R1_001.fastq.gz"],
      "T2_T_T"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i4_S4_R1_001.fastq.gz"],
      "B2_T_B"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i5_S5_R1_001.fastq.gz"],
      "M2_T_M"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i6_S6_R1_001.fastq.gz"],
      "T3_T_T"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i7_S7_R1_001.fastq.gz"],
      "B3_T_B"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i8_S8_R1_001.fastq.gz"],
      "M3_T_M"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i9_S9_R1_001.fastq.gz"],
      "KW8_EC_Podocyte" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-82/3018/3018-KCV-82-i38_S10_R1_001.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, hg19_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
