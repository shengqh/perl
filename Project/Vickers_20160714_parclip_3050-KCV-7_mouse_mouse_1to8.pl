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
    task_name  => "3050-KCV-7_1to8",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160714_parclip_3050-KCV-7_mouse_mouse_1to8/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N      => 0,
    fastq_remove_random => 1,                                       #for nextflex kit
    adapter             => "TGGAATTCTCGGGTGCCAAGG",
    cqstools            => "/home/shengq1/cqstools/cqstools.exe",

    #Data
    files => {
      "UPRT_MIP_09_1"  => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i1_S1_R1_001.fastq.gz"],
      "UPRT_MIP_10_1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i2_S2_R1_001.fastq.gz"],
      "UPRT_MIP_10_2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i3_S3_R1_001.fastq.gz"],
      "MIP_10_3"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i4_S4_R1_001.fastq.gz"],
      "UPRT_23_1"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i5_S5_R1_001.fastq.gz"],
      "UPRT_23_2"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i6_S6_R1_001.fastq.gz"],
      "UPRT_25_1"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i7_S7_R1_001.fastq.gz"],
      "UPRT_25_2"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i8_S8_R1_001.fastq.gz"],
    },
  },
  mm10_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
