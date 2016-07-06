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
    task_name  => "3018-KCV-81",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160706_parclip_3018-KCV-81_mouse/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    adapter        => "TGGAATTCTCGGGTGCCAAGG",
    cqstools       => "/home/shengq1/cqstools/cqstools.exe",

    #Data
    files => {
      "S01_macrophages" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i31_S1_R1_001.fastq.gz"],
      "S02_macrophages" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i32_S2_R1_001.fastq.gz"],
      "S03_macrophages" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i33_S3_R1_001.fastq.gz"],
      "S04_macrophages" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i34_S4_R1_001.fastq.gz"],
      "S05_endothelial" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i35_S5_R1_001.fastq.gz"],
      "S06_endothelial" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i36_S6_R1_001.fastq.gz"],
      "S07_endothelial" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i37_S7_R1_001.fastq.gz"],
      "S08_endothelial" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i38_S8_R1_001.fastq.gz"],
      "S09_podocytes"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i39_S9_R1_001.fastq.gz"],
      "S10_podocytes"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i40_S10_R1_001.fastq.gz"],
    },
  },
  mm10_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
