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
    task_name  => "3050-KCV-7_9to20",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160714_parclip_3050-KCV-7_mouse_mouse_9to20/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N      => 0,
    fastq_remove_random => 1,                                       #for nextflex kit
    adapter             => "TGGAATTCTCGGGTGCCAAGG",
    cqstools            => "/home/shengq1/cqstools/cqstools.exe",

    #Data
    files => {
      "UPRT_25_1"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i9_S9_R1_001.fastq.gz"],
      "UPRT_MIP_25_2_Dead" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i10_S10_R1_001.fastq.gz"],
      "UPRT_MIP_25_3_Dead" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i11_S11_R1_001.fastq.gz"],
      "UPRT_28_1"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i12_S12_R1_001.fastq.gz"],
      "UPRT_28_2_Dead"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i13_S13_R1_001.fastq.gz"],
      "UPRT_28_5"          => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i14_S14_R1_001.fastq.gz"],
      "UPRT_MIP_29_2"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i15_S15_R1_001.fastq.gz"],
      "UPRT_MIP_29_2_Dead" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i16_S16_R1_001.fastq.gz"],
      "MIP_31_1"           => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i18_S17_R1_001.fastq.gz"],
      "MIP_31_2"           => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i19_S18_R1_001.fastq.gz"],
      "MIP_31_2_Dead"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3050-KCV-7/3050/3050-KCV-7-i20_S19_R1_001.fastq.gz"],
    },
  },
  mm10_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
