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
    task_name  => "3018-KCV-81-M-H",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160706_parclip_3018-KCV-81_mouse_human/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    adapter        => "TGGAATTCTCGGGTGCCAAGG",
    cqstools       => "/home/shengq1/cqstools/cqstools.exe",

    #Data
    files => {
      "S01_M_P" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i31_S1_R1_001.fastq.gz"],
      "S02_M_E" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i32_S2_R1_001.fastq.gz"],
      "S03_M_P" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i33_S3_R1_001.fastq.gz"],
      "S04_M_P" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-81/3018/3018-KCV-81-i34_S4_R1_001.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;
