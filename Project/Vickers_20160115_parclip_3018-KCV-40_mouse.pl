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
    task_name  => "2797_mouse",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20160115_parclip_3018-KCV-40_mouse/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    cqstools       => "/home/shengq1/cqstools/CQS.Tools.exe",

    #Data
    files => {
      "UPRTxCD4_LIVER" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i1_S1_R1_001.fastq.gz"],
      "UPRTxCD4_WAT"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i2_S2_R1_001.fastq.gz"],
      "UPRTxCD4_ECs"   => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i3_S3_R1_001.fastq.gz"],
      "UPRT_LIVER"     => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i4_S4_R1_001.fastq.gz"],
      "UPRT_WAT"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i5_S5_R1_001.fastq.gz"],
      "UPRT_ECs"       => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i6_S6_R1_001.fastq.gz"],
      "CD4_LIVER"      => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i7_S7_R1_001.fastq.gz"],
      "CD4_WAT"        => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i8_S8_R1_001.fastq.gz"],
      "CD4_ECs"        => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-50/3018/3018-KCV-50-i9_S9_R1_001.fastq.gz"],
    },
  },
  mm10_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;

