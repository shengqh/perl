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
    task_name  => "NIH",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20151113_parclip_NIH/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    adapter        => "TGGAATTCTCGGGTGCCAAGG",

    cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

    #Data
    files => {
      "Parclip_01" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_1_ATCACG_L002_R1.fastq.gz"],
      "Parclip_02" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_2_CGATGT_L002_R1.fastq.gz"],
      "Parclip_03" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_3_TTAGGC_L002_R1.fastq.gz"],
      "Parclip_04" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_4_TGACCA_L002_R1.fastq.gz"],
      "Parclip_05" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_5_ACAGTG_L002_R1.fastq.gz"],
      "Parclip_06" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_6_GCCAAT_L002_R1.fastq.gz"],
      "Parclip_07" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_7_CAGATC_L002_R1.fastq.gz"],
      "Parclip_08" => ["/scratch/cqs/shengq1/vickers/data/201312_parclip_NIH/Vickers_Parclip_8_ACTTGA_L002_R1.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;

