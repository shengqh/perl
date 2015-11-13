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
    task_name  => "2797_rat",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20150911_parclip_gsnap_2797-rat/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    adapter        => "TGGAATTCTCGGGTGCCAAGG",
    cqstools       => "/home/shengq1/cqstools/CQS.Tools.exe",

    #Data
    files => {
      "RPI40" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI40_Ago2INS1Huh7.fastq.gz"],
      "RPI41" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI41_Ago3INS1Huh7.fastq.gz"],
      "RPI42" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI42_Ago2INS1HCEAC.fastq.gz"],
      "RPI43" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI43_Ago3INS1HCEAC.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, rn5_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;

