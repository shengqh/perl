#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use Pipeline::SmallRNAUtils;
use Pipeline::ParclipSmallRNA;
use CQS::PerformSmallRNA;
use Data::Dumper;

my $userdef = {

  #General options
  task_name  => "parclip_NIH",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150925_parclip_3018-KCV-15/",
  max_thread => 8,
  cluster    => "slurm",

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,
  adapter        => "TGGAATTCTCGGGTGCCAAGG",

  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  search_3utr => 1,
  utr3_db     => "/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt",

  #Data
  files => {
    "3018-KCV-15-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-15_ATGTCA_L006_R1_001.fastq.gz"],
    "3018-KCV-15-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-36_CCAACA_L006_R1_001.fastq.gz"],
    "3018-KCV-15-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-37_CGGAAT_L006_R1_001.fastq.gz"],
    "3018-KCV-15-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-46_TCCCGA_L006_R1_001.fastq.gz"],
    "3018-KCV-15-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-15_parclip/3018-KCV-15-47_TCGAAG_L006_R1_001.fastq.gz"],
  },
};

my $def = getSmallRNADefinition( $userdef, hg19_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;

