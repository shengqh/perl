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
    task_name  => "parclip",
    email      => "quanhu.sheng\@vanderbilt.edu",
    target_dir => "/scratch/cqs/shengq1/vickers/20151016_parclip_3018-KCV-44/",
    max_thread => 8,
    cluster    => "slurm",

    #Default software parameter (don't change it except you really know it)
    fastq_remove_N => 0,
    adapter        => "TGGAATTCTCGGGTGCCAAGG",

    cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

    search_3utr   => 1,
    bowtie1_index => hg19_genome()->{bowtie1_index},
    utr3_db       => "/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt",
    fasta_file    => "/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa",
    refgene_file  => "/gpfs21/scratch/cqs/shengq1/references/hg19_16569_M/hg19_refgene.tsv",

    #Data
    files => {
      "3018-KCV-44-1" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-44/3018-KCV-44-1_S1_R1_001.fastq.gz"],
      "3018-KCV-44-2" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-44/3018-KCV-44-2_S2_R1_001.fastq.gz"],
      "3018-KCV-44-3" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-44/3018-KCV-44-3_S3_R1_001.fastq.gz"],
      "3018-KCV-44-4" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-44/3018-KCV-44-4_S4_R1_001.fastq.gz"],
      "3018-KCV-44-6" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/3018-KCV-44/3018-KCV-44-6_S5_R1_001.fastq.gz"],
    },
  },
  hg19_3utr()
);

my $def = getSmallRNADefinition( $userdef, mm10_genome() );

#print Dumper($def);

performParclipSmallRNA($def);

1;

