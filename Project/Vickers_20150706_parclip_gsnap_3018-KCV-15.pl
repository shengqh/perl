#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "3018-KCV-15",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/vickers/20150706_parclip_gsnap_3018-KCV-15/",
  max_thread => 8,
  cluster    => "slurm",

  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  gsnap_index_directory => "/scratch/cqs/shengq1/references/hg19_16569_MT/gsnap_index_k14_2015-06-23/",
  gsnap_index_name      => "hg19_16569_MT",

  #Data
  files => {
    "3018-KCV-15-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-15_clipped_identical_NTA.fastq.gz"],
    "3018-KCV-15-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-36_clipped_identical_NTA.fastq.gz"],
    "3018-KCV-15-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-37_clipped_identical_NTA.fastq.gz"],
    "3018-KCV-15-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-46_clipped_identical_NTA.fastq.gz"],
    "3018-KCV-15-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-47_clipped_identical_NTA.fastq.gz"],
  },

  counts => {
    "3018-KCV-15-15" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-15_clipped_identical_NTA.fastq.gz.dupcount"],
    "3018-KCV-15-36" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-36_clipped_identical_NTA.fastq.gz.dupcount"],
    "3018-KCV-15-37" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-37_clipped_identical_NTA.fastq.gz.dupcount"],
    "3018-KCV-15-46" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-46_clipped_identical_NTA.fastq.gz.dupcount"],
    "3018-KCV-15-47" => ["/gpfs21/scratch/cqs/shengq1/vickers/20150312_parclip_3018-KCV-15/identical_NTA/result/3018-KCV-15-47_clipped_identical_NTA.fastq.gz.dupcount"],
  },
};

my $config = {
  general => { "task_name" => "3018-KCV-15", },
  gsnap   => {
    class                 => "Alignment::Gsnap",
    perform               => 0,
    target_dir            => $def->{target_dir} . "/gsnap",
    option                => "-y 0 -z 0 -Y 0 -Z 0 -m 1 -Q --trim-mismatch-score 0 --trim-indel-score 0 --mode ttoc-nonstranded --gunzip",
    gsnap_index_directory => $def->{gsnap_index_directory},
    gsnap_index_name      => $def->{gsnap_index_name},
    source                => $def->{files},
    sh_direct             => 1,
    pbs                   => {
      "email"    => $def->{email},
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    }
  },
  'gsnap_smallRNA_count' => {
    'cluster'         => 'slurm',
    'fasta_file'      => '/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed.fa',
    'sh_direct'       => 1,
    'perform'         => 1,
    'target_dir'      => $def->{target_dir} . "/gsnap_smallRNA_count",
    'fastq_files'     => $def->{files},
    'coordinate_file' => '/scratch/cqs/shengq1/references/smallrna/hg19_miRBase20_ucsc-tRNA_ensembl75.bed',
    'source_ref'      => 'gsnap',
    'cqs_tools'       => $def->{cqstools},
    'seqcount'        => $def->{counts},
    'samtools'        => '/scratch/cqs/shengq1/local/bin/samtools',
    'class'           => 'CQS::SmallRNACount',
    'option'          => '-s -e 4',
    'pbs'             => {
      'email'    => 'quanhu.sheng@vanderbilt.edu',
      'walltime' => '72',
      'mem'      => '40gb',
      'nodes'    => '1:ppn=1'
    },
  },
};

performConfig($config);

1;

