#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $input = "/scratch/pingj1/chipseq/tagAlign/InputDNAAll.tagAlign.gz";

my $config = {

  #General options
  task_name            => "chipseq_pipeline",
  email                => "quanhu.sheng\@vanderbilt.edu",
  target_dir           => "/scratch/cqs/shengq1/chipseq/20160826_pingjie_chipseq_test/",
  max_thread           => 8,
  min_read_length      => 16,
  cluster              => "slurm",
  search_not_identical => 0,
  blast_unmapped_reads => 1,
  cqstools             => "/home/shengq1/cqstools/cqstools.exe",

  #Data
  files => {
    "3364-WPT-14" => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-14_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-14_2.fq.gz" ],
    "3364-WPT-15" => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-15_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-15_2.fq.gz" ],
    "3364-WPT-16" => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-16_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-16_2.fq.gz" ],
    "3364-WPT-17" => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-17_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-17_2.fq.gz" ],
    "3364-WPT-18" => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-18_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-18_2.fq.gz" ],
    "3364-WPT-19" => [ "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-19_1.fq.gz", "/workspace/pingj1/chipseq/rawdata/REPLICATE_1/3364-WPT-19_2.fq.gz" ],
  },
  groups => {
    "3364-WPT-14" => ["3364-WPT-14"],
    "3364-WPT-15" => ["3364-WPT-15"],
    "3364-WPT-16" => ["3364-WPT-16"],
    "3364-WPT-17" => ["3364-WPT-17"],
    "3364-WPT-18" => ["3364-WPT-18"],
    "3364-WPT-19" => ["3364-WPT-19"],
  },
  inputs => {
    "3364-WPT-14" => $input,
    "3364-WPT-15" => $input,
    "3364-WPT-16" => $input,
    "3364-WPT-17" => $input,
    "3364-WPT-18" => $input,
    "3364-WPT-19" => $input,
  }
};

performChIPSeq($config);

1;

