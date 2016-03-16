#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $def = {

  #General options
  task_name  => "smallrna",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/20151016_guoyan_smallRNA_2501_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N        => 0,
  cqstools              => "/home/shengq1/cqstools/CQS.Tools.exe",
  search_not_identical  => 0,
  search_unmapped_reads => 0,
  blast_unmapped_reads  => 0,
  #Data
  files => {
    "2501-ASK-01" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-1_1.fastq.gz"],
    "2501-ASK-02" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-2_1.fastq.gz"],
    "2501-ASK-03" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-3ab_1.fastq.gz"],
    "2501-ASK-04" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-4_1.fastq.gz"],
    "2501-ASK-05" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-5_1.fastq.gz"],
    "2501-ASK-06" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-6_1.fastq.gz"],
    "2501-ASK-07" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-7_1.fastq.gz"],
    "2501-ASK-08" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-8_1.fastq.gz"],
    "2501-ASK-09" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-9ab_1.fastq.gz"],
    "2501-ASK-10" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-10_1.fastq.gz"],
    "2501-ASK-11" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-11_1.fastq.gz"],
    "2501-ASK-12" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-12_1.fastq.gz"],
    "2501-ASK-13" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-13_1.fastq.gz"],
    "2501-ASK-14" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-14_1.fastq.gz"],
    "2501-ASK-15" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-15_1.fastq.gz"],
    "2501-ASK-16" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-16_1.fastq.gz"],
    "2501-ASK-17" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-17_1.fastq.gz"],
    "2501-ASK-18" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-18_1.fastq.gz"],
    "2501-ASK-19" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-19ab_1.fastq.gz"],
    "2501-ASK-20" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-20_1.fastq.gz"],
    "2501-ASK-21" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-21_1.fastq.gz"],
    "2501-ASK-22" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-22_1.fastq.gz"],
    "2501-ASK-23" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-23_1.fastq.gz"],
    "2501-ASK-24" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-24_1.fastq.gz"],
    "2501-ASK-25" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-25_1.fastq.gz"],
    "2501-ASK-27" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-27_1.fastq.gz"],
    "2501-ASK-28" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-28ab_1.fastq.gz"],
    "2501-ASK-30" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-30_1.fastq.gz"],
    "2501-ASK-32" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-32ab_1.fastq.gz"],
    "2501-ASK-33" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-33_1.fastq.gz"],
    "2501-ASK-34" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-34ab_1.fastq.gz"],
    "2501-ASK-35" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-35ab_1.fastq.gz"],
    "2501-ASK-36" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-36ab_1.fastq.gz"],
    "2501-ASK-37" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-37_1.fastq.gz"],
    "2501-ASK-38" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-38_1.fastq.gz"],
    "2501-ASK-39" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-39ab_1.fastq.gz"],
    "2501-ASK-40" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-40_1.fastq.gz"],
    "2501-ASK-41" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-41_1.fastq.gz"],
    "2501-ASK-42" => ["/data/cqs/guom1/2501_rawdata_1/2501-ASK-42_1.fastq.gz"],
    "2501-ASK-44" => ["/data/cqs/guom1/2501_rawdata_2/2501-ASK-44ab_1.fastq.gz"],
  },
  tRNA_vis_group => {

    #"AML" => [
    "C" => [ "2501-ASK-01", "2501-ASK-04", "2501-ASK-15", "2501-ASK-17", "2501-ASK-21", "2501-ASK-23", "2501-ASK-24", "2501-ASK-25", "2501-ASK-27", "2501-ASK-37", "2501-ASK-40" ],

    #    "MDS" => [
    "A" => [
      "2501-ASK-02", "2501-ASK-03", "2501-ASK-05", "2501-ASK-08", "2501-ASK-10", "2501-ASK-11", "2501-ASK-12", "2501-ASK-13", "2501-ASK-14", "2501-ASK-18",
      "2501-ASK-28", "2501-ASK-30", "2501-ASK-32", "2501-ASK-33", "2501-ASK-35", "2501-ASK-41", "2501-ASK-42", "2501-ASK-44"
    ],

    #"MDS-AML" => [
    "B" => [ "2501-ASK-06", "2501-ASK-07", "2501-ASK-09", "2501-ASK-16", "2501-ASK-19", "2501-ASK-20", "2501-ASK-22", "2501-ASK-34", "2501-ASK-36", "2501-ASK-38", "2501-ASK-39" ]
  }
};

my $config = performSmallRNA_hg19( $def, 1 );

#performTask( $config, "tRNA_PositionVis" );

1;

