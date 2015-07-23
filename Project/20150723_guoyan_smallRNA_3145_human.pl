#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "20150723_guoyan_smallRNA_3145_human",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/20150723_guoyan_smallRNA_3145_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "3145-YG-01" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-1_1_sequence.txt.gz"],
    "3145-YG-02" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-2_1_sequence.txt.gz"],
    "3145-YG-03" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-3_1_sequence.txt.gz"],
    "3145-YG-04" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-4_1_sequence.txt.gz"],
    "3145-YG-05" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-5_1_sequence.txt.gz"],
    "3145-YG-06" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-6_1_sequence.txt.gz"],
    "3145-YG-07" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-7_1_sequence.txt.gz"],
    "3145-YG-08" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-8_1_sequence.txt.gz"],
    "3145-YG-09" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-9_1_sequence.txt.gz"],
    "3145-YG-10" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-10_1_sequence.txt.gz"],
    "3145-YG-11" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-11_1_sequence.txt.gz"],
    "3145-YG-12" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-12_1_sequence.txt.gz"],
    "3145-YG-13" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-13_1_sequence.txt.gz"],
    "3145-YG-14" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-14_1_sequence.txt.gz"],
    "3145-YG-15" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-15_1_sequence.txt.gz"],
    "3145-YG-16" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-16_1_sequence.txt.gz"],
    "3145-YG-17" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-17_1_sequence.txt.gz"],
    "3145-YG-18" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-18_1_sequence.txt.gz"],
    "3145-YG-19" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-19_1_sequence.txt.gz"],
    "3145-YG-20" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-20_1_sequence.txt.gz"],
    "3145-YG-21" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-21_1_sequence.txt.gz"],
    "3145-YG-22" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-22_1_sequence.txt.gz"],
    "3145-YG-23" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-23_1_sequence.txt.gz"],
    "3145-YG-24" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-24_1_sequence.txt.gz"],
    "3145-YG-25" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-25_1_sequence.txt.gz"],
    "3145-YG-26" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-26_1_sequence.txt.gz"],
  },
};

performSmallRNA_hg19($def);

1;

