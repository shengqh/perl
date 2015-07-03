#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "20150703_guoyan_smallRNA_2868_human",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/20150703_guoyan_smallRNA_2868_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 1,

  #Data
  files => {
    "2868-CMS-01" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-1_1_sequence.txt.gz"],
    "2868-CMS-02" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-2_1_sequence.txt.gz"],
    "2868-CMS-03" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-3_1_sequence.txt.gz"],
    "2868-CMS-04" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-4_1_sequence.txt.gz"],
    "2868-CMS-05" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-5_1_sequence.txt.gz"],
    "2868-CMS-06" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-6_1_sequence.txt.gz"],
    "2868-CMS-07" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-7_1_sequence.txt.gz"],
    "2868-CMS-08" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-8_1_sequence.txt.gz"],
    "2868-CMS-09" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-9_1_sequence.txt.gz"],
    "2868-CMS-10" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-10_1_sequence.txt.gz"],
    "2868-CMS-25" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-25_1_sequence.txt.gz"],
    "2868-CMS-26" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-26_1_sequence.txt.gz"],
    "2868-CMS-27" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-27_1_sequence.txt.gz"],
    "2868-CMS-28" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-28_1_sequence.txt.gz"],
    "2868-CMS-29" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-29_1_sequence.txt.gz"],
    "2868-CMS-30" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-30_1_sequence.txt.gz"],
    "2868-CMS-31" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-31_1_sequence.txt.gz"],
    "2868-CMS-32" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-32_1_sequence.txt.gz"],
    "2868-CMS-33" => ["/gpfs21/scratch/cqs/guoy1/2868/2868-CMS-33_1_sequence.txt.gz"],
  }
};

performSmallRNA_hg19($def);

1;

