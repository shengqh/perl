#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "20150723_guoyan_smallRNA_2687_human",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/20150723_guoyan_smallRNA_2687_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "2687-GG-129" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-129_1_sequence.txt.gz"],
    "2687-GG-138" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-138_1_sequence.txt.gz"],
    "2687-GG-142" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-142_1_sequence.txt.gz"],
    "2687-GG-149" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-149_1_sequence.txt.gz"],
    "2687-GG-152" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-152_1_sequence.txt.gz"],
    "2687-GG-155" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-155_1_sequence.txt.gz"],
    "2687-GG-158" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-158_1_sequence.txt.gz"],
    "2687-GG-159" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-159_1_sequence.txt.gz"],
    "2687-GG-161" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-161_1_sequence.txt.gz"],
    "2687-GG-162" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-162_1_sequence.txt.gz"],
    "2687-GG-163" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-163_1_sequence.txt.gz"],
    "2687-GG-164" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-164_1_sequence.txt.gz"],
    "2687-GG-165" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-165_1_sequence.txt.gz"],
    "2687-GG-166" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-166_1_sequence.txt.gz"],
    "2687-GG-167" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-167_1_sequence.txt.gz"],
    "2687-GG-168" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-168_1_sequence.txt.gz"],
    "2687-GG-169" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-169_1_sequence.txt.gz"],
  },
};

performSmallRNA_hg19($def);

1;

