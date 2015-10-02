#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "2687",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/20151002_guoyan_smallRNA_2687_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,
  
  #Tools
  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  #Data
  files => {
    "2687-GG-170" => ["/gpfs21/scratch/cqs/guoy1/GuGuoQiang/2687-GG-170_1_sequence.txt.gz"],
    "2687-GG-172" => ["/gpfs21/scratch/cqs/guoy1/GuGuoQiang/2687-GG-172_1_sequence.txt.gz"],
    "2687-GG-173" => ["/gpfs21/scratch/cqs/guoy1/GuGuoQiang/2687-GG-173_1_sequence.txt.gz"],
    "2687-GG-176" => ["/gpfs21/scratch/cqs/guoy1/GuGuoQiang/2687-GG-176_1_sequence.txt.gz"],
    "2687-GG-178" => ["/gpfs21/scratch/cqs/guoy1/GuGuoQiang/2687-GG-178_1_sequence.txt.gz"],
  },
};

performSmallRNA_hg19($def);

1;

