#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;
use CQS::ClassFactory;

my $def_human = {

  #General options
  task_name            => "3116_human",
  email                => "quanhu.sheng\@vanderbilt.edu",
  target_dir           => "/scratch/cqs/shengq1/smallRNA/20160215_guoyan_smallRNA_3116/",
  max_thread           => 8,
  min_read_length      => 16,
  cluster              => "slurm",
  search_not_identical => 1,
  blast_unmapped_reads => 1,
  cqstools             => "/home/shengq1/cqstools/CQS.Tools.exe",

  #Data
  files => {
    "3116-DCs-1" => ["/scratch/cqs/guoy1/3116/3116-DCs-1_1_sequence.txt.gz"],
    "3116-DCs-2" => ["/scratch/cqs/guoy1/3116/3116-DCs-2_1_sequence.txt.gz"],
    "3116-DCs-3" => ["/scratch/cqs/guoy1/3116/3116-DCs-3_1_sequence.txt.gz"],
    "3116-DCs-4" => ["/scratch/cqs/guoy1/3116/3116-DCs-4_1_sequence.txt.gz"],
  },
  groups => {
    "TEST1" => [ "3116-DCs-1", "3116-DCs-2" ],
    "TEST2" => [ "3116-DCs-3", "3116-DCs-4" ],
  },
  pairs => {
    "TEST" => [ "TEST1", "TEST2" ],
  },
};

my $config = performSmallRNA_hg19($def_human, 0);
performTask($config, "bowtie1_unmapped_sequence_count_table");
performTask($config, "bowtie1_unmapped_sequence_blast");
performTask($config, "sequencetask");

1;

