#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def_human = {

  #General options
  task_name       => "human_cellline_3145",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/smallRNA/20150615_guoyan_smallRNA_3145_human_cellline/",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",
  fastq_remove_N  => 0,
  run_cutadapt    => 1,

  #Data
  files => {
    "P3145-YG-1" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-1_1_sequence.txt.gz"],
    "P3145-YG-2" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-2_1_sequence.txt.gz"],
    "P3145-YG-3" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-3_1_sequence.txt.gz"],
    "P3145-YG-4" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-4_1_sequence.txt.gz"],
    "P3145-YG-5" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-5_1_sequence.txt.gz"],
    "P3145-YG-6" => ["/gpfs21/scratch/cqs/guoy1/3145/3145-YG-6_1_sequence.txt.gz"],
  }
};

performSmallRNA_hg19($def_human);
1;
