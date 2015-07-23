#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def = {

  #General options
  task_name  => "20150723_guoyan_smallRNA_3166_human",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => "/scratch/cqs/shengq1/smallRNA/20150723_guoyan_smallRNA_3166_human",
  max_thread => 8,

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,

  #Data
  files => {
    "3166-CRF-1" => ["/gpfs21/scratch/cqs/guoy1/3145/3166-CRF-1_1_sequence.txt.gz"],
    "3166-CRF-2" => ["/gpfs21/scratch/cqs/guoy1/3145/3166-CRF-2_1_sequence.txt.gz"],
  },
};

performSmallRNA_hg19($def);

1;

