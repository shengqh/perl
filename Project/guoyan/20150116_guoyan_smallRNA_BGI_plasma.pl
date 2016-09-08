#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def = {

  #General options
  task_name       => "BGI_Plasma",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/smallRNA/20150116_guoyan_smallRNA_BGI_plasma/",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",

  #Data
  files => {
  "BGI_Plasma-40" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/110831_I189_FCC01EUABXX_L8_WATiveSATTAASE-X-40_1.fq.gz"],
  "BGI_Plasma-41" => ["/gpfs21/scratch/cqs/guoy1/BGI_miRNA/110831_I189_FCC01EUABXX_L8_WATiveSAUTAASE-X-41_1.fq.gz"],
  }
};

performSmallRNAHuman($def);

1;

