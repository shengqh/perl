#!/usr/bin/perl
use strict;
use warnings;

use CQS::SmallRNA;

my $def_human = {

  #General options
  task_name       => "2687_human",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/smallRNA/20150206_guoyan_2687_human_mouse/human/",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",

  #Data
  files => {
    "2687-GG-39" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-39_1_sequence.txt.gz"],
    "2687-GG-40" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-40_1_sequence.txt.gz"],
    "2687-GG-41" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-41_1_sequence.txt.gz"],
    "2687-GG-61" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-61_1_sequence.txt.gz"],
    "2687-GG-65" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-65_1_sequence.txt.gz"],
    "2687-GG-87" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-87_1_sequence.txt.gz"],
    "2687-GG-88" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-88_1_sequence.txt.gz"],
    "2687-GG-90" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-90_1_sequence.txt.gz"],
    "2687-GG-91" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-91_1_sequence.txt.gz"],
    "2687-GG-92" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-92_1_sequence.txt.gz"],
    "2687-GG-94" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-94_1_sequence.txt.gz"],
  }
};

performSmallRNAHuman($def_human);

my $def_mouse = {

  #General options
  task_name       => "2687_mouse",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/smallRNA/20150206_guoyan_2687_human_mouse/mouse/",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",

  #Data
  files => {
    "2687-GG-35" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-35_1_sequence.txt.gz"],
    "2687-GG-44" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-44_1_sequence.txt.gz"],
    "2687-GG-47" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-47_1_sequence.txt.gz"],
    "2687-GG-52" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-52_1_sequence.txt.gz"],
    "2687-GG-55" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-55_1_sequence.txt.gz"],
    "2687-GG-59" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-59_1_sequence.txt.gz"],
    "2687-GG-60" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-60_1_sequence.txt.gz"],
    "2687-GG-61" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-61_1_sequence.txt.gz"],
    "2687-GG-62" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-62_1_sequence.txt.gz"],
    "2687-GG-63" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-63_1_sequence.txt.gz"],
    "2687-GG-65" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-65_1_sequence.txt.gz"],
    "2687-GG-66" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-66_1_sequence.txt.gz"],
    "2687-GG-67" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-67_1_sequence.txt.gz"],
    "2687-GG-69" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-69_1_sequence.txt.gz"],
    "2687-GG-70" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-70_1_sequence.txt.gz"],
    "2687-GG-72" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-72_1_sequence.txt.gz"],
    "2687-GG-73" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-73_1_sequence.txt.gz"],
    "2687-GG-75" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-75_1_sequence.txt.gz"],
    "2687-GG-76" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-76_1_sequence.txt.gz"],
    "2687-GG-79" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-79_1_sequence.txt.gz"],
    "2687-GG-81" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-81_1_sequence.txt.gz"],
    "2687-GG-85" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-85_1_sequence.txt.gz"],
    }
};

performSmallRNAMouse($def_mouse);

1;

