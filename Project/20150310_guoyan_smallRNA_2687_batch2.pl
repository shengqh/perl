#!/usr/bin/perl
use strict;
use warnings;

use CQS::PerformSmallRNA;

my $def_human = {

  #General options
  task_name       => "2687_b2",
  email           => "quanhu.sheng\@vanderbilt.edu",
  target_dir      => "/scratch/cqs/shengq1/smallRNA/20150310_guoyan_2687_human_mouse_batch2/human",
  max_thread      => 8,
  min_read_length => 16,
  cluster         => "slurm",

  #Data
  files => {
  "2687-GG-095" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-95_1_sequence.txt.gz"],
  "2687-GG-096" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-96_1_sequence.txt.gz"],
  "2687-GG-100" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-100_1_sequence.txt.gz"],
  "2687-GG-101" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-101_1_sequence.txt.gz"],
  "2687-GG-102" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-102_1_sequence.txt.gz"],
  "2687-GG-103" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-103_1_sequence.txt.gz"],
  "2687-GG-104" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-104_1_sequence.txt.gz"],
  "2687-GG-105" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-105_1_sequence.txt.gz"],
  "2687-GG-106" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-106_1_sequence.txt.gz"],
  "2687-GG-107" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-107_1_sequence.txt.gz"],
  "2687-GG-113" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-113_1_sequence.txt.gz"],
  "2687-GG-114" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-114_1_sequence.txt.gz"],
  "2687-GG-115" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-115_1_sequence.txt.gz"],
  "2687-GG-116" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-116_1_sequence.txt.gz"],
  "2687-GG-117" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-117_1_sequence.txt.gz"],
  "2687-GG-118" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-118_1_sequence.txt.gz"],
  "2687-GG-119" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-119_1_sequence.txt.gz"],
  "2687-GG-120" => ["/gpfs21/scratch/cqs/guoy1/2687/2687-GG-120_1_sequence.txt.gz"],
  }
};

performSmallRNA_hg19($def_human);

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

#performSmallRNA_mm10($def_mouse);


1;

