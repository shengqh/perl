#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/papers/20130723_guoyan_QCpaper");

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "QC" },
  fastqfiles => {
    "1055QC0003-Sanger" => ["/scratch/cqs/shengq1/papers/20130723_guoyan_QCpaper/raw/1055QC0003-Sangerfq2.txt"],
    "2110-GOOD"         => ["/gpfs21/scratch/cqs/guoy1/2110/rawdata/2110-JP-1_2_sequence.txt"],
    "2110-BAD-1"        => ["/gpfs21/scratch/cqs/guoy1/2110/rawdata/2110-JP-43_2_sequence.txt"],
    "2110-BAD-2"        => ["/gpfs21/scratch/cqs/guoy1/2110/rawdata/2110-JP-15_1_sequence.txt"],
  },
  fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
