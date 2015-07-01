#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20150630_bojana_tnbc";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => { "samples" => ["/scratch/cqs/shengq1/rnaseq/20150630_bojana_tnbc/raw/result/fastq.list"], },
  wget    => {
    class      => "Data::Wget",
    perform    => 1,
    target_dir => "${target_dir}/raw",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;

