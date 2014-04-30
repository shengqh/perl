#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20140415_bojana_FFPE_FF";

my $target_dir = "H:/shengquanhu/projects/Jennifer/20140415_bojana_MiSeq_FFPE_FF/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  groups  => {
    "FF"     => [ "IG-1",  "IG-7",  "IG-11", "IG-16", "IG-21", "IG-42", "IG-43", "IG-60", "IG-61" ],
    "FFPE"       => [ "IG-2",  "IG-8",  "IG-12", "IG-17", "IG-22", "IG-51", "IG-52", "IG-59", "IG-58" ],
    "FF_OLD" => [ "IG-1",  "IG-7",  "IG-11", "IG-16", "IG-21" ],
    "FFPE_OLD"   => [ "IG-2",  "IG-8",  "IG-12", "IG-17", "IG-22" ],
    "FF_NEW" => [ "IG-42", "IG-43", "IG-60", "IG-61" ],
    "FFPE_NEW"   => [ "IG-51", "IG-52", "IG-59", "IG-58" ],
    "FF_NO_IG1"     => [ "IG-7", "IG-11", "IG-16", "IG-21", "IG-42", "IG-43", "IG-60", "IG-61" ],
    "FFPE_NO_IG2"       => [ "IG-8", "IG-12", "IG-17", "IG-22", "IG-51", "IG-52", "IG-59", "IG-58" ],
    "FF_OLD_NO_IG1" => [ "IG-7", "IG-11", "IG-16", "IG-21" ],
    "FFPE_OLD_NO_IG2"   => [ "IG-8", "IG-12", "IG-17", "IG-22" ],
  },
  pairs => {
    "FFPE_VS_FF" => {
      groups => [ "FFPE", "FF" ],
      paired => 1
    },
    "FFPE_VS_FF_OLD" => {
      groups => [ "FFPE_OLD", "FF_OLD" ],
      paired => 1
    },
    "FFPE_VS_FF_NEW" => {
      groups => [ "FFPE_NEW", "FF_NEW" ],
      paired => 1
    },
    "FFPE_VS_FF_NO_IG2" => {
      groups => [ "FFPE_NO_IG2", "FF_NO_IG1" ],
      paired => 1
    },
    "FFPE_VS_FF_OLD_NO_IG2" => {
      groups => [ "FFPE_OLD_NO_IG2", "FF_OLD_NO_IG1" ],
      paired => 1
    },
  },
  deseq2 => {
    class      => "DESeq2",
    perform    => 1,
    target_dir => "${target_dir}/deseq2",
    option     => "",
    source_ref => "pairs",
    groups_ref => "groups",
    countfile  => "H:/shengquanhu/projects/Jennifer/20140415_bojana_MiSeq_FFPE_FF/IG-20140404-count_table-gene_symbol_FF_FFPE_matched.csv",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  }
};

performConfig($config);

1;
