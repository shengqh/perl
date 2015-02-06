#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20150130_bojana_HiSeq_FFPE_FF";

my $target_dir = "H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  groups  => {
    "FF" => [
      "2059-JP-10-1", "2059-JP-10-2", "2059-JP-16-1", "2059-JP-16-2", "2059-JP-17-1", "2059-JP-17-2", "2059-JP-10-3", "2059-JP-10-4", "2059-JP-10-5", "2059-JP-17-3",
      "2059-JP-17-4", "2059-JP-11-1", "2059-JP-11-2", "2059-JP-11-3", "2059-JP-16-3", "2059-JP-16-4", "2059-JP-16-5", "2059-JP-15-4", "2059-JP-15-5"
    ],
    "FFPE" => [
      "2059-JP-12-1", "2059-JP-12-2", "2059-JP-12-3", "2059-JP-12-4", "2059-JP-12-5", "2059-JP-13-1", "2059-JP-13-2", "2059-JP-13-3", "2059-JP-13-4", "2059-JP-13-5",
      "2059-JP-14-1", "2059-JP-14-3", "2059-JP-14-4", "2059-JP-14-5", "2059-JP-15-1", "2059-JP-15-2", "2059-JP-15-3", "2059-JP-11-4", "2059-JP-11-5"
    ],
    "FF2" => [
      "2059-JP-10-1", "2059-JP-10-2", "2059-JP-16-1", "2059-JP-16-2", "2059-JP-17-1", "2059-JP-17-2", "2059-JP-10-3", "2059-JP-10-4", "2059-JP-10-5", "2059-JP-17-3",
      "2059-JP-17-4", "2059-JP-11-1", "2059-JP-11-2", "2059-JP-11-3"
    ],
    "FFPE2" => [
      "2059-JP-12-1", "2059-JP-12-2", "2059-JP-12-3", "2059-JP-12-4", "2059-JP-12-5", "2059-JP-13-1", "2059-JP-13-2", "2059-JP-13-3", "2059-JP-13-4", "2059-JP-13-5",
      "2059-JP-14-1", "2059-JP-14-3", "2059-JP-14-4", "2059-JP-14-5"
    ],
  },
  pairs => {
    #"FFPE_VS_FF" => {
    #  groups => [ "FF",  "FFPE" ],
    #  paired => [ "P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20" ]
    #},
    "FFPE2_VS_FF2" => {
      groups => [ "FF2",  "FFPE2" ],
      paired => [ "P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P13", "P14", "P15" ]
    },
  },
  deseq2 => {
    class      => "DESeq2",
    perform    => 1,
    target_dir => "${target_dir}/deseq2",
    option     => "",
    source_ref => "pairs",
    groups_ref => "groups",
    countfile  => "H:/shengquanhu/projects/Jennifer/20150130_bojana_HiSeq_FFPE_FF/2059-JP-count_table-gene_symbol_FF_FFPE_matched.txt",
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
