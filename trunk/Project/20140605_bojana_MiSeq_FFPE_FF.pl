#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20140605_bojana_FFPE_FF";

my $target_dir =
  "H:/shengquanhu/projects/Jennifer/20140605_bojana_MiSeq_FFPE_FF/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  groups  => {
    "FF" => [
      "IG-40", "IG-41", "IG-42", "IG-43", "IG-44", "IG-45", "IG-46", "IG-47", "IG-61", "IG-60"
    ],
    "FFPE" => [
      "IG-49", "IG-50", "IG-51", "IG-52", "IG-53", "IG-54", "IG-55", "IG-56", "IG-58", "IG-59"
    ],
    "FF_NO_P11" => [
      "IG-40", "IG-41", "IG-42", "IG-43", "IG-44", "IG-45", "IG-46", "IG-61", "IG-60"
    ],
    "FFPE_NO_P11" => [
      "IG-49", "IG-50", "IG-51", "IG-52", "IG-53", "IG-54", "IG-55", "IG-58", "IG-59"
    ],
    "FF_NEW_SAMPLE" => [
      "IG-40", "IG-41", "IG-44", "IG-45", "IG-46"
    ],
    "FFPE_NEW_SAMPLE" => [
      "IG-49", "IG-50", "IG-53", "IG-54", "IG-55"
    ],
  },
  pairs => {
#    "FFPE_VS_FF" => {
#      groups => [ "FFPE", "FF" ],
#      paired => [ "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P13", "P14" ]
#    },
    "FFPE_VS_FF_NO_P11" => {
      groups => [ "FFPE_NO_P11", "FF_NO_P11" ],
      paired => [ "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P13", "P14" ]
    },
    "FFPE_VS_FF_NEWSAMPLE" => {
      groups => [ "FFPE_NEW_SAMPLE", "FF_NEW_SAMPLE" ],
      paired => [ "P04", "P05", "P08", "P09", "P10" ]
    },
  },
  deseq2 => {
    class      => "DESeq2",
    perform    => 1,
    target_dir => "${target_dir}/deseq2",
    option     => "",
    source_ref => "pairs",
    groups_ref => "groups",
    countfile  => "H:/shengquanhu/projects/Jennifer/20140605_bojana_MiSeq_FFPE_FF/20140605_bojana_MiSeq_FFPE_FF.csv",
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  }
};

performConfig($config);

1;
