#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "phos";

my $target_dir = "E:/shengquanhu/projects/20150807_VU_Phospho_Labelfree/deseq2";

#my $target_dir = "e:/temp";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  groups  => {
    "AM" => [ "AM291", "AM297", "AM298", "AM301", "AM303", "AM309", ],
    "UM" => [ "UM291", "UM297", "UM298", "UM301", "UM303", "UM309" ],
  },
  pairs => {
    "AM_VS_UM" => {
      groups => [ "UM",  "AM" ],
      paired => [ "P291", "P297", "P298", "P301", "P303", "P309" ]
    },
  },
  deseq2_rawcount => {
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/deseq2_rawcount",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile            => "E:/shengquanhu/projects/20150807_VU_Phospho_Labelfree/deseq2/raw_count.csv",
    sh_direct            => 1,
    show_DE_gene_cluster => 0,
    pvalue               => 0.05,
    fold_change          => 2.0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  deseq2_adjustedcount => {
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/deseq2_adjustedcount",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile            => "E:/shengquanhu/projects/20150807_VU_Phospho_Labelfree/deseq2/adjusted_count.csv",
    sh_direct            => 1,
    show_DE_gene_cluster => 0,
    pvalue               => 0.05,
    fold_change          => 2.0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

#performTask( $config, "star_deseq2" );
#performTask( $config, "star_deseq2_strict_criteria" );

#performTask( $config, "tophat2_genetable" );

1;

