#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "motif";

my $target_dir = "H:/shengquanhu/projects/20160504_phospho";

my $config = {
  general => { task_name => $task },
  groups  => {
    "NORM" => [ "S17", "S23", "S29", "S35", "S40", "S45", "S49", "S52" ],
    "NAFL" => [ "S01", "S02", "S27", "S32", "S36", "S46", "S50", "S51" ],
    "NASH" => [ "S04", "S07", "S14", "S18", "S20", "S21", "S47", "S48" ]
  },
  pairs => {
    "NAFL_VS_NORM" => {
      groups => [ "NORM", "NAFL" ],
    },
    "NASH_VS_NORM" => {
      groups => [ "NORM", "NASH" ],
    },
    "NAFL_VS_NASH" => {
      groups => [ "NASH", "NAFL" ],
    },
  },
  site_deseq2 => {
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/site_deseq2",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile            => "H:/shengquanhu/projects/20160504_phospho/Motif_site_count.tsv",
    sh_direct            => 1,
    show_DE_gene_cluster => 1,
    pvalue               => 0.05,
    fold_change          => 1.5,
    min_median_read      => 2,
    pbs                  => {
      "email"    => "quanhu.sheng\@vanderbilt.edu",
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  spectrum_deseq2 => {
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/spectrum_deseq2",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile            => "H:/shengquanhu/projects/20160504_phospho/Motif_spectrum_count.tsv",
    sh_direct            => 1,
    show_DE_gene_cluster => 1,
    pvalue               => 0.05,
    fold_change          => 1.5,
    min_median_read      => 2,
    pbs                  => {
      "email"    => "quanhu.sheng\@vanderbilt.edu",
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
