#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

#my $target_dir      = "E:/temp";
my $target_dir      = "/scratch/cqs/shengq1/proteomics/20151103_shilin_labelfree_elution_order_algorithm";
my $msgf_jar        = "/scratch/cqs/shengq1/local/bin/MSGFPlus/MSGFPlus.jar";
my $proteomicstools = "/home/shengq1/proteomicstools/ProteomicsTools.exe";
my $mod_file        = "${target_dir}/config/msgf_mods.txt";
my $database_yeast  = "${target_dir}/database/yeast_20150113_orf_trans.REVERSED.fasta";
my $email           = "quanhu.sheng\@vanderbilt.edu";

my $buildsummary_msgf_yeast_file = "${target_dir}/config/buildsummary_msgf_yeast.param";

my $datasets = {
  Elite_HCDOT_Yeast => {
    source => {
      "21July2013_TWR_Yeast" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20151103_shilin_labelfree_elution_order_algorithm/Elite_HCDOT_Yeast/jcoon_figure1_mgf/21July2013_TWR_Yeast.mgf"],
      "22July2013_TWR_Yeast" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20151103_shilin_labelfree_elution_order_algorithm/Elite_HCDOT_Yeast/jcoon_figure1_mgf/22July2013_TWR_Yeast.mgf"],
      "24July2013_TWR_Yeast" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20151103_shilin_labelfree_elution_order_algorithm/Elite_HCDOT_Yeast/jcoon_figure1_mgf/24July2013_TWR_Yeast_Night.mgf"],
      "25July2013_TWR_Yeast" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20151103_shilin_labelfree_elution_order_algorithm/Elite_HCDOT_Yeast/jcoon_figure1_mgf/25July2013_TWR_Yeast.mgf"],
      "27July2013_TWR_Yeast" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20151103_shilin_labelfree_elution_order_algorithm/Elite_HCDOT_Yeast/jcoon_figure1_mgf/27July2013_TWR_Yeast.mgf"],
      "31July2013_TWR_yeast" => ["/gpfs21/scratch/cqs/shengq1/proteomics/20151103_shilin_labelfree_elution_order_algorithm/Elite_HCDOT_Yeast/jcoon_figure1_mgf/31July2013_TWR_yeast.mgf"],
    },
    groups => { "TWR" => [ "21July2013_TWR_Yeast", "22July2013_TWR_Yeast", "24July2013_TWR_Yeast", "25July2013_TWR_Yeast", "27July2013_TWR_Yeast", "31July2013_TWR_yeast" ] },
    msgf_option                   => "-t 10ppm -ti \"0,1\" -tda 0 -m 3 -inst 1 -e 1 -protocol 5 -ntt 2 -n 2 -addFeatures 1",
    buildsummary_msgf_config_file => $buildsummary_msgf_yeast_file,
    database                      => $database_yeast
  },
};

my $configall = { general => { task_name => "yeast" }, };

for my $dataset ( sort keys %{$datasets} ) {
  my $def    = $datasets->{$dataset};
  my $curdir = create_directory_or_die("${target_dir}/$dataset/");

  my $config = {
    "${dataset}_msgf" => {
      task_name  => "${dataset}_msgf",
      class      => "Proteomics::Engine::MSGFPlus",
      perform    => 1,
      target_dir => "${curdir}/msgf",
      option     => $def->{msgf_option},
      source     => $def->{source},
      msgf_jar   => $msgf_jar,
      mod_file   => $mod_file,
      database   => $def->{database},
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=8",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    },
    "${dataset}_msgf_buildsummary" => {
      task_name       => "${dataset}_msgf",
      class           => "Proteomics::Summary::BuildSummary",
      perform         => 1,
      target_dir      => "${curdir}/buildsummary",
      option          => "",
      source_ref      => "${dataset}_msgf",
      groups          => $def->{groups},
      parameter_file  => $def->{buildsummary_msgf_config_file},
      proteomicstools => $proteomicstools,
      sh_direct       => 1,
      pbs             => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    }
  };

  $configall = merge( $configall, $config );
}

#print Dumper($configall);
performConfig($configall);

1;
