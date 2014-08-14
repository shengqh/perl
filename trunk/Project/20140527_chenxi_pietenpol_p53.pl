#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $target_dir = "/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general   => { task_name => "20140527_chenxi_pietenpol_p53" },
  ped_files => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.1.ped"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.3.ped"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.6.ped"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.7.ped"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.9.ped"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.10.ped"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.11.ped"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.12.ped"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.15.ped"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.16.ped"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.17.ped"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.18.ped"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.19.ped"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.20.ped"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/data_by_chrom/Pietenpol_p53.22.ped"]
  },
  map_files => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr1_combined_b37.txt"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr3_combined_b37.txt"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr6_combined_b37.txt"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr7_combined_b37.txt"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr9_combined_b37.txt"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr10_combined_b37.txt"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr11_combined_b37.txt"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr12_combined_b37.txt"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr15_combined_b37.txt"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr16_combined_b37.txt"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr17_combined_b37.txt"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr18_combined_b37.txt"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr19_combined_b37.txt"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr20_combined_b37.txt"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/genetic_map_chr22_combined_b37.txt"],
  },
  shapeit => {
    class        => "Imputation::Shapeit",
    perform      => 1,
    path_file => "/home/shengq1/local/bin/path_glibc2.14.txt",
    target_dir   => "${target_dir}/shapeit",
    option       => "",
    source_ref   => "ped_files",
    map_file_ref => "map_files",
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  }
};

performConfig($config);

1;