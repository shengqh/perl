#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $target_dir = "/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53";

my $seed = "1414591741";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $impute2_option = "-seed ". $seed;
my $impute2_option_filter = $impute2_option . " -filt_rules_l 'eur.maf<0.01' 'afr.maf<0.01' ";

my $config = {
  general   => { task_name => "20140527_chenxi_pietenpol_p53" },
  ped_files => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.1.ped"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.3.ped"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.6.ped"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.7.ped"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.9.ped"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.10.ped"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.11.ped"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.12.ped"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.15.ped"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.16.ped"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.17.ped"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.18.ped"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.19.ped"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.20.ped"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.22.ped"]
  },
  gen_files => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.1.gen"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.3.gen"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.6.gen"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.7.gen"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.9.gen"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.10.gen"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.11.gen"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.12.gen"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.15.gen"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.16.gen"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.17.gen"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.18.gen"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.19.gen"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.20.gen"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.22.gen"]
  },
  haps_files => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.1.haps"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.3.haps"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.6.haps"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.7.haps"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.9.haps"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.10.haps"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.11.haps"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.12.haps"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.15.haps"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.16.haps"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.17.haps"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.18.haps"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.19.haps"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.20.haps"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/Pietenpol_p53.22.haps"]
  },
  genetic_map_files => {
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
  haplo_files => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr3.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr6.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr7.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr9.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr10.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr11.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr12.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr15.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr16.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr17.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr18.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr19.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr20.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/ref_panel/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes"],
  },
  shapeit => {
    class                => "Imputation::Shapeit",
    perform              => 0,
    path_file            => "/home/shengq1/local/bin/path_glibc2.14.txt",
    target_dir           => "${target_dir}/shapeit_gen",
    option               => "--aligned -T 8 --seed " . $seed,
    source_ref           => "gen_files",
    genetic_map_file_ref => "genetic_map_files",
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  shapeit_impute2 => {
    class                 => "Imputation::Impute2",
    perform               => 0,
    target_dir            => "${target_dir}/shapeit_impute2",
    option                => $impute2_option_filter . " -use_prephased_g",
    max_chromosome_length => "250000000",
    interval              => "5000000",
    source_ref            => [ "shapeit", "haps\$" ],
    genetic_map_file_ref  => "genetic_map_files",
    haplo_file_ref        => "haplo_files",
    isPhased              => 1,
    sh_direct             => 0,
    pbs                   => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  impute2_direct => {
    class                 => "Imputation::Impute2",
    perform               => 1,
    target_dir            => "${target_dir}/impute2_direct",
    option                => $impute2_option,
    max_chromosome_length => "250000000",
    interval              => "5000000",
    source_ref            => "gen_files",
    genetic_map_file_ref  => "genetic_map_files",
    haplo_file_ref        => "haplo_files",
    isPhased              => 0,
    sh_direct             => 0,
    pbs                   => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  impute2_direct_filter => {
    class                 => "Imputation::Impute2",
    perform               => 1,
    target_dir            => "${target_dir}/impute2_direct_filter",
    option                => $impute2_option_filter,
    max_chromosome_length => "250000000",
    interval              => "5000000",
    source_ref            => "gen_files",
    genetic_map_file_ref  => "genetic_map_files",
    haplo_file_ref        => "haplo_files",
    isPhased              => 0,
    sh_direct             => 0,
    pbs                   => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  impute2_direct_missing_filter => {
    class                 => "Imputation::Impute2",
    perform               => 1,
    target_dir            => "${target_dir}/impute2_direct_filter_missing",
    option                => $impute2_option_filter,
    max_chromosome_length => "250000000",
    interval              => "5000000",
    source_ref            => "haps_files",
    genetic_map_file_ref  => "genetic_map_files",
    haplo_file_ref        => "haplo_files",
    isPhased              => 1,
    sh_direct             => 0,
    pbs                   => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  }
};

#performConfig($config);
performTask( $config, "shapeit" );


1;

