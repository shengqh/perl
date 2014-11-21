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

my $impute2_option        = "-seed " . $seed;
my $impute2_option_filter = $impute2_option . " -filt_rules_l 'eur.maf<0.01' 'afr.maf<0.01' ";

my $gens = {
  "660w" => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.01.gen"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.03.gen"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.06.gen"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.07.gen"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.09.gen"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.10.gen"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.11.gen"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.12.gen"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.15.gen"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.16.gen"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.17.gen"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.18.gen"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.19.gen"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.20.gen"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.660W.22.gen"],
  },
  "1M" => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.01.gen"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.03.gen"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.06.gen"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.07.gen"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.09.gen"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.10.gen"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.11.gen"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.12.gen"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.15.gen"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.16.gen"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.17.gen"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.18.gen"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.19.gen"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.20.gen"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.1M.22.gen"],
  },
  "OMNI" => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.01.gen"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.03.gen"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.06.gen"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.07.gen"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.09.gen"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.10.gen"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.11.gen"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.12.gen"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.15.gen"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.16.gen"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.17.gen"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.18.gen"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.19.gen"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.20.gen"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.OMNI.22.gen"],
  },
  "Omni5" => {
    "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.01.gen"],
    "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.03.gen"],
    "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.06.gen"],
    "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.07.gen"],
    "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.09.gen"],
    "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.10.gen"],
    "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.11.gen"],
    "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.12.gen"],
    "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.15.gen"],
    "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.16.gen"],
    "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.17.gen"],
    "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.18.gen"],
    "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.19.gen"],
    "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.20.gen"],
    "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/preimputation/byplatforms/Pietenpol_p53.Omni5.22.gen"],
  },
};

for my $platform ( sort keys %{$gens} ) {
  my $config = {
    general    => { task_name => "20140527_chenxi_pietenpol_p53" },
    gen_files  => $gens->{$platform},
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

    #for i in $(ls -d */); do
    #        for j in $(ls ${i}*.haps); do
    #                grep "ref:" ${j} >> ${i}/${i%%/}.ref.haps
    #        done
    #done
    #cqstools file_def -i . -r -f .ref.haps$ -n \(.+\).ref
    haps_ref_files => {
      "Pietenpol_p53.01" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.01/Pietenpol_p53.01.ref.haps"],
      "Pietenpol_p53.03" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.03/Pietenpol_p53.03.ref.haps"],
      "Pietenpol_p53.06" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.06/Pietenpol_p53.06.ref.haps"],
      "Pietenpol_p53.07" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.07/Pietenpol_p53.07.ref.haps"],
      "Pietenpol_p53.09" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.09/Pietenpol_p53.09.ref.haps"],
      "Pietenpol_p53.10" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.10/Pietenpol_p53.10.ref.haps"],
      "Pietenpol_p53.11" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.11/Pietenpol_p53.11.ref.haps"],
      "Pietenpol_p53.12" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.12/Pietenpol_p53.12.ref.haps"],
      "Pietenpol_p53.15" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.15/Pietenpol_p53.15.ref.haps"],
      "Pietenpol_p53.16" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.16/Pietenpol_p53.16.ref.haps"],
      "Pietenpol_p53.17" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.17/Pietenpol_p53.17.ref.haps"],
      "Pietenpol_p53.18" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.18/Pietenpol_p53.18.ref.haps"],
      "Pietenpol_p53.19" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.19/Pietenpol_p53.19.ref.haps"],
      "Pietenpol_p53.20" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.20/Pietenpol_p53.20.ref.haps"],
      "Pietenpol_p53.22" => ["/gpfs21/scratch/cqs/shengq1/chenxi/20140527_chenxi_pietenpol_p53/shapeit_gen/result/Pietenpol_p53.22/Pietenpol_p53.22.ref.haps"],
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
      perform              => 1,
      path_file            => "/home/shengq1/local/bin/path_glibc2.14.txt",
      target_dir           => "${target_dir}/shapeit_gen_" . $platform,
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
      target_dir            => "${target_dir}/shapeit_gen_impute2_" . $platform,
      option                => $impute2_option_filter . " -use_prephased_g",
      max_chromosome_length => "250000000",
      interval              => "5000000",
      source_ref            => "haps_ref_files",
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
      perform               => 0,
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
      perform               => 0,
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
      perform               => 0,
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

  performConfig($config);

  #performTask( $config, "shapeit_impute2" );
}

1;

