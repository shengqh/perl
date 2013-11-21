#!/usr/bin/perl
use strict;
use warnings;

my  @rna_groups = ("TCGA-A7-A0D9-RNA-TP-NT",
    "TCGA-BH-A0B3-RNA-TP-NT", 
    "TCGA-BH-A0B8-RNA-TP-NT",
    "TCGA-BH-A0BJ-RNA-TP-NT",
    "TCGA-BH-A0BM-RNA-TP-NT",
    "TCGA-BH-A0C0-RNA-TP-NT",
    "TCGA-BH-A0DK-RNA-TP-NT",
    "TCGA-BH-A0DP-RNA-TP-NT",
    "TCGA-BH-A0E0-RNA-TP-NT",
    "TCGA-BH-A0H7-RNA-TP-NT");
    
my $dir = "/scratch/cqs/shengq1/somaticmutation_comparison/16569_rsmc/result";
  
for my $key (@rna_groups){
  print($key . "\n");
  system("mono-sgen /home/shengq1/rsmc/rsmc.exe filter -i ${dir}/${key}/${key}.bases -o ${dir}/${key}/${key}.tsv 1> ${dir}/log/${key}.log");
}  
  
1;
