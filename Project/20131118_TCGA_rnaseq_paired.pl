#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20131118_rnaseq_paired");

my $transcript_gtf       = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.gtf";
my $hg19_map = "/scratch/cqs/shengq1/references/hg19/Homo_sapiens.GRCh37.73.map";
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "TCGA";

my $config = {
  general => { task_name => $task },
  files => {
    "TCGA-A7-A0D9-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-A7-A0D9-RNA-NT/TCGA-A7-A0D9-RNA-NT.bam"],
    "TCGA-A7-A0D9-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-A7-A0D9-RNA-TP/TCGA-A7-A0D9-RNA-TP.bam"],
    "TCGA-BH-A0B3-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0B3-RNA-NT/TCGA-BH-A0B3-RNA-NT.bam"],
    "TCGA-BH-A0B3-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0B3-RNA-TP/TCGA-BH-A0B3-RNA-TP.bam"],
    "TCGA-BH-A0B8-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0B8-RNA-NT/TCGA-BH-A0B8-RNA-NT.bam"],
    "TCGA-BH-A0B8-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0B8-RNA-TP/TCGA-BH-A0B8-RNA-TP.bam"],
    "TCGA-BH-A0BJ-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0BJ-RNA-NT/TCGA-BH-A0BJ-RNA-NT.bam"],
    "TCGA-BH-A0BJ-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0BJ-RNA-TP/TCGA-BH-A0BJ-RNA-TP.bam"],
    "TCGA-BH-A0BM-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0BM-RNA-NT/TCGA-BH-A0BM-RNA-NT.bam"],
    "TCGA-BH-A0BM-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0BM-RNA-TP/TCGA-BH-A0BM-RNA-TP.bam"],
    "TCGA-BH-A0C0-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0C0-RNA-NT/TCGA-BH-A0C0-RNA-NT.bam"],
    "TCGA-BH-A0C0-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0C0-RNA-TP/TCGA-BH-A0C0-RNA-TP.bam"],
    "TCGA-BH-A0DK-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0DK-RNA-NT/TCGA-BH-A0DK-RNA-NT.bam"],
    "TCGA-BH-A0DK-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0DK-RNA-TP/TCGA-BH-A0DK-RNA-TP.bam"],
    "TCGA-BH-A0DP-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0DP-RNA-NT/TCGA-BH-A0DP-RNA-NT.bam"],
    "TCGA-BH-A0DP-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0DP-RNA-TP/TCGA-BH-A0DP-RNA-TP.bam"],
    "TCGA-BH-A0E0-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0E0-RNA-NT/TCGA-BH-A0E0-RNA-NT.bam"],
    "TCGA-BH-A0E0-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0E0-RNA-TP/TCGA-BH-A0E0-RNA-TP.bam"],
    "TCGA-BH-A0H7-RNA-NT" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0H7-RNA-NT/TCGA-BH-A0H7-RNA-NT.bam"],
    "TCGA-BH-A0H7-RNA-TP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/rna_tophat2/result/TCGA-BH-A0H7-RNA-TP/TCGA-BH-A0H7-RNA-TP.bam"],
  },
  sortbam => {
    class         => "Sortbam",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "files",
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseqcount => {
    class      => "HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/htseqcount",
    option     => "",
    source_ref => "sortbam",
    gff_file   => $transcript_gtf,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  genetable => {
    class         => "CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/genetable",
    option        => "-p ENS --noheader -o ${task}_gene.count",
    source_ref    => "htseqcount",
    name_map_file => $hg19_map,
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
