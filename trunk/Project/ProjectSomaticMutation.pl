#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/somaticmutation_comparison");

my $email    = "quanhu.sheng\@vanderbilt.edu";
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";
my $samtools = "/home/shengq1/local/bin/samtools/samtools";

##hg19.16569###
my $fasta_file    = "/data/cqs/shengq1/reference/hg19.16569/bowtie2_index_2.1.0/hg19_rCRS.fa";
my $cosmic_file   = "/data/cqs/shengq1/reference/cosmic/cosmic_v66_20130725.hg19.16569.vcf";
my $snp_file      = "/data/cqs/shengq1/reference/snp137/hg19.16569/dbsnp_137.b37.vcf";
my $bowtie2_index = "/data/cqs/shengq1/reference/hg19.16569/bowtie2_index_2.1.0/hg19_rCRS";
my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $config = {
  general => {
    task_name => "somaticmutation"
  },
  rnafiles => {
    "TCGA-A7-A0D9-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-A7-A0D9-RNA_NT_sorted.bam"],
    "TCGA-A7-A0D9-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-A7-A0D9-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B3-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B3-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B3-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B3-RNA_TP_sorted.bam"],
    "TCGA-BH-A0B8-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B8-RNA_NT_sorted.bam"],
    "TCGA-BH-A0B8-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0B8-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BJ-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BJ-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BJ-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BJ-RNA_TP_sorted.bam"],
    "TCGA-BH-A0BM-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BM-RNA_NT_sorted.bam"],
    "TCGA-BH-A0BM-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0BM-RNA_TP_sorted.bam"],
    "TCGA-BH-A0C0-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0C0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0C0-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0C0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DK-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DK-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DK-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DK-RNA_TP_sorted.bam"],
    "TCGA-BH-A0DP-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DP-RNA_NT_sorted.bam"],
    "TCGA-BH-A0DP-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0DP-RNA_TP_sorted.bam"],
    "TCGA-BH-A0E0-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0E0-RNA_NT_sorted.bam"],
    "TCGA-BH-A0E0-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0E0-RNA_TP_sorted.bam"],
    "TCGA-BH-A0H7-NT" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0H7-RNA_NT_sorted.bam"],
    "TCGA-BH-A0H7-TP" => ["/scratch/cqs/shengq1/somaticmutation_comparison/bam/TCGA-BH-A0H7-RNA_TP_sorted.bam"],
  },
  rnagroups => {
    "TCGA-A7-A0D9" => [ "TCGA-A7-A0D9-NT", "TCGA-A7-A0D9-TP" ],
    "TCGA-BH-A0B3" => [ "TCGA-BH-A0B3-NT", "TCGA-BH-A0B3-TP" ],
    "TCGA-BH-A0B8" => [ "TCGA-BH-A0B8-NT", "TCGA-BH-A0B8-TP" ],
    "TCGA-BH-A0BJ" => [ "TCGA-BH-A0BJ-NT", "TCGA-BH-A0BJ-TP" ],
    "TCGA-BH-A0BM" => [ "TCGA-BH-A0BM-NT", "TCGA-BH-A0BM-TP" ],
    "TCGA-BH-A0C0" => [ "TCGA-BH-A0C0-NT", "TCGA-BH-A0C0-TP" ],
    "TCGA-BH-A0DK" => [ "TCGA-BH-A0DK-NT", "TCGA-BH-A0DK-TP" ],
    "TCGA-BH-A0DP" => [ "TCGA-BH-A0DP-NT", "TCGA-BH-A0DP-TP" ],
    "TCGA-BH-A0E0" => [ "TCGA-BH-A0E0-NT", "TCGA-BH-A0E0-TP" ],
    "TCGA-BH-A0H7" => [ "TCGA-BH-A0H7-NT", "TCGA-BH-A0H7-TP" ],
  },
  rsmc_nofilter => {
    class            => "RSMC",
    perform          => 0,
    target_dir       => "${target_dir}/rsmc_nofilter",
    option           => "-c 8 -n 20 -q 0 --max_normal_percentage 1.0 -g 0.0 -d 5 --not_filter_position --not_filter_strand ",    #thread mode
    source_ref       => "rnafiles",
    groups_ref       => "rnagroups",
    source_type      => "bam",                                                                                                   #source_type can be bam/mpileup
    fasta_file       => $fasta_file,
    annovar_buildver => "hg19",
    rnaediting_db    => "/data/cqs/shengq1/reference/rnaediting/hg19.txt",
    sh_direct        => 1,
    execute_file     => "/home/shengq1/rsmc/rsmc.exe",
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  muTect => {
    class        => "MuTect",
    perform      => 0,
    target_dir   => "${target_dir}/muTect",
    option       => "",
    source_ref   => "rnafiles",
    groups_ref   => "rnagroups",
    java_option  => "-Xmx40g",
    fasta_file   => $fasta_file,
    cosmic_file  => $cosmic_file,
    dbsnp_file   => $snp_file,
    bychromosome => 0,
    sh_direct    => 0,
    muTect_jar   => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  annovar_mutect => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/muTect",
    option     => $annovar_param,
    source_ref => [ "muTect", ".pass.vcf\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    cqstools   => $cqstools,
    affy_file  => "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  varscan2 => {
    class           => "VarScan2",
    perform         => 0,
    target_dir      => "${target_dir}/varscan2",
    option          => "",
    source_ref      => "rnafiles",
    groups_ref      => "rnagroups",
    fasta_file      => $fasta_file,
    min_coverage    => 10,
    somatic_p_value => 0.01,
    sh_direct       => 0,
    VarScan2_jar    => "/home/shengq1/local/bin/VarScan.v2.3.5.jar",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  annovar_varscan2 => {
    class      => "Annovar",
    perform    => 1,
    target_dir => "${target_dir}/varscan2",
    option     => $annovar_param,
    source_ref => [ "varscan2", ".Somatic.hc\$" ],
    annovar_db => $annovar_db,
    buildver   => "hg19",
    cqstools   => $cqstools,
    affy_file  => "/data/cqs/shengq1/reference/affy/HG-U133_Plus_2.na33.annot.csv",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
