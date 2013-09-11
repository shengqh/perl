#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;
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

###hg19.16571###
#my $fasta_file  = "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa";
#my $cosmic_file = "/data/cqs/shengq1/reference/cosmic/cosmic_v65_28052013.hg19.16571.vcf";
#my $snp_file    = "/data/cqs/shengq1/reference/snp137/hg19.16571/00-All.vcf";
#my $bowtie2_index = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";

my $config = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
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
  bam2fastq => {
    class      => "Bam2Fastq",
    perform    => 0,
    target_dir => "${target_dir}/bam2fastq",
    option     => "",
    source_ref => "rnafiles",
    cqstools   => $cqstools,
    samtools   => $samtools,
    ispaired   => 1,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  tophat2 => {
    class         => "Tophat2",
    perform       => 0,
    target_dir    => "${target_dir}/tophat2",
    option        => "--segment-length 25 -r 0 -p 6",
    source_ref    => "bam2fastq",
    bowtie2_index => $bowtie2_index,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  bwa => {
    class      => "BWA",
    perform    => 0,
    target_dir => "${target_dir}/bwa",
    option     => "-q 15 -t 8",
    fasta_file => $fasta_file,
    source_ref => "bam2fastq",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  refine => {
    class              => "GATKRefine",
    perform            => 0,
    target_dir         => "${target_dir}/refine",
    option             => "-Xmx40g",
    fasta_file         => $fasta_file,
    source_ref         => "bwa",
    thread_count       => 8,
    vcf_files          => [$snp_file],
    gatk_jar           => "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar",
    markDuplicates_jar => "/home/shengq1/local/bin/picard/MarkDuplicates.jar",
    sh_direct          => 1,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  rsmc_tophat2 => {
    class            => "RSMC",
    perform          => 1,
    target_dir       => "${target_dir}/rsmc_tophat2",
    option           => "-c 8",                                              #thread mode
    source_ref       => "tophat2",
    groups_ref       => "rnagroups",
    source_type      => "bam",                                               #source_type can be bam/mpileup
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
  rsmc_bwaRefine => {
    class            => "RSMC",
    perform          => 1,
    target_dir       => "${target_dir}/rsmc_bwaRefine",
    option           => "-c 8",                                              #thread mode
    source_ref       => "refine",
    groups_ref       => "rnagroups",
    source_type      => "bam",                                               #source_type can be bam/mpileup
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
  rsmc_nps => {
    class            => "RSMC",
    perform          => 1,
    target_dir       => "${target_dir}/rsmc_tophat2_not_filter_position_strand",
    option           => "-c 8 -n 0 -q 0 --max_normal_percentage 1.0 -g 0.0 -d 5 --not_filter_position --not_filter_strand ",          #thread mode
    source_ref       => "tophat2",
    groups_ref       => "rnagroups",
    source_type      => "bam",                                                     #source_type can be bam/mpileup
    fasta_file       => $fasta_file,
    annovar_buildver => "hg19",
    rnaediting_db    => "/data/cqs/shengq1/reference/rnaediting/hg19.txt",
    sh_direct        => 0,
    execute_file     => "/home/shengq1/rsmc/rsmc.exe",
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  muTect_tophat2 => {
    class        => "MuTect",
    perform      => 1,
    target_dir   => "${target_dir}/muTect_tophat2",
    option       => "",
    source_ref   => "tophat2",
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
  muTect_bwaRefine => {
    class        => "MuTect",
    perform      => 1,
    target_dir   => "${target_dir}/muTect_bwaRefine",
    option       => "",
    source_ref   => "refine",
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
  varscan2_tophat2 => {
    class           => "VarScan2",
    perform         => 1,
    target_dir      => "${target_dir}/varscan2_tophat2",
    option          => "",
    source_ref      => "tophat2",
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
  varscan2_bwaRefine => {
    class           => "VarScan2",
    perform         => 1,
    target_dir      => "${target_dir}/varscan2_bwaRefine",
    option          => "",
    source_ref      => "refine",
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
};

performConfig($config);

#rsmc( $config, "rsmc" );

1;
