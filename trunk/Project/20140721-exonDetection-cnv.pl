#!/usr/bin/perl
use strict;
use warnings;

use CQS::DNASeq;
use CQS::CNV;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir     = create_directory_or_die("/scratch/cqs/shengq1/exon_deletion");
my $bedfile        = "/scratch/cqs/lij17/cnv/SureSelect_XT_Human_All_Exon_V4_withoutchr_withoutY_lite.bed";
my $chromosome_dir = "/data/cqs/shengq1/reference/hg19_16569_MT_TCGA/chromosomes";
my $chromosome_len = "/data/cqs/shengq1/reference/hg19_16569_MT_TCGA/GRCh37-lite.len";
my $genome_fasta   = "/data/cqs/shengq1/reference/hg19_16569_MT_TCGA/GRCh37-lite.fa";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => "TCGA"
  },
  bamfiles => {
    "TCGA-A7-A0D9" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-A7-A0D9-01A-31W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B3" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0B3-01A-11W-A071-09_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0B8" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0B8-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BJ" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0BJ-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0BM" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0BM-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0C0" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0C0-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DK" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0DK-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0DP" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0DP-01A-21W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0E0" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0E0-01A-11W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
    "TCGA-BH-A0H7" => ["/gpfs21/scratch/cqs/shengq1/somaticmutation_comparison/bam/DNA_TP/TCGA-BH-A0H7-01A-13W-A071-09_HOLD_QC_PENDING_IlluminaGA-DNASeq_exome.bam"],
  },
  samtoolsindex => {
    class => "Samtools::Index",
    perform => 1,
    target_dir  => "${target_dir}/samtoolsindex",
    option      => "",
    source_ref  => "bamfiles",
    isbamsorted => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  cnvnator1000 => {
    target_dir     => "${target_dir}/cnvnator1000",
    option         => "",
    source_ref     => "bamfiles",
    binsize        => 1000,
    isbamsorted    => 1,
    chromosome_dir => "/scratch/cqs/shengq1/references/hg19chromosome",
    genome         => $genome_fasta,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  conifer => {
    class =>"CNV::Conifer",
    perform => 1,
    target_dir  => "${target_dir}/conifer",
    option      => "",
    source_ref  => "bamfiles",
    conifer     => "/home/shengq1/pylibs/bin/conifer.py",
    bedfile     => $bedfile,
    isbamsorted => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  cnmops => {
    class => "CNV::cnMops",
    perform => 1,
    target_dir  => "${target_dir}/cnmops",
    option      => "",
    source_ref  => "bamfiles",
    bedfile     => $bedfile,
    pairmode    => "paired",
    isbamsorted => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  freec => {
    class => "CNV::Freec",
    perform => 1,
    target_dir             => "${target_dir}/freec",
    option                 => "",
    source_ref             => "bamfiles",
    chrLenFile             => $chromosome_len,
    ploidy                 => 2,
    coefficientOfVariation => 0.01,
    chrFiles               => $chromosome_dir,
    inputFormat            => "BAM",
    mateOrientation        => "FR",

    #bedfile                => $bedfile, #provide bed file will make freec asking for the control samples
    pbs => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

samtools_index($config, "samtoolsindex");

cnvnator( $config, "cnvnator1000" );

conifer( $config, "conifer" );

cnmops( $config, "cnmops" );

freec( $config, "freec" );

1;
