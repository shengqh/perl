#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = "/scratch/cqs/shengq1/dnaseq/20160715_exomeseq_vcf_pipeline";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools   = "/home/shengq1/cqstools/CQS.Tools.exe";
my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";

my $bwa_fasta = "/scratch/cqs/shengq1/references/hg19_16569_MT/bwa_index_0.7.12/hg19_16569_MT.fa";
my $dbsnp     = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $hapmap    = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni      = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000     = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills     = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $covered_bed = "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/0699701_Covered.bed";

my $config = {
  general => { task_name => "exomeseq" },

  files => {
    "SA001" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_1/1_CATCAAGT_L007_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_1/1_CATCAAGT_L007_R2_001.fastq.gz"
    ],
    "SA003BP" => [
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_3BP/3BP_TCTTCACA_L008_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/raw/Sample_3BP/3BP_TCTTCACA_L008_R2_001.fastq.gz"
    ],
  },
  bwa => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "",
    bwa_index  => $bwa_fasta,
    source_ref => "files",
    picard_jar => $picard_jar,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine => {
    class            => "GATK::Refine",
    perform          => 1,
    target_dir       => "${target_dir}/bwa_refine",
    option           => "-Xmx40g",
    gatk_option      => "--fix_misencoded_quality_scores",
    fasta_file       => $bwa_fasta,
    source_ref       => "bwa",
    indel_vcf_files  => [ $mills, $g1000 ],
    vcf_files        => [$dbsnp],
    bed_file         => $covered_bed,
    interval_padding => 100,
    gatk_jar         => $gatk_jar,
    picard_jar       => $picard_jar,
    sh_direct        => 0,
    sorted           => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf => {
    class            => "GATK::HaplotypeCallerGVCF",
    perform          => 1,
    target_dir       => "${target_dir}/bwa_refine_hc_gvcf",
    option           => "",
    source_ref       => "bwa_refine",
    java_option      => "",
    fasta_file       => $bwa_fasta,
    dbsnp_vcf        => $dbsnp,
    bed_file         => $covered_bed,
    interval_padding => 100,
    gatk_jar         => $gatk_jar,
    extension        => ".gvcf",
    sh_direct        => 0,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_vqsr",
    option      => "",
    source_ref  => "bwa_refine_hc_gvcf",
    java_option => "",
    fasta_file  => $bwa_fasta,
    dbsnp_vcf   => $dbsnp,
    hapmap_vcf  => $hapmap,
    omni_vcf    => $omni,
    g1000_vcf   => $g1000,
    mills_vcf   => $mills,
    bed_file    => $covered_bed,
    gatk_jar    => $gatk_jar,
    cqstools    => $cqstools,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "bwa", "bwa_refine", "bwa_refine_hc_gvcf" ],
      step2 => ["bwa_refine_hc_gvcf_vqsr"],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;

