#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/home/cylee/tiger/wgs");
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools   = "/home/cylee/tools/cqstools/cqstools.exe";
my $gatk_jar   = "/home/cylee/gatk/GenomeAnalysisTK.jar";
my $picard_jar = "/home/cylee/tools/picard.jar";

my $bwa_fasta = "/home/cylee/tiger/bundle/b37/bwa_index_0.7.12/human_g1k_v37.fasta";
my $dbsnp     = "/home/cylee/tiger/bundle/dbsnp/human_GRCh37_v147_16569_MT.vcf";
my $hapmap    = "/home/cylee/tiger/bundle/b37/hapmap_3.3.b37.vcf";
my $omni      = "/home/cylee/tiger/bundle/b37/1000G_omni2.5.b37.vcf";
my $g1000     = "/home/cylee/tiger/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills     = "/home/cylee/tiger/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $config = {
  general => { task_name => "wgs" },

  files => {
    "Control" => ["/home/cylee/tiger/wgs/data/SRR3713934.sra"],
    "Sample"  => ["/home/cylee/tiger/wgs/data/SRR3713939.sra"],
  },
  sra2fastq => {
    class      => "SRA::FastqDump",
    perform    => 1,
    ispaired   => 0,
    target_dir => "${target_dir}/sra2fastq",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  bwa => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "",
    bwa_index  => $bwa_fasta,
    source_ref => "sra2fastq",
    picard_jar => $picard_jar,
    sh_direct  => 1,
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
    indel_vcf_files  => [$mills],
    known_vcf_files  => [ $dbsnp, $g1000 ],
    interval_padding => 100,
    gatk_jar         => $gatk_jar,
    picard_jar       => $picard_jar,
    sh_direct        => 1,
    sorted           => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf => {
    class            => "GATK::HaplotypeCaller",
    perform          => 1,
    target_dir       => "${target_dir}/bwa_refine_hc_gvcf",
    option           => "",
    source_ref       => "bwa_refine",
    java_option      => "",
    fasta_file       => $bwa_fasta,
    dbsnp_vcf        => $dbsnp,
    interval_padding => 100,
    gatk_jar         => $gatk_jar,
    vqsr_mode        => 1,
    extension        => ".gvcf",
    sh_direct        => 1,
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
    gatk_jar    => $gatk_jar,
    vqsr_mode   => 1,
    cqstools    => $cqstools,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_hardFilter => {
    class       => "GATK::VariantFilter",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf_hardFilter",
    option      => "",
    source_ref  => "bwa_refine_hc_gvcf",
    java_option => "",
    fasta_file  => $bwa_fasta,
    dbsnp_vcf   => $dbsnp,
    is_rna      => 0,
    gatk_jar    => $gatk_jar,
    cqstools    => $cqstools,
    vqsr_mode   => 0,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=12",
      "walltime" => "4",
      "mem"      => "70gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "sra2fastq",               "bwa", "bwa_refine", "bwa_refine_hc_gvcf" ],
      step2 => [ "bwa_refine_hc_gvcf_vqsr", "bwa_refine_hc_gvcf_hardFilter" ],
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


