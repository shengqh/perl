#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $target_dir = "/scratch/cqs/shengq1/dnaseq/20160829_liuqi_gene_panel";
my $cqstools   = "/home/shengq1/cqstools/cqstools.exe";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $bwa_fasta      = "/scratch/cqs/shengq1/references/gatk/b37/bwa_index_0.7.12/human_g1k_v37.fasta";
my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.map";

my $dbsnp  = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $hapmap = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni   = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000  = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills  = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";

my $cosmic    = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_MT.vcf";

my $annovar_param = "-protocol refGene,snp138,cosmic70 -operation g,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";
my $gatk_jar      = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar    = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";

my $covered_bed   = "/gpfs21/scratch/cqs/shengq1/dnaseq/20150406_bojana_dnaseq_selectedgenes/0699701_Covered.bed";

my $qc3_perl = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $cluster = "slurm";

my $config = {
  general => { task_name => "liuqi_gene" },

  files => {
  },
  groups => {
  },

  fastqc => {
    class      => "QC::FastQC",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    cluster    => $cluster,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    cluster    => $cluster,
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
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
    class       => "GATK::Refine",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine",
    option      => "-Xmx40g",
    gatk_option => "--fix_misencoded_quality_scores",
    fasta_file  => $bwa_fasta,
    source_ref  => "bwa",
    vcf_files   => [$dbsnp, $mills],
    gatk_jar    => $gatk_jar,
    picard_jar  => $picard_jar,
    sh_direct   => 0,
    sorted      => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf => {
    class       => "GATK::HaplotypeCallerGVCF",
    perform     => 1,
    target_dir  => "${target_dir}/bwa_refine_hc_gvcf",
    option      => "",
    source_ref  => "bwa_refine",
    java_option => "",
    fasta_file  => $bwa_fasta,
    dbsnp_vcf   => $dbsnp,
    gatk_jar    => $gatk_jar,
    extension   => ".gvcf",
    sh_direct   => 0,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr => {
    class       => "GATK::VariantFilterVQSR",
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
    cqstools    => $cqstools,
    sh_direct   => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=24",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_refine_hc_gvcf_vqsr_annovar => {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/bwa_refine_hc_gvcf_vqsr_annovar",
    source_ref => [ "bwa_refine_hc_gvcf_vqsr", "snp" ],
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => "hg19",
    sh_direct  => 1,
    isvcf      => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => ["fastqc"],
      step2 => ["fastqc_summary"],
      step3 => [ "bwa", "bwa_refine", "bwa_refine_hc_gvcf" ],
      step4 => ["bwa_refine_hc_gvcf_vqsr"],
      step5 => ["bwa_refine_hc_gvcf_vqsr_annovar"],
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

