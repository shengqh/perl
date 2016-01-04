#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Hash::Merge qw( merge );

my $target_dir = "/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human";
my $email      = "quanhu.sheng\@vanderbilt.edu";

my $cqstools       = "/home/shengq1/cqstools/CQS.Tools.exe";
my $bwa_fasta      = "/scratch/cqs/shengq1/references/hg19_16569_MT/bwa_index_0.7.12/hg19_16569_MT.fa";
my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v75/Homo_sapiens.GRCh37.75.MT.map";

my $dbsnp      = "/scratch/cqs/shengq1/references/dbsnp/human_GRCh37_v142_16569_MT.vcf";
my $hapmap     = "/scratch/cqs/shengq1/references/gatk/b37/hapmap_3.3.b37.vcf";
my $omni       = "/scratch/cqs/shengq1/references/gatk/b37/1000G_omni2.5.b37.vcf";
my $g1000      = "/scratch/cqs/shengq1/references/gatk/b37/1000G_phase1.snps.high_confidence.b37.vcf";
my $mills      = "/scratch/cqs/shengq1/references/gatk/b37/Mills_and_1000G_gold_standard.indels.b37.vcf";
my $cosmic     = "/scratch/cqs/shengq1/references/cosmic/cosmic_v71_hg19_16569_M.vcf";
my $gatk_jar   = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_jar = "/scratch/cqs/shengq1/local/bin/picard/picard.jar";
my $qc3_perl   = "/scratch/cqs/shengq1/local/bin/qc3/qc3.pl";

my $config = {
  general => { task_name => "3307" },
  files   => {
    "JDB1_1K_NoTNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-1_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-1_2_sequence.txt.gz"
    ],
    "JDB2_50K_NoTNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-2_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-2_2_sequence.txt.gz"
    ],
    "JDB3_1K_TNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-3_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-3_2_sequence.txt.gz"
    ],
    "JDB4_50K_TNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-4_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-4_2_sequence.txt.gz"
    ],
    "JDB6_1K_NoTNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-6_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-6_2_sequence.txt.gz"
    ],
    "JDB7_1K_TNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-7_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-7_2_sequence.txt.gz"
    ],
    "JDB8_50K_NoTNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-8_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-8_2_sequence.txt.gz"
    ],
    "JDB9_50K_TNF" => [
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-9_1_sequence.txt.gz",
      "/gpfs21/scratch/cqs/shengq1/atacseq/20160104_janathan_atacseq_3307_human/data/PE-150/3307-JDB-9_2_sequence.txt.gz"
    ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "40gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  }
};

performConfig($config);

1;

