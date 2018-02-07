#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq2/test");
my $email      = "quanhu.sheng.1\@vanderbilt.edu";

my $config = {
  general => { task_name => "mutect2indel" },
  files   => {
    "tumor"  => ["/data/cqs/liuq6/Beauchamp/exomeseq_PXD/02somatic/bwa_refine/result/WD50232_P1.rmdup.indel.recal.bam.PTEN.bam"],
    "normal" => ["/data/cqs/liuq6/Beauchamp/exomeseq_PXD/02somatic/bwa_refine/result/WD50232_Normal.rmdup.indel.recal.bam.PTEN.bam"],
  },
  mutect2indel => {
    class      => "GATK::MuTect2Indel",
    perform    => 1,
    target_dir => "${target_dir}/mutect2Indel",
    option     => "-maxNumHaplotypesInPopulation 3",
    gatk_jar   => "/scratch/cqs/shengq2/local/bin/gatk/GenomeAnalysisTK.jar",
    fasta_file => "/scratch/cqs/shengq2/references/gatk/b37/bwa_index_0.7.12/human_g1k_v37.fasta",
    dbsnp_file => "/scratch/cqs/shengq2/references/gatk/b37/dbsnp_150.b37.vcf",
    source_ref => "files",
    groups     => {
      "TEST" => [ "normal", "tumor" ],
    },
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
