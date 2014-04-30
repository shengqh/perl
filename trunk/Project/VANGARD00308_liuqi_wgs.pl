#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $vangard = "VANGARD00308";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${vangard}_liuqi_wgs");

my $fasta_file = "/data/cqs/shengq1/reference/hg19_illumina/hg19_illumina_16571_chr.fa";
my $chrLenFile = "/data/cqs/shengq1/reference/hg19_illumina/hg19_illumina_16571_chr.len";
my $chrFiles   = "/data/cqs/shengq1/reference/hg19_illumina/chromosomes";

my $mutect = "/home/shengq1/local/bin/muTect-1.1.4.jar";
my $cosmic = "/data/cqs/shengq1/reference/hg19_illumina/illumina-cosmicv67.vcf";
my $dbsnp  = "/data/cqs/shengq1/reference/hg19_illumina/illumina-dbsnp138.vcf";

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $annovar_param = "-protocol refGene,snp137,cosmic64,esp6500si_all,1000g2012apr_all -operation g,f,f,f,f --remove";
my $annovar_db    = "/scratch/cqs/shengq1/references/annovar/humandb/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => "${vangard}" },
  files   => {
    "LP6005748-NT" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00308_liuqi_wgs/bam/LP6005748-DNA_A01/Assembly/LP6005748-DNA_A01.bam"],
    "LP6005749-TP" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00308_liuqi_wgs/bam/LP6005749-DNA_A01/Assembly/LP6005749-DNA_A01.bam"],
  },
  groups => { "VANGARD00308" => [ "LP6005748-NT", "LP6005749-TP" ], },
  muTect => {
    class       => "GATK::MuTect",
    perform     => 1,
    target_dir  => "${target_dir}/muTect",
    option      => "",
    source_ref  => "files",
    groups_ref  => "groups",
    java_option => "-Xmx40g",
    fasta_file  => $fasta_file,
    cosmic_file => $cosmic,
    dbsnp_file  => $dbsnp,
    sh_direct   => 0,
    muTect_jar  => $mutect,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
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
  freec => {
    class                  => "CNV::Freec",
    perform                => 1,
    target_dir             => "${target_dir}/freec",
    option                 => "",
    source_ref             => "files",
    groups_ref             => "groups",
    chrLenFile             => "",
    source_type            => "bam",                   #source_type can be bam/mpileup
    chrLenFile             => $chrLenFile,
    chrFiles               => $chrFiles,
    ploidy                 => 2,
    coefficientOfVariation => 0.05,
    inputFormat            => "bam",
    mateOrientation        => "FR",
    maxThreads             => "8",
    sh_direct              => 0,
    pbs                    => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
