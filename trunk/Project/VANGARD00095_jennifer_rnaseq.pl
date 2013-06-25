#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;
use CQS::ClassFactory;

my $vangard = "VANGARD00095";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${vangard}_jennifer_rnaseq");

my $bwa_dir = "${target_dir}/bwa_refine";

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "${vangard}" },
  fastqfiles => {
    "2562-JP-1" => [ "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-1_1.fastq.gz", "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-1_2.fastq.gz" ],
    "2562-JP-2" => [ "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-2_1.fastq.gz", "/blue/sequencer/Runs/projects/2562-JP/2013-06-20/2562-JP-2_2.fastq.gz" ]
  },
  fastqc => {
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  bwa => {
    class      => "BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa",
    option     => "-q 15 -t 8",
    fasta_file => "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  refine => {
    class              => "GATK::Refine",
    perform            => 1,
    target_dir         => "${target_dir}/refine",
    option             => "-Xmx40g",
    fasta_file         => "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa",
    source_ref         => "bwa",
    thread_count       => 8,
    vcf_files          => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
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
  bwa_refine => {
    target_dir         => $bwa_dir,
    option             => "-q 15 -t 8",
    option_samse       => "",
    option_sampe       => "",
    option_gatk        => "-Xmx40g",
    fasta_file         => "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa",
    source_ref         => "fastqfiles",
    thread_count       => 8,
    vcf_files          => ["/data/cqs/shengq1/reference/snp137/human/00-All.vcf"],
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
  muTect => {
    target_dir    => "${target_dir}/muTect",
    option        => "-nt 8",
    source_ref    => "bwa_refine",
    fasta_file    => "/data/cqs/guoy1/reference/hg19/bwa_index_0.7.4/hg19_chr.fa",
    cosmic_file   => "/data/cqs/shengq1/reference/cosmic/cosmic_v65_28052013.hg19.16571.vcf",
    dbsnp_file    => "/data/cqs/shengq1/reference/snp137/human/00-All.vcf",
    annovar_param => "--buildver hg19 --verdbsnp 137 --ver1000g 1000g2012apr --veresp 6500si --genetype refgene --alltranscript --remove",
    annovar_db    => "/scratch/cqs/shengq1/references/annovar/humandb/",
    sh_direct     => 1,
    execute_file  => "/home/shengq1/local/bin/muTect-1.1.4.jar",
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

#fastqc_by_pbs( $config, "fastqc" );
#bwa_refine( $config, "bwa_refine" );
#muTect( $config, "muTect" );
performConfig($config);

1;
