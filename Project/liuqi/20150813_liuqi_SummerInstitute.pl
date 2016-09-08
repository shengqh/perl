#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $task_name  = "Qi";
my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20150813_liuqi_SummerInstitute");

my $transcript_gtf       = "/scratch/liuq6/reference/Homo_sapiens.GRCh37.75_chr1-22-X-Y-M.gtf";
my $transcript_gtf_index = "/scratch/liuq6/reference/gtfindex/Homo_sapiens.GRCh37.75";
my $bowtie2_index        = "/scratch/liuq6/reference/bowtie2_index/hg19";
my $fasta_file           = "/scratch/liuq6/reference/bowtie2_index/hg19.fa";

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task_name, },
  files   => {
    "BRD4_ctrl1" => ["/gpfs21/scratch/liuq6/SummerInstitute/data/BRD4_ctrl1.fastq"],
    "BRD4_ctrl2" => ["/gpfs21/scratch/liuq6/SummerInstitute/data/BRD4_ctrl2.fastq"],
    "BRD4_ER1"   => ["/gpfs21/scratch/liuq6/SummerInstitute/data/BRD4_ER1.fastq"],
    "BRD4_ER2"   => ["/gpfs21/scratch/liuq6/SummerInstitute/data/BRD4_ER2.fastq"],
    "ctrl1"      => ["/gpfs21/scratch/liuq6/SummerInstitute/data/ctrl1.fastq"],
    "ctrl2"      => ["/gpfs21/scratch/liuq6/SummerInstitute/data/ctrl2.fastq"],
    "ER1"        => ["/gpfs21/scratch/liuq6/SummerInstitute/data/ER1.fastq"],
    "ER2"        => ["/gpfs21/scratch/liuq6/SummerInstitute/data/ER2.fastq"],
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
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Alignment::Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "-p 8 --transcriptome-only",
    source_ref           => "files",
    transcript_gtf_index => $transcript_gtf_index,
    bowtie2_index        => $bowtie2_index,
    rename_bam           => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
