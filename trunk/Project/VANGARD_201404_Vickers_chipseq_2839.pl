#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;
use CQS::CNV;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201404_chipseq_2839");

my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $fasta_file = "/data/cqs/shengq1/reference/hg19.16569/bwa_index_0.7.4/hg19_rCRS.fa";
my $dbsnp      = "/data/cqs/guoy1/reference/dbsnp138/00-All.vcf";
my $gatk       = "/home/shengq1/local/bin/GATK/GenomeAnalysisTK.jar";
my $picard_dir = "/home/shengq1/local/bin/picard/";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "201404_chipseq" },
  fastqfiles => {
    "2839-KCV-1" => ["/gpfs21/scratch/vantage_repo/Vickers/2839/2839-KCV-1_1_sequence.txt.gz"],
    "2839-KCV-2" => ["/gpfs21/scratch/vantage_repo/Vickers/2839/2839-KCV-2_1_sequence.txt.gz"],
    "2839-KCV-3" => ["/gpfs21/scratch/vantage_repo/Vickers/2839/2839-KCV-3_1_sequence.txt.gz"],
  },
  groups => {
    "2839-KCV-1" => ["2839-KCV-1"],
    "2839-KCV-2" => ["2839-KCV-2"],
    "Control"    => ["2839-KCV-3"],
  },
  pairs => {
    "2839-KCV-1_vs_Control" => [ "Control", "2839-KCV-1" ],
    "2839-KCV-2_vs_Control" => [ "Control", "2839-KCV-2" ],
  },
  fastq_trimmer => {
    class      => "CQS::FastqTrimmer",
    perform    => 1,
    target_dir => "${target_dir}/FastqTrimmer",
    option     => "-n -z",
    source_ref => "fastqfiles",
    extension  => "_trimmed.fastq.gz",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastq_trimmer",
    sh_direct  => 1,
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
    option     => "-t 8",
    fasta_file => $fasta_file,
    source_ref => "fastq_trimmer",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  filter => {
    class      => "Samtools::View",
    perform    => 1,
    target_dir => "${target_dir}/filter",
    option     => "-b -q 20",
    source_ref => "bwa",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  markdup => {
    class              => "Picard::MarkDuplicates",
    perform            => 1,
    target_dir         => "${target_dir}/markdup",
    option             => "-Xmx40g",
    source_ref         => "filter",
    thread_count       => 8,
    markDuplicates_jar => "${picard_dir}/MarkDuplicates.jar",
    sh_direct          => 0,
    pbs                => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  MACS => {
    class      => "Chipseq::MACS",
    perform    => 1,
    target_dir => "${target_dir}/MACS",
    option     => "",
    source_ref => "markdup",
    groups_ref => "groups",
    pairs_ref  => "pairs",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
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
      individual => [ "fastq_trimmer", "fastqc", "bwa", "filter", "markdup" ],
      pair       => ["MACS"]
    },
    sh_direct => 1,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

#performTask( $config, "homerMakeTagDirectory" );
#performTask( $config, "homerFindPeaks" );

1;
