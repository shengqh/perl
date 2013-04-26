#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/miRNA/VANGARD00055");
my $email      = "quanhu.sheng\@vanderbilt.edu";
my $task_name  = "VANGARD00055";

my $config_rat = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name . "_rat"
  },
  fastqfiles => {
    "2516-01" => ["2516-KCV-1_1_sequence.txt.gz"],
    "2516-02" => ["2516-KCV-2_1_sequence.txt.gz"],
    "2516-03" => ["2516-KCV-3_1_sequence.txt.gz"],
    "2516-04" => ["2516-KCV-4_1_sequence.txt.gz"],
    "2516-05" => ["2516-KCV-5_1_sequence.txt.gz"],
    "2516-06" => ["2516-KCV-6_1_sequence.txt.gz"],
    "2516-07" => ["2516-KCV-7_1_sequence.txt.gz"],
    "2516-08" => ["2516-KCV-8_1_sequence.txt.gz"],
    "2516-09" => ["2516-KCV-9_1_sequence.txt.gz"],
  },
  bwa => {
    target_dir      => "${target_dir}/bwa_genome",
    option          => "-q 15 -t 8",
    option_samse    => "",
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/rn5/rn5.fa",
    estimate_insert => 1,
    source_ref      => "fastqfiles",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
};

my $config_human = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name . "_human"
  },
  fastqfiles => {
    "2516-10" => ["2516-KCV-10_1_sequence.txt.gz"],
    "2516-11" => ["2516-KCV-11_1_sequence.txt.gz"],
    "2516-12" => ["2516-KCV-12_1_sequence.txt.gz"],
    "2516-13" => ["2516-KCV-13_1_sequence.txt.gz"],
    "2516-14" => ["2516-KCV-14_1_sequence.txt.gz"],
    "2516-15" => ["2516-KCV-15_1_sequence.txt.gz"],
    "2516-16" => ["2516-KCV-16_1_sequence.txt.gz"],
    "2516-17" => ["2516-KCV-17_1_sequence.txt.gz"],
    "2516-18" => ["2516-KCV-18_1_sequence.txt.gz"],
    "2516-19" => ["2516-KCV-19_1_sequence.txt.gz"],
    "2516-20" => ["2516-KCV-20_1_sequence.txt.gz"],
    "2516-21" => ["2516-KCV-21_1_sequence.txt.gz"],
    "2516-22" => ["2516-KCV-22_1_sequence.txt.gz"],
    "2516-23" => ["2516-KCV-23_1_sequence.txt.gz"],
    "2516-24" => ["2516-KCV-24_1_sequence.txt.gz"],
    "2516-25" => ["2516-KCV-25_1_sequence.txt.gz"],
    "2516-26" => ["2516-KCV-26_1_sequence.txt.gz"],
    "2516-27" => ["2516-KCV-27_1_sequence.txt.gz"],
    "2516-28" => ["2516-KCV-28_1_sequence.txt.gz"],
    "2516-29" => ["2516-KCV-29_1_sequence.txt.gz"],
    "2516-30" => ["2516-KCV-30_1_sequence.txt.gz"],
    "2516-31" => ["2516-KCV-31_1_sequence.txt.gz"],
    "2516-32" => ["2516-KCV-32_1_sequence.txt.gz"],
    "2516-33" => ["2516-KCV-33_1_sequence.txt.gz"],
    "2516-34" => ["2516-KCV-34_1_sequence.txt.gz"],
    "2516-35" => ["2516-KCV-35_1_sequence.txt.gz"],
    "2516-36" => ["2516-KCV-36_1_sequence.txt.gz"],
    "2516-37" => ["2516-KCV-37_1_sequence.txt.gz"],
    "2516-38" => ["2516-KCV-38_1_sequence.txt.gz"],
    "2516-39" => ["2516-KCV-39_1_sequence.txt.gz"],
    "2516-40" => ["2516-KCV-40_1_sequence.txt.gz"],
    "2516-41" => ["2516-KCV-41_1_sequence.txt.gz"],
    "2516-42" => ["2516-KCV-42_1_sequence.txt.gz"],
    "2516-43" => ["2516-KCV-43_1_sequence.txt.gz"],
    "2516-44" => ["2516-KCV-44_1_sequence.txt.gz"],
    "2516-45" => ["2516-KCV-45_1_sequence.txt.gz"],
    "2516-46" => ["2516-KCV-46_1_sequence.txt.gz"],
    "2516-47" => ["2516-KCV-47_1_sequence.txt.gz"],
    "2516-48" => ["2516-KCV-48_1_sequence.txt.gz"],
  },
  bwa => {
    target_dir      => "${target_dir}/bwa_genome",
    option          => "-q 15 -t 8",
    option_samse    => "",
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/guoy1/reference/hg19/hg19_chr.fa",
    estimate_insert => 1,
    source_ref      => "fastqfiles",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
};

my $config_mirna = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name . "_mirna"
  },
  fastqfiles => {
    "2516-01" => ["2516-KCV-1_1_sequence.txt.gz"],
    "2516-02" => ["2516-KCV-2_1_sequence.txt.gz"],
    "2516-03" => ["2516-KCV-3_1_sequence.txt.gz"],
    "2516-04" => ["2516-KCV-4_1_sequence.txt.gz"],
    "2516-05" => ["2516-KCV-5_1_sequence.txt.gz"],
    "2516-06" => ["2516-KCV-6_1_sequence.txt.gz"],
    "2516-07" => ["2516-KCV-7_1_sequence.txt.gz"],
    "2516-08" => ["2516-KCV-8_1_sequence.txt.gz"],
    "2516-09" => ["2516-KCV-9_1_sequence.txt.gz"],
    "2516-10" => ["2516-KCV-10_1_sequence.txt.gz"],
    "2516-11" => ["2516-KCV-11_1_sequence.txt.gz"],
    "2516-12" => ["2516-KCV-12_1_sequence.txt.gz"],
    "2516-13" => ["2516-KCV-13_1_sequence.txt.gz"],
    "2516-14" => ["2516-KCV-14_1_sequence.txt.gz"],
    "2516-15" => ["2516-KCV-15_1_sequence.txt.gz"],
    "2516-16" => ["2516-KCV-16_1_sequence.txt.gz"],
    "2516-17" => ["2516-KCV-17_1_sequence.txt.gz"],
    "2516-18" => ["2516-KCV-18_1_sequence.txt.gz"],
    "2516-19" => ["2516-KCV-19_1_sequence.txt.gz"],
    "2516-20" => ["2516-KCV-20_1_sequence.txt.gz"],
    "2516-21" => ["2516-KCV-21_1_sequence.txt.gz"],
    "2516-22" => ["2516-KCV-22_1_sequence.txt.gz"],
    "2516-23" => ["2516-KCV-23_1_sequence.txt.gz"],
    "2516-24" => ["2516-KCV-24_1_sequence.txt.gz"],
    "2516-25" => ["2516-KCV-25_1_sequence.txt.gz"],
    "2516-26" => ["2516-KCV-26_1_sequence.txt.gz"],
    "2516-27" => ["2516-KCV-27_1_sequence.txt.gz"],
    "2516-28" => ["2516-KCV-28_1_sequence.txt.gz"],
    "2516-29" => ["2516-KCV-29_1_sequence.txt.gz"],
    "2516-30" => ["2516-KCV-30_1_sequence.txt.gz"],
    "2516-31" => ["2516-KCV-31_1_sequence.txt.gz"],
    "2516-32" => ["2516-KCV-32_1_sequence.txt.gz"],
    "2516-33" => ["2516-KCV-33_1_sequence.txt.gz"],
    "2516-34" => ["2516-KCV-34_1_sequence.txt.gz"],
    "2516-35" => ["2516-KCV-35_1_sequence.txt.gz"],
    "2516-36" => ["2516-KCV-36_1_sequence.txt.gz"],
    "2516-37" => ["2516-KCV-37_1_sequence.txt.gz"],
    "2516-38" => ["2516-KCV-38_1_sequence.txt.gz"],
    "2516-39" => ["2516-KCV-39_1_sequence.txt.gz"],
    "2516-40" => ["2516-KCV-40_1_sequence.txt.gz"],
    "2516-41" => ["2516-KCV-41_1_sequence.txt.gz"],
    "2516-42" => ["2516-KCV-42_1_sequence.txt.gz"],
    "2516-43" => ["2516-KCV-43_1_sequence.txt.gz"],
    "2516-44" => ["2516-KCV-44_1_sequence.txt.gz"],
    "2516-45" => ["2516-KCV-45_1_sequence.txt.gz"],
    "2516-46" => ["2516-KCV-46_1_sequence.txt.gz"],
    "2516-47" => ["2516-KCV-47_1_sequence.txt.gz"],
    "2516-48" => ["2516-KCV-48_1_sequence.txt.gz"],
  },
  bwa => {
    target_dir      => "${target_dir}/bwa_miRBase",
    option          => "-q 15 -t 1",
    option_samse    => "",
    source_ref      => "fastqfiles",
    fasta_file      => "/data/cqs/shengq1/reference/miRBase19/mature.fa",
    estimate_insert => 1,
    source_ref      => "fastqfiles",
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
};

bwa_by_pbs_single( $config_rat, "bwa" );
bwa_by_pbs_single( $config_human, "bwa" );
bwa_by_pbs_single( $config_mirna, "bwa" );

1;
