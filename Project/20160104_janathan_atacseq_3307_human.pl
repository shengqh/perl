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

my $cqstools  = "/home/shengq1/cqstools/CQS.Tools.exe";
my $bwa_fasta = "/scratch/cqs/shengq1/references/hg19_tmp/bwa_index_0.7.12/hg19_16569_M.fa";

my $config = {
  general => { task_name => "3307" },
  files   => {
    "JDB1_1K_NoTNF"  => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-1_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-1_2_sequence.txt.gz" ],
    "JDB2_50K_NoTNF" => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-2_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-2_2_sequence.txt.gz" ],
    "JDB3_1K_TNF"    => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-3_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-3_2_sequence.txt.gz" ],
    "JDB4_50K_TNF"   => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-4_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-4_2_sequence.txt.gz" ],
    "JDB6_1K_NoTNF"  => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-6_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-6_2_sequence.txt.gz" ],
    "JDB7_1K_TNF"    => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-7_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-7_2_sequence.txt.gz" ],
    "JDB8_50K_NoTNF" => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-8_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-8_2_sequence.txt.gz" ],
    "JDB9_50K_TNF"   => [ "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-9_1_sequence.txt.gz", "/data/cqs/shengq1/data/3307/PE-150/3307-JDB-9_2_sequence.txt.gz" ],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 0,
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
    perform    => 0,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    cqstools   => $cqstools,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  cutadapt => {
    class      => "Trimmer::Cutadapt",
    perform    => 0,
    target_dir => "${target_dir}/cutadapt",
    option     => "-O 10 -m 50",
    source_ref => "files",
    adapter    => "CTGTCTCTTATA",
    extension  => "_clipped.fastq.gz",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  fastqlen => {
    class      => "CQS::FastqLen",
    perform    => 0,
    target_dir => "${target_dir}/fastqlen",
    option     => "",
    source_ref => "cutadapt",
    cqstools   => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "24",
      "mem"      => "20gb"
    },
  },
  bwa_pretrim => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa_pretrim",
    option     => "",
    bwa_index  => $bwa_fasta,
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  bwa_posttrim => {
    class      => "Alignment::BWA",
    perform    => 1,
    target_dir => "${target_dir}/bwa_posttrim",
    option     => "",
    bwa_index  => $bwa_fasta,
    source_ref => "cutadapt",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

};

performConfig($config);

1;

