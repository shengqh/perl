#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20160419_janathan_chipseq_gse53999_hg19");

my $fasta_file   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools     = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "gse53998_hg19";

my $config = {
  general => { task_name => $task },
  files   => {
    "EC_H3K27AC_CON" => ["/scratch/cqs/shengq1/chipseq/20151208_gse53998/cutadapt/result/EC_H3K27AC_CON_clipped.fastq.gz"],
    "EC_WCE_CON"     => ["/scratch/cqs/shengq1/chipseq/20151208_gse53998/cutadapt/result/EC_WCE_CON_clipped.fastq.gz"],
  },
  treatments => {
    "EC_H3K27AC_CON" => ["EC_H3K27AC_CON"],
  },
  inputs => {
    "EC_H3K27AC_CON" => ["EC_WCE_CON"],
  },
  bowtie1 => {
    class         => "Alignment::Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/bowtie1",
    option        => "-v 1 -m 1 --best --strata",
    fasta_file    => $fasta_file,
    source_ref    => "files",
    bowtie1_index => $bowtie_index,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs1callpeak => {
    class      => "Chipseq::MACS",
    perform    => 1,
    target_dir => "${target_dir}/macs1callpeak",
    option     => "-w",
    source_ref => "bowtie1",
    groups_ref => "treatments",
    inputs_ref => "inputs",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2callpeak => {
    class        => "Chipseq::MACS2Callpeak",
    perform      => 1,
    target_dir   => "${target_dir}/macs2callpeak",
    option       => "-g hs --broad -B -p 1e-9",
    source_ref   => "bowtie1",
    groups_ref   => "treatments",
    controls_ref => "inputs",
    sh_direct    => 0,
    pbs          => {
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
      step_1 => ["bowtie1"],
      step_2 => [ "macs1callpeak", "macs2callpeak" ],
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
