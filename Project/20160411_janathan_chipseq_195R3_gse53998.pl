#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20160411_janathan_chipseq_195R3_gse53998");

my $fasta_file     = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools       = "/home/shengq1/cqstools/cqstools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "Comparison";

my $config = {
  general   => { task_name => $task },
  bam_files => {
    "EC_BRD4_CON" => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/bowtie1/result/EC_BRD4_CON/EC_BRD4_CON.bam"],
  },
  treatments => {
    "EC_BRD4_CON" => ["EC_BRD4_CON"],
  },
  macs1callpeak => {
    class      => "Chipseq::MACS",
    perform    => 1,
    target_dir => "${target_dir}/macs1callpeak",
    option     => "-p 1e-9 -w -S --space=20",
    source_ref => "bam_files",
    groups_ref => "treatments",
    sh_direct  => 0,
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
      step_2 => ["macs1callpeak"],
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

#performTask( $config, "macs2callpeak_bradner_rose2" );

1;
