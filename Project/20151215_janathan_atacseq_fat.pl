#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/atacseq/20151215_janathan_atacseq_fat");

my $fasta_file   = "/scratch/cqs/shengq1/references/hg18_chrM/bowtie_index_1.1.2/hg18_chrM.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/hg18_chrM/bowtie_index_1.1.2/hg18_chrM";
my $cqstools     = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "fat";

my $config = {
  general => { task_name => $task },
  files   => {
    "SQ1_CHOW"   => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4895_mm9.noChrM.fix.rmdup.sorted.bam"],
    "SQ2_CHOW"   => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4896_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc1_CHOW" => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4897_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc2_CHOW" => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4898_mm9.noChrM.fix.rmdup.sorted.bam"],
    "SQ1_HFD"    => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4899_mm9.noChrM.fix.rmdup.sorted.bam"],
    "SQ2_HFD"    => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4900_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc1_HFD"  => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4901_mm9.noChrM.fix.rmdup.sorted.bam"],
    "Visc2_HFD"  => ["/gpfs21/scratch/cqs/shengq1/atacseq/data/fat/20150911_4902_mm9.noChrM.fix.rmdup.sorted.bam"],
  },
  groups => {
    "SQ_CHOW"   => [ "SQ1_CHOW",   "SQ2_CHOW" ],
    "Visc_CHOW" => [ "Visc1_CHOW", "Visc2_CHOW" ],
    "SQ_HFD"    => [ "SQ1_HFD",    "SQ2_HFD" ],
    "Visc_HFD"  => [ "Visc1_HFD",  "Visc2_HFD" ],
  },
  pairs => {
    "SQ_Visc_CHOW" => [ "SQ_CHOW", "Visc_CHOW" ],
    "SQ_Visc_HFD"  => [ "SQ_HFD",  "Visc_HFD" ],
  },
  bam2bed => {
    class      => "ATACseq::BamToBed",
    perform    => 0,
    target_dir => "${target_dir}/bam2bed",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2callpeak => {
    class      => "Chipseq::MACS2Callpeak",
    perform    => 0,
    target_dir => "${target_dir}/macs2callpeak",
    option     => "-f BED -g mm -B -q 0.01",
    source_ref => "bam2bed",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  macs2bdgdiff => {
    class      => "Chipseq::MACS2Bdgdiff",
    perform    => 1,
    target_dir => "${target_dir}/macs2bdgdiff",
    option     => "",
    source_ref => "macs2callpeak",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
