#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::DNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::SomaticMutation;
use CQS::ClassFactory;

my $vangard = "VANGARD00223";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${vangard}_liuqi_rnaseq_bacteria");

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "${vangard}" },
  fastqfiles => {
    "2732-EPS-01" => ["${target_dir}/raw/2732-EPS-1_1_sequence.txt"],
    "2732-EPS-02" => ["${target_dir}/raw/2732-EPS-2_1_sequence.txt"],
    "2732-EPS-03" => ["${target_dir}/raw/2732-EPS-3_1_sequence.txt"],
    "2732-EPS-04" => ["${target_dir}/raw/2732-EPS-4_1_sequence.txt"],
    "2732-EPS-05" => ["${target_dir}/raw/2732-EPS-5_1_sequence.txt"],
    "2732-EPS-06" => ["${target_dir}/raw/2732-EPS-6_1_sequence.txt"],
    "2732-EPS-07" => ["${target_dir}/raw/2732-EPS-7_1_sequence.txt"],
    "2732-EPS-08" => ["${target_dir}/raw/2732-EPS-8_1_sequence.txt"],
    "2732-EPS-09" => ["${target_dir}/raw/2732-EPS-9_1_sequence.txt"],
    "2732-EPS-10" => ["${target_dir}/raw/2732-EPS-10_1_sequence.txt"],
    "2732-EPS-11" => ["${target_dir}/raw/2732-EPS-11_1_sequence.txt"],
    "2732-EPS-12" => ["${target_dir}/raw/2732-EPS-12_1_sequence.txt"],
  },
  groups => {
    "DMSO"                 => [ "2732-EPS-01", "2732-EPS-02", "2732-EPS-03" ],
    "50uM_0070"            => [ "2732-EPS-04", "2732-EPS-05", "2732-EPS-06" ],
    "20uM_NDGA"            => [ "2732-EPS-07", "2732-EPS-08", "2732-EPS-09" ],
    "100uM_chlorpromazine" => [ "2732-EPS-10", "2732-EPS-11", "2732-EPS-12" ],
  },
  pairs => {
    "50uM_0070_vs_DMSO"            => [ "50uM_0070",            "DMSO" ],
    "20uM_NDGA_vs_DMSO"            => [ "20uM_NDGA",            "DMSO" ],
    "100uM_chlorpromazine_vs_DMSO" => [ "100uM_chlorpromazine", "DMSO" ]
  },
  rockhopper => {
    class          => "Bacteria::RockHopper",
    perform        => 1,
    target_dir     => "${target_dir}/rockhopper",
    source_ref     => "fastqfiles",
    groups_ref     => "groups",
    pairs_ref      => "pairs",
    java_option    => "-Xmx10g",
    rockhopper_jar => "/scratch/cqs/shengq1/local/bin/Rockhopper.jar",
    genome_dir     => "${target_dir}/genome",
    option         => "-p 8 -TIME -v true",
    sh_direct      => 1,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  bowtie2 => {
    class          => "Bowtie2",
    perform        => 1,
    target_dir     => "${target_dir}/bowtie2",
    source_ref     => "fastqfiles",
    bowtie2_index =>"${target_dir}/genome/bowtie2-index",
    option         => "-p 8 --very-sensitive",
    sh_direct      => 1,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
};

performConfig($config);

1;
