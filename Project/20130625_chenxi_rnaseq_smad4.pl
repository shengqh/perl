#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20130625_chenxi_rnaseq_smad4");

my $transcript_gtf = "/data/cqs/guoy1/reference/mm10/mm10_annotation/Mus_musculus.GRCm38.68.gtf";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index => "/data/cqs/guoy1/reference/mm10/bowtie2_index/mm10",
    path_file     => "/home/shengq1/bin/path.txt",
    task_name     => "2379"
  },
  fastqfiles => {
    "2288-RDB-35" => ["/data/cqs/chenx/smad4/2013/2288-RDB-35_1_sequence.txt.gz"],
    "2288-RDB-39" => ["/data/cqs/chenx/smad4/2013/2288-RDB-39_1_sequence.txt.gz"],
    "2288-RDB-40" => ["/data/cqs/chenx/smad4/2013/2288-RDB-40_1_sequence.txt.gz"],
    "2288-RDB-47" => ["/data/cqs/chenx/smad4/2013/2288-RDB-47_1_sequence.txt.gz"],
    "2288-RDB-48" => ["/data/cqs/chenx/smad4/2013/2288-RDB-48_1_sequence.txt.gz"],
    "2288-RDB-51" => ["/data/cqs/chenx/smad4/2013/2288-RDB-51_1_sequence.txt.gz"],
  },
  groups => {
    "WT" => [ "2288-RDB-35", "2288-RDB-47", "2288-RDB-48" ],
    "KO" => [ "2288-RDB-39", "2288-RDB-40", "2288-RDB-51" ],
  },
  pairs  => { "KO_vs_WT" => [ "WT", "KO" ], },
  fastqc => {
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class      => "Tophat2",
    perform    => 1,
    target_dir => "${target_dir}/tophat2",
    option     => "--segment-length 25 -r 0 -p 6",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  cuffdiff => {
    class                => "Cuffdiff",
    perform              => 1,
    target_dir           => "${target_dir}/cuffdiff",
    option               => "-p 8 -u -N",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => "/scratch/cqs/shengq1/gtfindex/mm10_GRCm38_68",
    source_ref           => "tophat2",
    groups_ref           => "groups",
    pairs_ref            => "pairs",
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
};

fastqc_by_pbs( $config, "fastqc" );
performConfig($config);

1;
