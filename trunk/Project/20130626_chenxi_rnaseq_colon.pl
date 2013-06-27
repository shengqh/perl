#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20130626_chenxi_rnaseq_colon");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "colon" },
  fastqfiles => {
    "2280-RDB-22" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-22_1_sequence.txt.gz"],
    "2280-RDB-23" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-23_1_sequence.txt.gz"],
    "2280-RDB-24" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-24_1_sequence.txt.gz"],
    "2280-RDB-25" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-25_1_sequence.txt.gz"],
    "2280-RDB-26" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-26_1_sequence.txt.gz"],
    "2280-RDB-27" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-27_1_sequence.txt.gz"],
    "2280-RDB-28" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-28_1_sequence.txt.gz"],
    "2280-RDB-29" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-29_1_sequence.txt.gz"],
    "2280-RDB-30" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-30_1_sequence.txt.gz"],
    "2280-RDB-31" => ["/gpfs20/data/cqs/chenx/colon/06_2013/2280-RDB-31_1_sequence.txt.gz"],
  },
  groups => {
    "S1" => [ "2280-RDB-22", "2280-RDB-23", "2280-RDB-24", "2280-RDB-25", "2280-RDB-26" ],
    "S2" => [ "2280-RDB-27", "2280-RDB-28", "2280-RDB-29", "2280-RDB-30", "2280-RDB-31" ],
  },
  pairs  => { "S2_vs_S1" => [ "S1", "S2" ], },
  fastqc => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class         => "Tophat2",
    perform       => 1,
    target_dir    => "${target_dir}/tophat2",
    option        => "--segment-length 25 -r 0 -p 6",
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    sh_direct     => 1,
    pbs           => {
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
    transcript_gtf_index => $transcript_gtf_index,
    fasta_file           => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa",
    source_ref           => "tophat2",
    groups_ref           => "groups",
    pairs_ref            => "pairs",
    sh_direct            => 1,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
};

performConfig($config);

1;
