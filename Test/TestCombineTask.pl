#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/temp/test");

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => "testCombineTask" },
  files   => {
    "MiSeqSample1" => [
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20140218_bojana_MiSeq_HiSeq/rawdata/IG-33_408637_S4_L001_R1_001.fastq.gz",
      "/gpfs21/scratch/cqs/shengq1/rnaseq/20140218_bojana_MiSeq_HiSeq/rawdata/IG-33_408637_S4_L001_R2_001.fastq.gz"
    ],
  },
  preprocessing => {
    class      => "CQS::CombineTask",
    source     => [ "preprocessing::remove_contamination_sequences", "preprocessing::cutadapt" ],
    target_dir => $target_dir,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "20gb"
    },
    "remove_contamination_sequences" => {
      class      => "CQS::Perl",
      perform    => 1,
      target_dir => "${target_dir}/remove_contamination_sequences",
      option     => "ABCDEFG",
      output_ext => "_removeSeq.fastq.gz",
      perlFile   => "removeSequenceInFastq.pl",
      source_ref => "files",
      sh_direct  => 0,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "2",
        "mem"      => "20gb"
      },
    },
    cutadapt => {
      class      => "Trimming::Cutadapt",
      perform    => 1,
      target_dir => "${target_dir}/cutadapt",
      option     => "-m 12 -O 10 -e 0.083",
      source_ref => "preprocessing::remove_contamination_sequences",
      adaptor    => "TGGAATTCTCGGGTGCCAAGG",
      extension  => "_clipped.fastq",
      sh_direct  => 0,
      gzipped    => 1,
      pbs        => {
        "email"    => $email,
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "20gb"
      },
    },
  }
};

performConfig($config);

1;
