#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $task_name  = "bcsubtype";
my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20130708_jennifer_rnaseq_subtype");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name,
  },
  fastqfiles => { "JP2588"           => [ "${target_dir}/data/2588-JP-1_1.fastq.gz", "${target_dir}/data/2588-JP-1_2.fastq.gz", ], },
  groups     => { "JP2588"           => ["JP2588"], },
  pairs      => { "JP2588_vs_JP2588" => [ "JP2588",                                  "JP2588" ], },
  fastqc     => {
    class      => "FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "-p 8",
    source_ref           => "fastqfiles",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  cuffdiff => {
    class          => "Cuffdiff",
    perform        => 1,
    target_dir     => "${target_dir}/cuffdiff",
    option         => "-p 8 -u -N",
    fasta_file     => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa",
    transcript_gtf => $transcript_gtf,
    source_ref     => "tophat2",
    groups_ref     => "groups",
    pairs_ref      => "pairs",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;
