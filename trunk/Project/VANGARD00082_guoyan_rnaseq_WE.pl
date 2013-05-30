#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $task_name = "VANGARD00082";
my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${task_name}_guoyan_rnaseq_WE");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    path_file            => "/home/shengq1/local/bin/path.txt",
    task_name            => $task_name,
  },
  fastqfiles => {
    "WE3" => [ "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-3_1.fastq.gz", "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-3_2.fastq.gz", ],
    "WE4" => [ "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-4_1.fastq.gz", "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-4_2.fastq.gz", ],
    "WE5" => [ "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-5_1.fastq.gz", "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-5_2.fastq.gz", ],
    "WE6" => [ "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-6_1.fastq.gz", "/blue/sequencer/Runs/projects/2461-WE/2013-04-16/2461-WE-6_2.fastq.gz", ],
  },
  groups => {
    "WE3" => ["WE3"],
    "WE4" => ["WE4"],
    "WE5" => ["WE5"],
    "WE6" => ["WE6"],
  },
  pairs => {
    "WE4_vs_WE3" => [ "WE3", "WE4" ],
    "WE5_vs_WE6" => [ "WE6", "WE5" ],
  },
  tophat2 => {
    target_dir           => "${target_dir}/tophat2",
    option               => "-p 8",
    batchmode            => 0,
    sortbam              => 1,
    source_ref           => "fastqfiles",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "240",
      "mem"      => "40gb"
    },
  },
  cuffdiff => {
    target_dir     => "${target_dir}/cuffdiff",
    option         => "-p 8 -u -N",
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
  rename_diff => {
    target_dir => "${target_dir}/cuffdiff/result/comparison",
    root_dir   => "${target_dir}/cuffdiff/result",
    gene_only  => 1
  },
};

tophat2_by_pbs( $config, "tophat2" );

cuffdiff_by_pbs( $config, "cuffdiff" );

#copy_and_rename_cuffdiff_file( $config, "rename_diff" );

1;
