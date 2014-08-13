#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/VANGARD00054");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    path_file            => "/home/shengq1/local/bin/path.txt",
    task_name            => "VANGARD00054"
  },
  fastqfiles => {
    "DK-1" => ["/scratch/cqs/liuq6/ovariancancer/2503-DK-1_1_sequence.txt.gz"],
    "DK-3" => ["/scratch/cqs/liuq6/ovariancancer/2503-DK-3_1_sequence.txt.gz"],
    "DK-4" => ["/scratch/cqs/liuq6/ovariancancer/2503-DK-4_1_sequence.txt.gz"],
    "DK-5" => ["/scratch/cqs/liuq6/ovariancancer/2503-DK-5_1_sequence.txt.gz"],
    "DK-6" => ["/scratch/cqs/liuq6/ovariancancer/2503-DK-6_1_sequence.txt.gz"],
    "DK-7" => ["/scratch/cqs/liuq6/ovariancancer/2503-DK-7_1_sequence.txt.gz"],
  },
  groups => {
    "CONTROL" => [ "DK-1", "DK-3", "DK-4" ],
    "TREATED" => [ "DK-5", "DK-6", "DK-7" ],
  },
  pairs  => { "TREATED_vs_CONTROL" => [ "TREATED", "CONTROL" ], },
  fastqc => {
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqfiles",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
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
  rnaseqc => {
    class => "CQS::RNASeQC",
    perform => 1,
    target_dir     => "${target_dir}/RNASeQC_2",
    option         => "",
    transcript_gtf => $transcript_gtf,
    fasta_file   => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19.fa",
    jar    => "/home/shengq1/local/bin/RNA-SeQC_v1.1.7.jar",
    source_ref     => "tophat2",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cufflinks => {
    target_dir     => "${target_dir}/cufflinks",
    option         => "-p 8 -u -N",
    transcript_gtf => $transcript_gtf,
    source_ref     => "tophat2",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  cuffmerge => {
    target_dir => "${target_dir}/cuffmerge",
    option     => "-p 8",
    source_ref => "cufflinks",
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  cufflinks_cuffdiff => {
    target_dir         => "${target_dir}/cufflinks_cuffdiff",
    option             => "-p 8 -u -N",
    transcript_gtf_ref => "cuffmerge",
    source_ref         => "tophat2",
    groups_ref         => "groups",
    pairs_ref          => "pairs",
    pbs                => {
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

performTask($config, "rnaseqc");

#fastqc_by_pbs( $config, "fastqc" );

#tophat2_by_pbs( $config, "tophat2" );

#cuffdiff_by_pbs( $config, "cuffdiff" );

##call_RNASeQC($config, "rnaseqc");
#
##cufflinks_by_pbs( $config, "cufflinks" );
#
##cuffmerge_by_pbs( $config, "cuffmerge" );
#
##cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );
#
#copy_and_rename_cuffdiff_file($config, "rename_diff");

1;
