#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/VANGARD00057");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    bowtie2_index        => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    path_file            => "/home/shengq1/local/bin/path.txt",
    task_name            => "VANGARD00057"
  },
  fastqfiles => {
    "WE5" => [ "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-5_1.fastq.gz", "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-5_2.fastq.gz", ],
    "WE6" => [ "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-6_1.fastq.gz", "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-6_2.fastq.gz", ],
    "WE7" => [ "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-7_1.fastq.gz", "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-7_2.fastq.gz", ],
    "WE8" => [ "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-8_1.fastq.gz", "/blue/sequencer/Runs/projects/2203-WE/2013-03-21/2203-WE-8_2.fastq.gz", ],
  },
  groups => {
    "WE5" => ["WE5"],
    "WE6" => ["WE6"],
    "WE7" => ["WE7"],
    "WE8" => ["WE8"],
  },
  pairs => {
    "WE6_vs_WE5" => [ "WE6", "WE5" ],
    "WE7_vs_WE5" => [ "WE7", "WE5" ],
    "WE8_vs_WE7" => [ "WE8", "WE7" ],
  },
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
  rnaseqc => {
    target_dir     => "${target_dir}/RNASeQC",
    option         => "",
    transcript_gtf => $transcript_gtf,
    genome_fasta   => "/data/cqs/guoy1/reference/hg19/hg19_chr.fa",
    rnaseqc_jar    => "/home/shengq1/local/bin/RNA-SeQC_v1.1.7.jar",
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
    target_dir => "${target_dir}/cufflinks_cuffdiff/result/comparison",
    root_dir   => "${target_dir}/cufflinks_cuffdiff/result",
    gene_only  => 1
  },
};

#fastqc_by_pbs( $config, "fastqc" );

tophat2_by_pbs( $config, "tophat2" );

#call_RNASeQC($config, "rnaseqc");

cufflinks_by_pbs( $config, "cufflinks" );

cuffmerge_by_pbs( $config, "cuffmerge" );

cuffdiff_by_pbs( $config, "cufflinks_cuffdiff" );

#copy_and_rename_cuffdiff_file($config, "rename_diff");

1;
