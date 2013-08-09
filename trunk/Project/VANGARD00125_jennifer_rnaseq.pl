#!/usr/bin/perl
use strict;
use warnings;

use CQS::QC;
use CQS::RNASeq;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ClassFactory;

my $task_name  = "VANGARD00125";
my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/vangard/${task_name}_jennifer_rnaseq");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";
my $bowtie2_index        = "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19";
my $fasta_file           = "/data/cqs/guoy1/reference/hg19/hg19_chr.fa";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => $task_name,
  },
  fastqfiles => {
    "2059-JP-0" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-0_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-0_2.fastq.gz" ],
    "2059-JP-1" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-1_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-1_2.fastq.gz" ],
    "2059-JP-2" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-2_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-2_2.fastq.gz" ],
    "2059-JP-3" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-3_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-3_2.fastq.gz" ],
    "2059-JP-4" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-4_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-4_2.fastq.gz" ],
    "2059-JP-5" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-5_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-5_2.fastq.gz" ],
    "2059-JP-6" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-6_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-6_2.fastq.gz" ],
    "2059-JP-7" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-7_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-7_2.fastq.gz" ],
    "2059-JP-8" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-8_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-8_2.fastq.gz" ],
    "2059-JP-9" => [ "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-9_1.fastq.gz", "/autofs/blue_sequencer/Runs/projects/2059-JP/2013-07-24/2059-JP-9_2.fastq.gz" ],
  },
  fastqc => {
    class      => "FastQC",
    perform    => 0,
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
    perform              => 0,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "fastqfiles",
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    bowtie2_index        => $bowtie2_index,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  reorderbam => {
    class         => "ReorderSam",
    perform       => 0,
    target_dir    => "${target_dir}/reordersam",
    option        => "",
    jar           => "/home/shengq1/local/bin/picard/ReorderSam.jar",
    fasta_file    => $fasta_file,
    source_ref    => "tophat2",
    sort_by_query => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  rnaseqc => {
    class          => "RNASeQC",
    perform        => 0,
    target_dir     => "${target_dir}/rnaseqc",
    option         => "",
    source_ref     => "tophat2",
    jar            => "/home/shengq1/local/bin/RNA-SeQC_v1.1.7.jar",
    fasta_file     => $fasta_file,
    transcript_gtf => $transcript_gtf,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  sortbam => {
    class         => "Sortbam",
    perform       => 1,
    target_dir    => "${target_dir}/sortname",
    option        => "",
    source_ref    => "tophat2",
    sort_by_query => 1,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  htseq => {
    class      => "HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/htseqcount",
    option     => "",
    source_ref => ["sortbam"],
    gff_file   => $transcript_gtf,
    sh_direct  => 1,
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
