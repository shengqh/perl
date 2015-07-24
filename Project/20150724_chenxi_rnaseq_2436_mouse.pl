#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "20150724_chenxi_rnaseq_2436_mouse";

my $target_dir = "/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse";

#my $target_dir = "e:/temp";

my $transcript_gtf       = "/scratch/cqs/shengq1/references/ensembl_gtf/v81/Mus_musculus.GRCm38.81.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/references/ensembl_gtf/v81/Mus_musculus.GRCm38.81.gtf.index";
my $name_map_file        = "/scratch/cqs/shengq1/references/ensembl_gtf/v81/Mus_musculus.GRCm38.81.map";
my $bowtie2_index        = "/scratch/cqs/shengq1/references/mm10/bowtie2_index_2.2.4/mm10";
my $cqstools             = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email                = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "2436-RM-09" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-9_1_sequence.fastq.gz"],
    "2436-RM-10" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-10_1_sequence.fastq.gz"],
    "2436-RM-11" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-11_1_sequence.fastq.gz"],
    "2436-RM-12" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-12_1_sequence.fastq.gz"],
    "2436-RM-13" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-13_1_sequence.fastq.gz"],
    "2436-RM-14" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-14_1_sequence.fastq.gz"],
    "2436-RM-17" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-17_1_sequence.fastq.gz"],
    "2436-RM-18" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-18_1_sequence.fastq.gz"],
    "2436-RM-19" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-19_1_sequence.fastq.gz"],
    "2436-RM-20" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-20_1_sequence.fastq.gz"],
    "2436-RM-21" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-21_1_sequence.fastq.gz"],
    "2436-RM-22" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-22_1_sequence.fastq.gz"],
  },
  fastqc => {
    class      => "QC::FastQC",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "files",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  fastqc_summary => {
    class      => "QC::FastQCSummary",
    perform    => 1,
    target_dir => "${target_dir}/fastqc",
    option     => "",
    source_ref => "fastqc",
    cqstools   => $cqstools,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class                => "Alignment::Tophat2",
    perform              => 1,
    target_dir           => "${target_dir}/tophat2",
    option               => "--segment-length 25 -r 0 -p 8",
    source_ref           => "files",
    bowtie2_index        => $bowtie2_index,
    transcript_gtf       => $transcript_gtf,
    transcript_gtf_index => $transcript_gtf_index,
    sh_direct            => 0,
    pbs                  => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  tophat2_sort => {
    class         => "Samtools::Sort",
    perform       => 1,
    target_dir    => "${target_dir}/tophat2_sort",
    option        => "",
    source_ref    => [ "tophat2", ".bam" ],
    sort_by_query => 1,
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  tophat2_htseqcount => {
    class      => "Count::HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/tophat2_htseqcount",
    option     => "-r name",
    source_ref => "tophat2_sort",
    gff_file   => $transcript_gtf,
    ispairend  => 1,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tophat2_genetable => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/tophat2_genetable",
    option        => "-p ENS --noheader -e -o ${task}_gene.count",
    source_ref    => "tophat2_htseqcount",
    name_map_file => $name_map_file,
    cqs_tools     => $cqstools,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "tophat2",           "tophat2_sort", "tophat2_htseqcount", "fastqc" ],
      step2 => [ "tophat2_genetable", "fastqc_summary" ],
    },
    sh_direct => 0,
    pbs       => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
};

performConfig($config);

1;

