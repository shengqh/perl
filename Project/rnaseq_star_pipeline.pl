#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "star_pipeline";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/star_pipeline");

#my $target_dir = "e:/temp";

my $transcript_gtf = "/scratch/cqs/shengq1/references/ucsc/mm10_ucsc.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ucsc/mm10_ucsc.map";
my $star_index     = "/scratch/cqs/shengq1/references/mm10_chr/STAR_index_ucsc_sjdb49";
my $fasta_file     = "/scratch/cqs/shengq1/references/mm10_chr/mm10.fa";
my $cqstools       = "/home/shengq1/cqstools/CQS.Tools.exe";
my $email          = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "2436-RM-09" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-9_1_sequence.fastq.gz"],
    "2436-RM-10" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-10_1_sequence.fastq.gz"],
    "2436-RM-11" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-11_1_sequence.fastq.gz"],
    "2436-RM-12" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-12_1_sequence.fastq.gz"],
    "2436-RM-13" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-13_1_sequence.fastq.gz"],
    "2436-RM-14" => ["/gpfs21/scratch/cqs/shengq1/rnaseq/20150724_chenxi_rnaseq_mouse/raw/2436-RM-14_1_sequence.fastq.gz"],
  },
  groups => {
    "G1" => [ "2436-RM-09", "2436-RM-10", "2436-RM-11" ],
    "G2" => [ "2436-RM-12", "2436-RM-13", "2436-RM-14" ],
  },
  pairs => {
    "G2_VS_G1_Paired" => {
      groups => [ "G1", "G2" ],
      paired => [ "P1", "P2", "P3" ]
    },
    "G2_VS_G1_UnPaired" => { groups => [ "G1", "G2" ], },
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
  star => {
    class      => "Alignment::STAR",
    perform    => 1,
    target_dir => "${target_dir}/star",
    option     => "--twopassMode Basic",
    source_ref => "files",
    genome_dir => $star_index,
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  star_htseqcount => {
    class      => "Count::HTSeqCount",
    perform    => 1,
    target_dir => "${target_dir}/star_htseqcount",
    option     => "",
    source_ref => [ "star", "_Aligned.out.bam" ],
    gff_file   => $transcript_gtf,
    ispairend  => 0,
    stranded   => "reverse",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  star_genetable => {
    class         => "CQS::CQSDatatable",
    perform       => 1,
    target_dir    => "${target_dir}/star_genetable",
    option        => "-p uc --noheader -e -o ${task}_gene.count",
    source_ref    => "star_htseqcount",
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
  star_deseq2 => {
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/star_deseq2",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile_ref        => "star_genetable",
    sh_direct            => 1,
    show_DE_gene_cluster => 1,
    pvalue               => 0.05,
    fold_change          => 2.0,
    min_median_read      => 5,
    pbs                  => {
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
      step1 => [ "fastqc",         "star",           "star_htseqcount" ],
      step2 => [ "fastqc_summary", "star_genetable", "star_deseq2" ],
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

