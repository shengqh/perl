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

my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl_gtf/v81/Mus_musculus.GRCm38.81.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl_gtf/v81/Mus_musculus.GRCm38.81.map";
my $star_index     = "/scratch/cqs/shengq1/references/mm10/STAR_index_v38.81_2.4.2a_sjdb49";
my $fasta_file     = "/scratch/cqs/shengq1/references/mm10/mm10.fa";
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
  star => {
    class      => "Alignment::STAR",
    perform    => 1,
    target_dir => "${target_dir}/star",
    option     => "",
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
  star_index => {
    class          => "Alignment::STARIndex",
    perform        => 1,
    target_dir     => "${target_dir}/star_index",
    option         => "--sjdbOverhang 49 --limitSjdbInsertNsj 2000000",
    source_ref     => [ "star", "tab\$" ],
    fasta_file     => $fasta_file,
    transcript_gtf => $transcript_gtf,
    sh_direct      => 1,
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=24",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  star_2nd_pass => {
    class           => "Alignment::STAR",
    perform         => 1,
    target_dir      => "${target_dir}/star_2nd_pass",
    option          => "",
    source_ref      => "files",
    genome_dir_ref  => "star_index",
    output_unsorted => 1,
    sh_direct       => 0,
    pbs             => {
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
    source_ref => [ "star_2nd_pass", "_Aligned.out.bam" ],
    gff_file   => $transcript_gtf,
    ispairend  => 0,
    isstranded => 1,
    sh_direct  => 1,
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
    option        => "-p ENS --noheader -e -o ${task}_gene.count",
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
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step1 => [ "star",          "fastqc" ],
      step2 => [ "star_index",    "fastqc_summary" ],
      step3 => [ "star_2nd_pass", "star_htseqcount" ],
      step4 => ["star_genetable"],
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

