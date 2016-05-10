#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;

my $task = "B3436";

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20160509_brown_3436");

#my $target_dir = "e:/temp";

my $transcript_gtf = "/scratch/cqs/shengq1/references/ensembl/v78/Mus_musculus.GRCm38.78.MT.gtf";
my $name_map_file  = "/scratch/cqs/shengq1/references/ensembl/v78/Mus_musculus.GRCm38.78.MT.map";
my $star_index     = "/scratch/cqs/shengq1/references/mm10_sorted_M/STAR_index_v38.78_2.5.0b_sjdb49";
my $fasta_file     = "/scratch/cqs/shengq1/references/mm10_sorted_M/mm10.fa";
my $cqstools       = "/home/shengq1/cqstools/cqstools.exe";
my $email          = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => { task_name => $task },
  files   => {
    "JDB-01" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-1_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-1_2_sequence.txt.gz" ],
    "JDB-02" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-2_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-2_2_sequence.txt.gz" ],
    "JDB-03" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-3_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-3_2_sequence.txt.gz" ],
    "JDB-04" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-4_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-4_2_sequence.txt.gz" ],
    "JDB-05" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-5_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-5_2_sequence.txt.gz" ],
    "JDB-06" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-6_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-6_2_sequence.txt.gz" ],
    "JDB-07" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-7_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-7_2_sequence.txt.gz" ],
    "JDB-08" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-8_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-8_2_sequence.txt.gz" ],
    "JDB-09" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-9_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-9_2_sequence.txt.gz" ],
    "JDB-11" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-11_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-11_2_sequence.txt.gz" ],
    "JDB-12" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-12_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-12_2_sequence.txt.gz" ],
    "JDB-13" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-13_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-13_2_sequence.txt.gz" ],
    "JDB-14" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-14_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-14_2_sequence.txt.gz" ],
    "JDB-16" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-16_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-16_2_sequence.txt.gz" ],
    "JDB-17" =>
      [ "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-17_1_sequence.txt.gz", "/gpfs21/scratch/cqs/shengq1/rnaseq/20160509_brown_3436/data/3436-JDB-17_2_sequence.txt.gz" ],
  },
  groups => {
    "FED"      => [ "JDB-01", "JDB-02", "JDB-03" ],
    "DMSO"     => [ "JDB-04", "JDB-05", "JDB-06" ],
    "JQ1"      => [ "JDB-07", "JDB-08", "JDB-09" ],
    "CAPTISOL" => [ "JDB-11", "JDB-12", "JDB-16" ],
    "dBET"     => [ "JDB-13", "JDB-14", "JDB-17" ],
  },
  pairs => {
    "DMSO_vs_FED"      => [ "FED",      "DMSO" ],
    "JQ1_vs_FED"       => [ "FED",      "JQ1" ],
    "CAPTISOL_vs_FED"  => [ "FED",      "CAPTISOL" ],
    "dBET_vs_FED"      => [ "FED",      "dBET" ],
    "JQ1_vs_DMSO"      => [ "DMSO",     "JQ1" ],
    "dBET_vs_CAPTISOL" => [ "CAPTISOL", "dBET" ],
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
    class                     => "Alignment::STAR",
    perform                   => 1,
    target_dir                => "${target_dir}/star",
    option                    => "--twopassMode Basic",
    source_ref                => "files",
    genome_dir                => $star_index,
    output_sort_by_coordinate => 1,
    sh_direct                 => 0,
    pbs                       => {
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
    ispairend  => 1,
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
  star_genetable_deseq2 => {
    class                => "Comparison::DESeq2",
    perform              => 1,
    target_dir           => "${target_dir}/star_genetable_deseq2",
    option               => "",
    source_ref           => "pairs",
    groups_ref           => "groups",
    countfile_ref        => "star_genetable",
    sh_direct            => 1,
    show_DE_gene_cluster => 0,
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
  star_genetable_correlation => {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/star_genetable_correlation",
    rtemplate                => "countTableVisFunctions.R,countTableCorrelation.R",
    output_file              => "parameterSampleFile1",
    output_file_ext          => ".Correlation.png",
    parameterSampleFile1_ref => [ "star_genetable", ".count\$" ],
    parameterSampleFile2_ref => "groups",
    sh_direct                => 1,
    pbs                      => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
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
      step2 => [ "fastqc_summary", "star_genetable", "star_genetable_correlation", "star_genetable_deseq2" ],
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

