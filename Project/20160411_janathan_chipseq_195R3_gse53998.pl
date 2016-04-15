#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use Data::Dumper;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/chipseq/20160411_janathan_chipseq_195R3_gse53998");

my $fasta_file   = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome.fa";
my $bowtie_index = "/scratch/cqs/shengq1/references/gencode/hg19/bowtie_index_1.1.2/GRCh37.p13.genome";
my $cqstools     = "/home/shengq1/cqstools/cqstools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "chipseq";

my $config = {
  general   => { task_name => $task },
  bam_files => {
    "EC_H3K27AC_CON"                 => ["/scratch/cqs/shengq1/chipseq/20151208_gse53998/bowtie1/result/EC_H3K27AC_CON/EC_H3K27AC_CON.bam"],
    "d17_static_iPSctrl2_H3K27ac" => ["/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/bowtie1/result/d17_static_iPSctrl2_H3K27ac/d17_static_iPSctrl2_H3K27ac.bam"],
  },
  peaks => {
    "EC_H3K27AC_CON"                 => ["/scratch/cqs/shengq1/chipseq/20151208_gse53998/MACS/result/EC_H3K27AC_CON/EC_H3K27AC_CON_peaks.name.bed"],
    "d17_static_iPSctrl2_H3K27ac" => ["/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/macs1callpeak/result/d17_static_iPSctrl2_H3K27ac/d17_static_iPSctrl2_H3K27ac_peaks.name.bed"],
  },
  groups => {
    "comparison" => [ "EC_H3K27AC_CON", "d17_static_iPSctrl2_H3K27ac" ]
  },
  merge_bed => {
    class      => "Bedtools::Merge",
    perform    => 1,
    target_dir => "${target_dir}/merge_bed",
    option     => "-d 10 -c 4 -o collapse -delim \"_UNION_\"",
    source_ref => "peaks",
    groups_ref => "groups",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  htseqcount => {
    class        => "Count::HTSeqCount",
    perform      => 1,
    target_dir   => "${target_dir}/htseqcount",
    option       => "",
    source_ref   => "bam_files",
    gff_file_ref => [ "merge_bed", ".gff\$" ],
    sh_direct    => 0,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  peaktable => {
    class      => "CQS::CQSDatatable",
    perform    => 1,
    target_dir => "${target_dir}/peaktable",
    option     => "--noheader -e -o ${task}_peak.count",
    source_ref => "htseqcount",
    cqs_tools  => $cqstools,
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  correlation => {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/peaktable",
    rtemplate                => "countTableCorrelation.R",
    output_file              => "parameterSampleFile1",
    output_file_ext          => ".Correlation.png",
    parameterSampleFile1_ref => ["peaktable", ".count\$"],
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
      step_1 => [ "merge_bed", "htseqcount" ],
      step_2 => [ "peaktable", "correlation" ],
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

#performTask( $config, "macs2callpeak_bradner_rose2" );

1;
