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
    "EC_BRD4_CON" => ["/gpfs21/scratch/cqs/shengq1/chipseq/20151208_gse53998/bowtie1/result/EC_BRD4_CON/EC_BRD4_CON.bam"],
  },
  treatments => {
    "EC_BRD4_CON" => ["EC_BRD4_CON"],
  },
  "peaks" => {
    "d17_static_iPSctrl2_H3K27ac" => ["/scratch/cqs/shengq1/chipseq/20160302_janathan_chipseq_195R3/macs1callpeak/result/d17_static_iPSctrl2_H3K27ac/d17_static_iPSctrl2_H3K27ac_peaks.name.bed"],
    "EC_BRD4_CON"                 => ["/scratch/cqs/shengq1/chipseq/20160411_janathan_chipseq_195R3_gse53998/macs1callpeak/result/EC_BRD4_CON/EC_BRD4_CON_peaks.name.bed"],
  },
  groups => {
    "comparison" => [ "EC_BRD4_CON", "d17_static_iPSctrl2_H3K27ac" ]
  },
  macs1callpeak => {
    class      => "Chipseq::MACS",
    perform    => 0,
    target_dir => "${target_dir}/macs1callpeak",
    option     => "-p 1e-9 -w -S --space=50",
    source_ref => "bam_files",
    groups_ref => "treatments",
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
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
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 0,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      step_2 => ["macs1callpeak"],
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
