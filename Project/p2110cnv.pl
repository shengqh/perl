#!/usr/bin/perl
use strict;
use warnings;

use CQS::DNASeq;
use CQS::CNV;
use CQS::FileUtils;
use CQS::SystemUtils;

my $target_dir     = create_directory_or_die("/scratch/cqs/shengq1/dnaseq/2110");
my $probefile      = "/scratch/cqs/lij17/cnv/SureSelect_XT_Human_All_Exon_V4_withoutchr_withoutY_lite.bed";
my $chromosome_dir = "/scratch/cqs/shengq1/references/hg19chromosome";
my $chromosome_len = "/scratch/cqs/shengq1/references/hg19chromosome/hg19.len";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general => {
    path_file => "/home/shengq1/local/bin/path.txt",
    task_name => "2110"
  },
  bamfiles => {
    "2110_JP_01" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-1_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_02" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-2_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_03" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-3_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_04" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-4_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_05" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-5_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_06" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-6_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_07" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-7_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_08" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-8_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_09" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-9_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_10" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-10_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_11" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-11_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_12" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-12_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_13" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-13_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_14" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-14_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_15" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-15_realigned_recal_rmdup.sorted.bam"],
    "2110_JP_16" => ["/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-16_realigned_recal_rmdup.sorted.bam"],
  },
  samtoolsindex => {
    target_dir  => "${target_dir}/samtoolsindex",
    option      => "",
    source_ref  => "bamfiles",
    isbamsorted => 0,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "10gb"
    },
  },
  cnvnator100 => {
    target_dir     => "${target_dir}/cnvnator100",
    option         => "",
    source_ref     => "bamfiles",
    binsize        => 100,
    probefile      => $probefile,
    isbamsorted    => 0,
    chromosome_dir => "/scratch/cqs/shengq1/references/hg19chromosome",
    pbs            => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
  conifer => {
    target_dir  => "${target_dir}/conifer",
    option      => "",
    source_ref  => "bamfiles",
    conifer     => "/home/shengq1/pylibs/bin/conifer.py",
    probefile   => $probefile,
    isbamsorted => 0,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "720",
      "mem"      => "10gb"
    },
  },
  cnmops => {
    target_dir  => "${target_dir}/cnmops",
    option      => "",
    source_ref  => "bamfiles",
    probefile   => $probefile,
    pairmode    => "paired",
    isbamsorted => 1,
    pbs         => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  },
  freec => {
    target_dir             => "${target_dir}/freec",
    option                 => "",
    source_ref             => "bamfiles",
    chrLenFile             => $chromosome_len,
    ploidy                 => 2,
    coefficientOfVariation => 0.01,
    chrFiles               => $chromosome_dir,
    inputFormat            => "BAM",
    mateOrientation        => "FR",
    pbs                    => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "720",
      "mem"      => "40gb"
    },
  }
};

#samtools_index($config, "samtoolsindex");

cnvnator( $config, "cnvnator100" );

conifer( $config, "conifer" );

cnmops( $config, "cnmops" );

freec( $config, "freec" );

1;
