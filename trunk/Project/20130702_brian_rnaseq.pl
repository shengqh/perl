#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/20130702_brian_rnaseq");

my $transcript_gtf       = "/data/cqs/guoy1/reference/annotation2/hg19/Homo_sapiens.GRCh37.68.gtf";
my $transcript_gtf_index = "/scratch/cqs/shengq1/gtfindex/hg19_GRCh37_68";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $config = {
  general    => { task_name => "brian" },
  fastqfiles => {
    "2526-TPS-01" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-1_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-1_2.fastq.gz" ],
    "2526-TPS-02" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-2_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-2_2.fastq.gz" ],
    "2526-TPS-03" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-3_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-3_2.fastq.gz" ],
    "2526-TPS-05" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-5_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-5_2.fastq.gz" ],
    "2526-TPS-06" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-6_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-6_2.fastq.gz" ],
    "2526-TPS-07" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-7_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-7_2.fastq.gz" ],
    "2526-TPS-08" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-8_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-8_2.fastq.gz" ],
    "2526-TPS-09" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-9_1.fastq.gz",  "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-9_2.fastq.gz" ],
    "2526-TPS-10" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-10_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-10_2.fastq.gz" ],
    "2526-TPS-11" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-11_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-11_2.fastq.gz" ],
    "2526-TPS-12" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-12_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-12_2.fastq.gz" ],
    "2526-TPS-13" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-13_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-13_2.fastq.gz" ],
    "2526-TPS-14" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-14_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-14_2.fastq.gz" ],
    "2526-TPS-15" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-15_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-15_2.fastq.gz" ],
    "2526-TPS-16" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-16_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-16_2.fastq.gz" ],
    "2526-TPS-17" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-17_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-17_2.fastq.gz" ],
    "2526-TPS-18" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-18_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-18_2.fastq.gz" ],
    "2526-TPS-20" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-20_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-20_2.fastq.gz" ],
    "2526-TPS-21" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-21_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-21_2.fastq.gz" ],
    "2526-TPS-22" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-22_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-22_2.fastq.gz" ],
    "2526-TPS-23" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-23_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-23_2.fastq.gz" ],
    "2526-TPS-24" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-24_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-24_2.fastq.gz" ],
    "2526-TPS-25" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-25_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-25_2.fastq.gz" ],
    "2526-TPS-26" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-26_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-26_2.fastq.gz" ],
    "2526-TPS-27" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-27_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-27_2.fastq.gz" ],
    "2526-TPS-28" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-28_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-28_2.fastq.gz" ],
    "2526-TPS-29" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-29_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-29_2.fastq.gz" ],
    "2526-TPS-30" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-30_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-30_2.fastq.gz" ],
    "2526-TPS-31" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-31_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-31_2.fastq.gz" ],
    "2526-TPS-32" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-32_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-32_2.fastq.gz" ],
    "2526-TPS-33" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-33_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-33_2.fastq.gz" ],
    "2526-TPS-34" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-34_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-34_2.fastq.gz" ],
    "2526-TPS-35" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-35_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-35_2.fastq.gz" ],
    "2526-TPS-36" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-36_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-36_2.fastq.gz" ],
    "2526-TPS-37" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-37_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-37_2.fastq.gz" ],
    "2526-TPS-38" => [ "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-38_1.fastq.gz", "/gpfs20/data/strickt2/TNBC_RNAseq/run1-38/2526-TPS-38_2.fastq.gz" ],
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
      "nodes"    => "1:ppn=2",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  },
  tophat2 => {
    class         => "Tophat2",
    perform       => 0,
    target_dir    => "${target_dir}/tophat2",
    option        => "--segment-length 25 -r 0 -p 6",
    source_ref    => "fastqfiles",
    bowtie2_index => "/data/cqs/guoy1/reference/hg19/bowtie2_index/hg19",
    sh_direct     => 0,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
};

performConfig($config);

1;
