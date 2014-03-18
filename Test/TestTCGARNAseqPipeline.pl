#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;

my $target_dir = create_directory_or_die("/scratch/cqs/shengq1/rnaseq/TCGApipeline");

my $picard_dir   = "/home/shengq1/local/bin/picard";
my $bedtools_dir = "/scratch/cqs/lij17/softwares/bin/bedtools-2.17.0/bin";
my $tcga_bin_dir = "/scratch/cqs/shengq1/local/bin/TCGA/webshare.bioinf.unc.edu/public/mRNAseq_TCGA";

my $email = "quanhu.sheng\@vanderbilt.edu";
my $task  = "pipeline";

my $config = {
  general => { task_name => $task },
  files   => {
    "S1" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s1_sequence.txt"],
    "S2" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s2_sequence.txt"],
    "S3" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s3_sequence.txt"],
    "S4" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s4_sequence.txt"],
    "S5" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s5_sequence.txt"],
    "S6" => ["/gpfs21/scratch/cqs/shengq1/report/rawdata/s6_sequence.txt"],
  },
  tcga => {
    class            => "TCGA::RNAseq",
    perform          => 1,
    target_dir       => "${target_dir}",
    option           => "",
    source_ref       => "files",
    picard_dir       => $picard_dir,
    bedtools_dir     => $bedtools_dir,
    tcga_bin_dir     => $tcga_bin_dir,
    mapslice_version => "2",
    rename_bam       => 1,
    sh_direct        => 1,
    pbs              => {
      "email"    => $email,
      "nodes"    => "1:ppn=6",
      "walltime" => "72",
      "mem"      => "30gb"
    },
  },
};

performConfig($config);

1;
