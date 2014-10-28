#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201409_smallRNA_3018_1_soybean_chicken/");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $email = "quanhu.sheng\@vanderbilt.edu";

my $samtools           = "/home/shengq1/local/bin/samtools/samtools";
my $bowtie1_option_1mm = "-a -m 100 --best --strata -v 1 -p 8";

my $count_option = "-s --gtf_key gene";

my $bowtie1_option_pm = "-a -m 100 --best --strata -v 0 -p 8";

my $target_dir = $root;
my $config     = {
  general => { "task_name" => "3018_soybean_chicken" },
  files   => {
    "3018-KCV-1-1" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201409_smallRNA_3018/human/identical/result/3018-KCV-1-1_clipped_identical.fastq.gz"],
    "3018-KCV-1-2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201409_smallRNA_3018/human/identical/result/3018-KCV-1-2_clipped_identical.fastq.gz"]
  },
  count_files => {
    "3018-KCV-1-1" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201409_smallRNA_3018/human/identical/result/3018-KCV-1-1_clipped_identical.fastq.dupcount"],
    "3018-KCV-1-2" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD_Vickers/201409_smallRNA_3018/human/identical/result/3018-KCV-1-2_clipped_identical.fastq.dupcount"]
  },
  soybean_bowtie1_pm => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/soybean_bowtie1_pm",
    option        => $bowtie1_option_pm,
    source_ref    => "files",
    bowtie1_index => "/scratch/cqs/shengq1/references/soybean/bowtie_index_1.1.0/gma_ref_V1.1",
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  soybean_bowtie1_pm_count => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => "${target_dir}/soybean_bowtie1_pm_count",
    option          => $count_option . " --unmapped_fastq",
    source_ref      => "soybean_bowtie1_pm",
    fastq_files_ref => "files",
    seqcount_ref    => "count_files",
    cqs_tools       => $cqstools,
    gff_file        => "/scratch/cqs/shengq1/references/soybean/ref_V1.1_top_level.gff3",
    samtools        => $samtools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  chicken_bowtie1 => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_dir}/chicken_bowtie1_pm",
    option        => $bowtie1_option_pm,
    source_ref    => ["soybean_bowtie1_pm_count", "fastq.gz"],
    bowtie1_index => "/scratch/cqs/shengq1/references/chicken/bowtie_index_1.1.0/gga_ref_Gallus_gallus-4.0",
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  chicken_bowtie1_pm_count => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => "${target_dir}/chicken_bowtie1_pm_count",
    option          => $count_option,
    source_ref      => "chicken_bowtie1",
    fastq_files_ref => "files",
    seqcount_ref    => "count_files",
    cqs_tools       => $cqstools,
    gff_file        => "/scratch/cqs/shengq1/references/soybean/ref_V1.1_top_level.gff3",
    samtools        => $samtools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
};

performConfig($config);

1;

