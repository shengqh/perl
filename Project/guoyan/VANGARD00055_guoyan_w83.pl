#!/usr/bin/perl
use strict;
use warnings;

use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;

my $root     = create_directory_or_die("/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2");
my $cqstools = "/home/shengq1/cqstools/CQS.Tools.exe";

my $target_dir = $root;

my $email     = "quanhu.sheng\@vanderbilt.edu";
my $task_name = "VANGARD00055";

my $samtools = "/home/shengq1/local/bin/samtools/samtools";

my $bowtie1_option_pm       = "-a -m 100 --best --strata -v 0 -l 12 -p 8";
my $bowtie1_option_1mm      = "-a -m 100 --best --strata -v 1 -l 12 -p 8";
my $bowtie1_option_1mm_trim = "-a -m 100 --best --strata -v 1 -l 12 -p 8 --trim3 3";
my $bowtie1_option_3mm      = "-a -m 100 --best --strata -v 3 -l 12 -p 8";

my $mirnacount_option          = "-s";                                                    #ignore score
my $trnacount_option           = "--length --sequence";
my $mirna_overlap_count_option = "-s --gtf_key miRNA";
my $mirna_fasta                = "/data/cqs/shengq1/reference/miRBase20/mature.dna.fa";

my $w83_bowtie1_index = "/scratch/cqs/shengq1/references/gingivalis_W83/bowtie_1.0.0_index/Gingivalis_W83";
my $w83_gtf_index     = "/scratch/cqs/shengq1/references/gingivalis_W83/Gingivalis_W83.gtf";

my $target_w83_dir = create_directory_or_die( $target_dir . "/w83" );

my $w83config = {
  general  => { "task_name" => "w83", },
  bamfiles => {
    "2572-KCV-1-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-19/2572-KCV-1-19.bam"],
    "2572-KCV-1-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-20/2572-KCV-1-20.bam"],
    "2572-KCV-1-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-21/2572-KCV-1-21.bam"],
    "2572-KCV-1-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-22/2572-KCV-1-22.bam"],
    "2572-KCV-1-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-23/2572-KCV-1-23.bam"],
    "2572-KCV-1-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-24/2572-KCV-1-24.bam"],
    "2572-KCV-1-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-25/2572-KCV-1-25.bam"],
    "2572-KCV-1-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-26/2572-KCV-1-26.bam"],
    "2572-KCV-1-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-27/2572-KCV-1-27.bam"],
    "2572-KCV-1-28" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-28/2572-KCV-1-28.bam"],
    "2572-KCV-1-29" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-29/2572-KCV-1-29.bam"],
    "2572-KCV-1-30" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/human/topN_bowtie1_genome_cutadapt_1mm/result/2572-KCV-1-30/2572-KCV-1-30.bam"],
  },
  countfiles => {
    "2572-KCV-1-19" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-19_clipped_identical.dupcount"],
    "2572-KCV-1-20" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-20_clipped_identical.dupcount"],
    "2572-KCV-1-21" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-21_clipped_identical.dupcount"],
    "2572-KCV-1-22" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-22_clipped_identical.dupcount"],
    "2572-KCV-1-23" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-23_clipped_identical.dupcount"],
    "2572-KCV-1-24" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-24_clipped_identical.dupcount"],
    "2572-KCV-1-25" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-25_clipped_identical.dupcount"],
    "2572-KCV-1-26" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-26_clipped_identical.dupcount"],
    "2572-KCV-1-27" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-27_clipped_identical.dupcount"],
    "2572-KCV-1-28" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-28_clipped_identical.dupcount"],
    "2572-KCV-1-29" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-29_clipped_identical.dupcount"],
    "2572-KCV-1-30" => ["/gpfs21/scratch/cqs/shengq1/vangard/VANGARD00055_guoyan_mirna_v2/identical/result/2572-KCV-1-30_clipped_identical.dupcount"],
  },
  bam2fastq => {
    class               => "Bam2Fastq",
    perform             => 1,
    target_dir          => "${target_w83_dir}/bam2fastq",
    option              => "",
    source_ref          => "bamfiles",
    cqstools            => $cqstools,
    ispaired            => 0,
    unmapped_only       => 1,
    sort_before_convert => 0,
    sort_thread         => 12,
    unzipped            => 1,
    sh_direct           => 1,
    pbs                 => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  bowtie1_genome_cutadapt_topN_1mm => {
    class         => "Bowtie1",
    perform       => 1,
    target_dir    => "${target_w83_dir}/topN_bowtie1_genome_cutadapt_1mm",
    option        => $bowtie1_option_1mm,
    source_ref    => "bam2fastq",
    bowtie1_index => $w83_bowtie1_index,
    samonly       => 0,
    sh_direct     => 1,
    pbs           => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  count_1mm => {
    class           => "CQSMappedCount",
    perform         => 1,
    target_dir      => "${target_w83_dir}/bowtie1_genome_cutadapt_1mm_count",
    option          => "",
    source_ref      => "bowtie1_genome_cutadapt_topN_1mm",
    fastq_files_ref => "bam2fastq",
    seqcount_ref    => "countfiles",
    cqs_tools       => $cqstools,
    gff_file        => $w83_gtf_index,
    samtools        => $samtools,
    sh_direct       => 1,
    pbs             => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },
  tRNA_1mm_table => {
    class      => "CQSMappedTable",
    perform    => 1,
    target_dir => "${target_w83_dir}/summary",
    option     => "",
    source_ref => [ "count_1mm", ".xml" ],
    cqs_tools  => $cqstools,
    prefix     => "gingivalis_W83_1mm_",
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "10gb"
    },
  },
  overall => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${target_w83_dir}/overall",
    option     => "",
    source     => { individual => [ "bam2fastq", "bowtie1_genome_cutadapt_topN_1mm", "count_1mm" ], },
    sh_direct  => 1,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  }
};

performConfig($w83config);

1;
