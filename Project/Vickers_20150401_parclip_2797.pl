#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::PerformSmallRNA;
use CQS::FileUtils;

my $root    = "/scratch/cqs/shengq1/vickers/20150401_parclip_2797_rat_mouse_human/";
my $rat_def = {

  #General options
  task_name  => "Parclip2797",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => $root . "rn5",
  max_thread => 8,
  cluster    => "slurm",

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,
  adapter        => "TGGAATTCTCGGGTGCCAAGG",

  samtools => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  #Data
  files => {
    "RPI40" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI40_Ago2INS1Huh7.fastq.gz"],
    "RPI41" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI41_Ago3INS1Huh7.fastq.gz"],
    "RPI42" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI42_Ago2INS1HCEAC.fastq.gz"],
    "RPI43" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI43_Ago3INS1HCEAC.fastq.gz"],
  },
};

my $rn5config = performSmallRNA_rn5( $rat_def, 1 );

my $mouse_def = {

  #General options
  task_name  => "Parclip2797",
  email      => "quanhu.sheng\@vanderbilt.edu",
  target_dir => $root . "mm10",
  max_thread => 8,
  cluster    => "slurm",

  #Default software parameter (don't change it except you really know it)
  fastq_remove_N => 0,
  adapter        => "TGGAATTCTCGGGTGCCAAGG",

  samtools => "/scratch/cqs/shengq1/local/bin/samtools",
  cqstools => "/home/shengq1/cqstools/CQS.Tools.exe",

  #Data
  files => {
    "RPI47" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI47_Ago2MIN6Huh7.fastq.gz"],
    "RPI48" => ["/gpfs21/scratch/cqs/shengq1/vickers/data/2797_demultiplexing/2797-KCV-1_RPI48_Ago3MIN6Huh7.fastq.gz"],
  },
};

my $mm10config = performSmallRNA_mm10( $mouse_def, 1 );

my $hg19genome   = hg19_genome();
my $samtools     = "/scratch/cqs/shengq1/local/bin/samtools";
my $cqstools     = "/home/shengq1/cqstools/CQS.Tools.exe";
my $utr3_db      = "/data/cqs/shengq1/reference/utr3/20140612_ucsc_hg19_3UTR.txt";
my $email        = "quanhu.sheng\@vanderbilt.edu";
my $bowtie1_pm   = "-a -m 100 --best --strata -v 0 -p 8";
my $hg19target   = create_directory_or_die("${root}hg19");
my $fasta_file   = "/scratch/cqs/shengq1/references/hg19_16569_M/hg19_16569_M.fa";
my $refgene_file = "/scratch/cqs/shengq1/references/hg19_16569_M/hg19_refgene.tsv";

my $hg19config = {
  general           => { "task_name" => "Parclip2797", },
  bowtie1_genome_pm => {
    class             => "Bowtie1",
    perform           => 1,
    target_dir        => "${hg19target}/bowtie1_genome_pm",
    option            => $bowtie1_pm,
    source_config_ref => [ $rn5config, "identical", ".fastq.gz\$", $mm10config, "identical", ".fastq.gz\$" ],
    bowtie1_index     => $hg19genome->{bowtie1_index},
    samonly           => 0,
    sh_direct         => 1,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  },

  t2c_finder => {
    class             => "CQS::ParclipT2CFinder",
    perform           => 1,
    target_dir        => "${hg19target}/t2c_finder",
    option            => "-p 0.05 -e 0.013",
    source_config_ref => [ $rn5config, "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml\$", $mm10config, "bowtie1_genome_1mm_NTA_smallRNA_count", ".mapped.xml\$" ],
    cqs_tools         => $cqstools,
    sh_direct         => 1,
    pbs               => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  utr3_count => {
    class                  => "CQS::SmallRNACount",
    perform                => 1,
    target_dir             => "${hg19target}/utr3_count",
    option                 => "-m 0",
    source_ref             => [ "bowtie1_genome_pm", ".bam\$" ],
    fastq_files_config_ref => [ $rn5config, "identical", ".fastq.gz\$", $mm10config, "identical", ".fastq.gz\$" ],
    seqcount_config_ref    => [ $rn5config, "identical", ".dupcount\$", $mm10config, "identical", ".dupcount\$" ],
    cqs_tools              => $cqstools,
    coordinate_file        => $utr3_db,
    samtools               => $samtools,
    sh_direct              => 1,
    pbs                    => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  mirna_utr3_target => {
    class        => "CQS::ParclipMirnaTarget",
    perform      => 1,
    target_dir   => "${hg19target}/mirna_utr3_target",
    option       => "",
    source_ref   => [ "t2c_finder", ".xml\$" ],
    target_ref   => [ "utr3_count", ".xml\$" ],
    fasta_file   => $fasta_file,
    refgene_file => $refgene_file,
    cqs_tools    => $cqstools,
    sh_direct    => 1,
    pbs          => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
  sequencetask => {
    class      => "CQS::SequenceTask",
    perform    => 1,
    target_dir => "${hg19target}/sequencetask",
    option     => "",
    source     => { one => [ "bowtie1_genome_pm", "t2c_finder", "utr3_count", "mirna_utr3_target" ], },
    sh_direct  => 0,
    pbs        => {
      "email"    => $email,
      "nodes"    => "1:ppn=1",
      "walltime" => "72",
      "mem"      => "20gb"
    },
  },
};

performConfig($hg19config);

1;

